#include "MemoryBudget.hpp"
#include <exception>

// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BayesOmics, a statistical genetics software package
// developed by Shouye Liu.
//
// Copyright (C) 2025 Shouye Liu <syliu.xue@foxmail.com>
//
// BayesOmics is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// BayesOmics is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with BayesOmics. If not, see <https://www.gnu.org/licenses/>.


#include "ParallelGibbs.hpp"
#include "ModelBayesCO.hpp"
#include "OmicsMixture.hpp"

void BayesCO::Intercept::sampleFromFC(VectorXd &wbcorr, const double &vare,const VectorXd &nGWAS){
    wbcorr = wbcorr.array() + value;
    double rhs = wbcorr.sum();
    double invLhs = 1/ ((double) nGWAS.mean());
    double muHat = invLhs * rhs;
    value = Normal::sample(muHat,invLhs * vare);
    wbcorr.array() = wbcorr.array() - value;
}

void BayesCO::InterceptEQTL::sampleFromFC(vector<VectorXd> &wAcorr, const map<string, vector<int> > &genePheIdxMap, const VectorXd &varEps, const vector<VectorDat> &neQTL){
    for (unsigned i = 0; i < numGenes; i++){
        wAcorr[i] = wAcorr[i].array() + values[i];
        double rhs = wAcorr[i].sum();
        double invLhs = 1/((double) neQTL[i].values.mean());
        double muwHat = invLhs * rhs;
        values[i] = Normal::sample(muwHat,invLhs * varEps[i]);
        wAcorr[i].array() = wAcorr[i].array() - values[i];
    }
}



// EIEO model

void BayesCO::SnpEffects::sampleFromFCEIEO(const Data &data, VectorXd &ycorr,
    vector<VectorXd> &wcorr, SnpEffectVec &beta, EQTLJointVec &alpha,
    SnpEffectVec &betaLatent, EQTLJointVec &alphaLatent, DeltaMat &components,
    const VectorXd &gamma, const VectorXd &piNonEqtl, const MatrixXd &piEqtl,
    double sigmaNonEqtl, const vector<MatrixXd> &covarianceInverse,
    const VectorXd &varEps, double vare) {

    // CO-EIEO retains sequential SNP Gibbs updates. Cache only immutable
    // molecular column norms; parallel preparation never draws random numbers.
    if(geneGenotypeNorms.size()!=alpha.size()) {
        geneGenotypeNorms.resize(alpha.size());
        vector<string> orderedGenes(alpha.size());
        for(const auto &item:data.geneID2IdxMap)orderedGenes.at(item.second)=item.first;
        for(unsigned g=0;g<alpha.size();++g)geneGenotypeNorms[g].resize(alpha[g]->values.size());
        int normWorkers=ParallelGibbs::workers(alpha.size());
        if(MemoryBudget::enabled()) {
            const auto perWorker=MemoryBudget::bytes(data.geneGenotype().rows(),1);
            normWorkers=std::min<int>(normWorkers,std::max<uint64_t>(1,MemoryBudget::status().available/4/std::max<uint64_t>(1,perWorker)));
        }
        vector<std::exception_ptr> normErrors(alpha.size());
        #pragma omp parallel for schedule(dynamic) num_threads(normWorkers) if(!data.geneGenotype().streamed() && normWorkers>1)
        for(unsigned g=0;g<alpha.size();++g) {
          try {
            const auto &individuals=data.genePheIdxMap.at(orderedGenes[g]);
            for(const auto &item:alpha[g]->name2index) {
                const VectorXd column=data.geneGenotype().col(data.cisSnpID2IdxMap.at(item.first))(individuals);
                geneGenotypeNorms[g][item.second]=ParallelGibbs::dot(column,column);
            }
          } catch(...) {normErrors[g]=std::current_exception();}
        }
        for(const auto &error:normErrors)if(error){geneGenotypeNorms.clear();std::rethrow_exception(error);}
    }
    if(genotypeNorms.size()!=data.genotype().cols())ParallelGibbs::columnNorms(data.genotype(),genotypeNorms);
    const unsigned K=gamma.size(), G=data.numKeptGenes;
    const OmicsMixture::Draws draws(data);
    numNonZeros=numNonZerosNonEqtl=numNonNullEqtl=numNonNullSnpTot=0;
    numNonZerosEqtlVec.setZero(); ssqNonEqtl=ssqBetaEqtl=ssqAlphaEqtl=0;
    ssqEqtlMat.assign(G,Matrix2d::Zero());
    numNonZerosEqtlVecAcrossGenesPostIW.assign(G,0);
    componentCountsNonEqtl.setZero(K); componentCountsEqtl.setZero(2,K);
    for (auto *p : components) p->values.setZero();
    for (unsigned j=0;j<data.numIncdSnps;++j) {
        panelScores.begin(j,data.genotype(),ycorr,genotypeNorms);
        const string &id=data.snpEffectNames[j];
        std::size_t draw=draws.offset[j];
        const auto x=data.genotype().col(j);
        const double xx=genotypeNorms[j];
        auto genes=data.gwasSnpID2geneIDMap.find(id);
        if (genes==data.gwasSnpID2geneIDMap.end()) {
            const double old=beta[0]->getValue(id), score=((panelScores.enabled()?panelScores.score(j):ParallelGibbs::dot(x,ycorr))+xx*old)/vare;
            auto conditional=OmicsMixture::eieo(xx/vare,score,1/sigmaNonEqtl,0,0,gamma,piNonEqtl);
            const unsigned k=OmicsMixture::choose(conditional.probability,draws.uniform[draw]);
            const double effect=k ? sqrt(gamma[k])*(conditional.mean[k]+sqrt(conditional.variance[k])*draws.normal[draw]) : 0;
            beta[0]->setValue(id,effect); ParallelGibbs::axpy(x,old-effect,ycorr);panelScores.update(j,old-effect);
            ++componentCountsNonEqtl[k]; components[0]->values(j,0)=k;
            if (k) {ssqNonEqtl+=effect*effect/gamma[k]; ++numNonZerosNonEqtl; ++numNonZeros; ++numNonNullSnpTot;}
            continue;
        }
        bool activeBeta=false, activeAlpha=false;
        for (const string &gene : genes->second) {
            const unsigned g=data.geneID2IdxMap.at(gene);
            const int ci=data.cisSnpID2IdxMap.at(id);
            const VectorXd z=data.geneGenotype().col(ci)(data.genePheIdxMap.at(gene));
            const double zz=geneGenotypeNorms[g][alpha[g]->name2index.at(id)];
            const Vector2d old(beta[g+1]->getValue(id),alpha[g]->getValue(id));
            Vector2d latent(betaLatent[g+1]->getValue(id),alphaLatent[g]->getValue(id));
            const Vector2d score((ParallelGibbs::dot(x,ycorr)+xx*old[0])/vare,(ParallelGibbs::dot(z,wcorr[g])+zz*old[1])/varEps[g]);
            const Vector2d precision(xx/vare,zz/varEps[g]);
            for (unsigned t=0;t<2;++t) {
                const auto &inv=covarianceInverse[g];
                auto conditional=OmicsMixture::eieo(precision[t],score[t],inv(t,t),inv(t,1-t),latent[1-t],gamma,piEqtl.row(t));
                const unsigned k=OmicsMixture::choose(conditional.probability,draws.uniform[draw+t]);
                latent[t]=(conditional.mean[k]+sqrt(conditional.variance[k])*draws.normal[draw+t]);
                const double effect=sqrt(gamma[k])*latent[t];
                ++componentCountsEqtl(t,k); components[t]->values(j,g+1)=k;
                if (t==0) {
                    beta[g+1]->setValue(id,effect); betaLatent[g+1]->setValue(id,latent[t]);
                    ParallelGibbs::axpy(x,old[t]-effect,ycorr);panelScores.update(j,old[t]-effect);
                    if(k) {activeBeta=true; ++numNonZerosEqtlVec[t]; ssqBetaEqtl+=latent[t]*latent[t];}
                } else {
                    alpha[g]->setValue(id,effect); alphaLatent[g]->setValue(id,latent[t]);
                    ParallelGibbs::axpy(z,old[t]-effect,wcorr[g]);
                    if(k) {activeAlpha=true; ++numNonZerosEqtlVec[t]; ssqAlphaEqtl+=latent[t]*latent[t];}
                }
            }
            draw+=2;
            // All latent pairs inform the covariance, including zero components.
            ssqEqtlMat[g].noalias()+=latent*latent.transpose();
            ++numNonZerosEqtlVecAcrossGenesPostIW[g];
        }
        if(activeBeta) {++numNonZeros; ++numNonNullSnpTot;}
        if(activeAlpha) ++numNonNullEqtl;
    }
    rebuildPredictions(data,beta,alpha);
}


void BayesCO::GenotypicVar::compute(VectorXd &betaTotal,const MatrixXd &Z){
    value = Gadget::calcVariance(Z * betaTotal); 
}

void BayesCO::GenotypicVarGene::compute(VectorXd &geneEffects,  MatrixXd &eQTLJointMat, const MatrixXd &ZGene, const map<string, vector<int> > &genePheIdxMap, const vector<string> &geneEffectNames,const map<int,vector<int>> &gene2cisSnpMap){
        vector<int> snpIdxInGene;
        vector<int> genePheIdxPerGene;
        vector<int> gene2cisSnpIdxPerGene;
        value = 0;
        int numGenes = gene2cisSnpMap.size();
        for(unsigned i = 0; i < numGenes; i++ ){
            gene2cisSnpIdxPerGene = gene2cisSnpMap.at(i);
            genePheIdxPerGene = genePheIdxMap.at(geneEffectNames[i]);
            VectorXd ZE =  ZGene(genePheIdxPerGene,gene2cisSnpIdxPerGene)* eQTLJointMat(gene2cisSnpIdxPerGene,i);
            value = value +  Gadget::calcVariance(ZE * geneEffects[i]);
        }                           
}

void BayesCO::GenotypicVarGeneCis::compute( MatrixXd &eQTLJointMat, const MatrixXd &ZGene, 
                const map<string, vector<int> > &genePheIdxMap, const vector<string> &geneEffectNames, const map<int,vector<int>> &gene2cisSnpMap){
        vector<int> snpIdxInGene;
        vector<int> genePheIdxPerGene;
        vector<int> gene2cisSnpIdxPerGene;
        int numGenes = gene2cisSnpMap.size();
        for(unsigned i = 0; i < numGenes; i++ ){
            gene2cisSnpIdxPerGene = gene2cisSnpMap.at(i);
            genePheIdxPerGene = genePheIdxMap.at(geneEffectNames[i]);
            // values[i] = Gadget::calcVariance( ZGene(genePheIdxPerGene,gene2cisSnpIdxPerGene) * eQTLJointMat(snpIdxInGene,i) ) ;
            values[i] = Gadget::calcVariance( ZGene(genePheIdxPerGene,gene2cisSnpIdxPerGene) * eQTLJointMat(gene2cisSnpIdxPerGene,i) ) ;

        }                        
}

void BayesCO::ResidualVar::sampleFromFC( const VectorXd &ycorr, const VectorXd &nGWAS){
    // double sse = (y - X * snpEffect ).dot((y - X * snpEffect ));
    double sse = ycorr.dot(ycorr);
    double dfTilde = df + nGWAS.mean();
    double scaleTilde = sse + df*scale;
    value = InvChiSq::sample(dfTilde, scaleTilde);
    
}

void BayesCO::ResidualVar::sampleFromFC( const VectorXd &y,const MatrixXd X, const VectorXd &snpEffect, const VectorXd &nGWAS){
    double sse = (y - X * snpEffect ).dot((y - X * snpEffect ));
    double dfTilde = df + nGWAS.mean();
    double scaleTilde = sse + df*scale;
    value = InvChiSq::sample(dfTilde, scaleTilde);
}

void BayesCO::ResidualVareEQTL::sampleFromFC(const vector<VectorXd> &wAcorr, const vector<VectorDat> &neQTL){
    int numGenes = wAcorr.size();
    for (int i = 0; i < numGenes; i++){
        double sseEps = wAcorr[i].dot(wAcorr[i]);
        double dfTilde = (*this)[i]->df  + (double) neQTL[i].values.mean();
        double scaleTilde = sseEps + (*this)[i]->df * (*this)[i]->scale;
        values[i] = Stat::InvChiSq::sample(dfTilde, scaleTilde);
    }

}

void BayesCO::ResidualVareEQTL::sampleFromFC( const vector<VectorXd> &w, const MatrixXd ZGene, EQTLJointVec &eQTLJointVec, const vector<VectorDat> &neQTL,
                                            const map<string, vector<int> > &genePheIdxMap, const vector<string> &geneEffectNames,
                                            const map<int,vector<int>> &gene2cisSnpMap){
    int numGenes = w.size();
    vector<int> genePheIdxPerGene;
    vector<int> gene2cisSnpIdxPerGene;
    for (int i = 0; i < numGenes; i++){
        gene2cisSnpIdxPerGene = gene2cisSnpMap.at(i);
        genePheIdxPerGene = genePheIdxMap.at(geneEffectNames[i]);
        VectorXd sseVec = w[i];
        sseVec = sseVec - ZGene(genePheIdxPerGene,gene2cisSnpIdxPerGene) * eQTLJointVec[i]->values;
        double sseEps = sseVec.dot(sseVec);
        double dfTilde = (*this)[i]->df  + (double) neQTL[i].values.mean();
        double scaleTilde = sseEps + (*this)[i]->df * (*this)[i]->scale;
        values[i] = Stat::InvChiSq::sample(dfTilde, scaleTilde);
    }
}

void BayesCO::SigmaSqBetaNonEqtl::sampleFromFC(const double snpEffSumSq, const unsigned numSnpEff){
    double dfTilde = df + numSnpEff;
    double scaleTilde = snpEffSumSq + df * scale;
    value = InvChiSq::sample(dfTilde, scaleTilde);
}

void BayesCO::SigmaSqAlphaVec::sampleFromFC(const MatrixXd &deltaMat, const MatrixXd &eQTLJointMat, const map<int,vector<int>> &gene2cisSnpMap){
    VectorXd snpEffSumSq(geneNames.size());
    VectorXd numSnpEff(geneNames.size());
    vector<int> eQTLCommonCis;
    // construct snpEffectSumSq and numSnpEff by using eQTLJointMat 
    for(int i =0; i < geneNames.size(); ++i){
        eQTLCommonCis.clear();
        vector<int> snpIdxInGeneSet = gene2cisSnpMap.at(i);
        for(int j = 0; j < snpIdxInGeneSet.size(); j++){
            if(eQTLJointMat(snpIdxInGeneSet[j],i) !=0 )  eQTLCommonCis.push_back(snpIdxInGeneSet[j]);
        }
        snpEffSumSq[i] = eQTLJointMat(eQTLCommonCis,i).sum();
        numSnpEff[i]  = eQTLCommonCis.size();
        (*this)[i]->sampleFromFC(snpEffSumSq[i], numSnpEff[i]);
        values[i] = (*this)[i]->value;
    }
}

void BayesCO::SigmaSqAlphaVec::sampleFromPrior(){
    for(unsigned i = 0; i < geneNames.size();++i){
        values[i] = InvChiSq::sample((*this)[i]->df, (*this)[i]->scale);
    }
}

void BayesCO::SigmaSqMat::setPrior(const double &sigmaSqBetaEqtl, const VectorXd &sigmaSqAlphaVec){
    Matrix2d varcov, geneCor;
    geneCor << 1, 0, 0, 1;
    for(int i =0; i  < numGenes; ++i){
        varcov << sigmaSqBetaEqtl, 0, 0, sigmaSqAlphaVec[i];
        varcovPriors[i] = varcov * 0.5;
        sigmaSqMats[i] = varcov;
        (*this)[i]->values=varcov;
        sigmaSqInvMats[i] = sigmaSqMats[i].inverse();
        sigmaSqDetLogVec[i] = log( sigmaSqMats[i].determinant());
    }

    if (numGenes != 0){
        sigmaSqBetaEqtlPM = sigmaSqBetaEqtl;
        sigmaSqAlphaAll = sigmaSqAlphaVec.mean();
        scaleBetaEqtl = (nub-2)/nub * sigmaSqBetaEqtl;
        scaleAlphaAll = ((nua-2)/nua * sigmaSqAlphaVec.array()).mean();
    }
}

void BayesCO::SigmaSqMat::sampleFromFCInvWishartGeneral(const double &ssqBetaEqtl, const double &ssqAlphaEqtl, vector<Matrix2d> ssqEqtlMat,
         const VectorXd &numNonZerosEqtlVec, const vector<unsigned> numNonZerosEqtlVecAcrossGenesPostIW,const bool messageBool){
    sigmaSqBetaEqtlPM = InvChiSq::sample(nub + numNonZerosEqtlVec(0), ssqBetaEqtl + nub * scaleBetaEqtl );
    sigmaSqAlphaAll  = InvChiSq::sample(nua + numNonZerosEqtlVec(1), ssqAlphaEqtl + nua * scaleAlphaAll );
    Matrix2d varcovPriorsChr;
    double dfIW = 0.0;
    varcovPriorsChr << sigmaSqBetaEqtlPM * 0.5 , 0.0, 0.0, sigmaSqAlphaAll * 0.5 ;
    Matrix2d   sampleEigen;
    for(int i =0; i < numGenes; ++i){

        MatrixXd effArmEigen = ssqEqtlMat[i] + varcovPriorsChr;
        dfIW = 4.0 + numNonZerosEqtlVecAcrossGenesPostIW[i];
        arma::dmat effArma = arma::dmat(effArmEigen.data(),effArmEigen.rows(),effArmEigen.cols(),false,false);
        // arma::dmat psiParam = arma::dmat(varcovPriorsChr.data(), varcovPriorsChr.rows(), varcovPriorsChr.cols(),false, false);
        arma::dmat psiParam = effArma ; 
        /////////////////////////////////////////////////////////////////////
        // if(!Gadget::checkScaleMatrix(psiParam,dfIW,false)) {
        //     continue;
        // } // psiParam is not symmetric/hermitian positive definite;
        arma::dmat sample = arma::iwishrnd(psiParam, dfIW);
        // arma::dmat sample = stats::rinvwish(psiParam,dfIW);
        MatrixXd   sampleEigen = Eigen::Map<Eigen::MatrixXd>(sample.memptr(),sample.n_rows, sample.n_cols);
        /////////////////////////////////////////////////////////////////////
        sigmaSqMats[i] = sampleEigen;
        (*this)[i]->values=sampleEigen;
        sigmaSqInvMats[i] = sigmaSqMats[i].inverse();
        sigmaSqDetLogVec[i] = logf( sigmaSqMats[i].determinant());
    }
}

void BayesCO::GeneEffects::sampleFromeFC(Data data, int iter , VectorXd &betaTotal, MatrixXd &eQTLJointMat,
                            const map<int,vector<int> > &gene2gwasSnpMap,const map<int,vector<int>> &gene2cisSnpMap, 
                            double &sigmaSqTheta, double &vareMed, const double piTheta, VectorXd deltaGene ){
    double sample, oldSample,varRes,rhs,uhat,invLhs,muHat;
    double logDelta0, logDelta1, probDelta1;
    double sigmaSqThetaInv, sigmaSqThetaLog;
    double muBeta;
    vector<int> eQTLCommon; // gwas-snp
    vector<int> eQTLCommonCis; // cis-snp
    deltaGene.setZero();
    VectorXd betaCorr, betaHat;
    MatrixXd H, HPH; 
    set<int> eQTLUniqSet; 
    set<int>::iterator iterSet;
    for(int i =0; i < numGenes; ++i){
        vector<int> snpIdxInGeneSet = gene2gwasSnpMap.at(i);
        vector<int> gene2cisSnpIdxPerGene = gene2cisSnpMap.at(i);
        for(int j = 0; j < snpIdxInGeneSet.size(); j++){
            int gwasIdx = snpIdxInGeneSet[j];
            int eqtlIdx = gene2cisSnpIdxPerGene[j];
            // iterSet = eQTLUniqSet.find(gwasIdx);
            if(eQTLUniqSet.find(gwasIdx) != eQTLUniqSet.end()){
                continue;
            }
            if (betaTotal[gwasIdx] !=0 && eQTLJointMat(eqtlIdx,i) != 0 ){
                eQTLCommon.push_back(gwasIdx);
                eQTLCommonCis.push_back(eqtlIdx);
                eQTLUniqSet.insert(gwasIdx); 
            }
        }
    }
    nnGene = 0;
    ssqGene = 0;
    // sigmaSqTheta = 0.0484;
    // piTheta = 0.000829633134322635;
    sigmaSqThetaInv = 1/sigmaSqTheta;
    sigmaSqThetaLog = logf(sigmaSqTheta);

    betaCorr = betaTotal(eQTLCommon);
    // cout << "betacorr: " << betaCorr.head(6) << betaCorr.tail(6) << endl;
    H = eQTLJointMat(eQTLCommonCis,Eigen::indexing::all);
    // cout << "H: " << H << endl;
    HPH = (H.transpose() * H).diagonal();
    // cout << "HPH: " << HPH << endl;
    muBeta = betaCorr.mean();
    betaCorr.array() = betaCorr.array() - muBeta - (H * values).array();
    /// Start mcmc process
    for(unsigned iteri =0; iteri < 1; iteri++){
        betaCorr.array() += muBeta;
        rhs = betaCorr.sum();
        invLhs = 1/(double) eQTLCommon.size();
        muHat = invLhs*rhs;
        muBeta = Normal::sample(muHat,invLhs*vareMed);
        betaCorr.array() -= muBeta;
        betaHat.setZero(eQTLCommon.size()); // betaHat = 0;

        for (unsigned i = 0; i < numGenes; i++) {
            double logPi = logf(piTheta);
            double logPiComp = logf(1-piTheta);
            oldSample = values[i];
            rhs = H.col(i).dot(betaCorr);
            rhs = rhs + HPH(i) * oldSample;
            rhs = rhs/vareMed;
            invLhs = 1.0/(HPH(i)/vareMed + sigmaSqThetaInv);
            uhat = invLhs*rhs;
            logDelta0 = logPiComp;
            logDelta1 = 0.5f*(logf(invLhs) - sigmaSqThetaLog + uhat*rhs) + logPi;
            probDelta1   = 1.0/(1.0 + expf(logDelta0-logDelta1));
            if (Stat::ranf() < probDelta1) {
            // if (Stat::ranf() < 1.0) {
                sample = uhat + Stat::snorm()* sqrt(invLhs);
                sample = Normal::sample(uhat,invLhs);
                values[i] = sample;
                ssqGene += sample * sample;
                betaCorr = betaCorr + H.col(i)*(oldSample - values[i]);
                betaHat = betaHat + H.col(i) * values[i];
                deltaGene[i] = 1;
                nnGene = nnGene + 1;
            } else {
                if (oldSample) {
                    betaCorr = betaCorr + H.col(i)*oldSample;
                }
                values[i] = 0;
            }
        } // End of gene loop
    } // End of mcmc sampling loop
    vareMed = Gadget::calcVariance(betaCorr);
    propMed = Gadget::calcVariance(betaHat)/(Gadget::calcVariance(betaHat) + vareMed);
}


void BayesCO::SigmaSqTheta::sampleFromFC(const double snpEffSumSq, const unsigned numSnpEff){
    double dfTilde = df + numSnpEff;
    double scaleTilde = snpEffSumSq + df * scale;
    value = InvChiSq::sample(dfTilde, scaleTilde);
}

void BayesCO::setStartVal(const Data &data){
    if(data.numKeptGenes)
        sigmaSqMats.setPrior(sigmaSqBetaEqtl.value, sigmaSqAlphaVec.values);
}

void BayesCO::sampleUnknowns(){
    /////////////////////////////////////
    ///// debug real simu situation /////
    /////////////////////////////////////
    unsigned outi = 0;

    static int iter = 0;
    unsigned cnt=0;
    intercept.sampleFromFC(wbcorr,vare.value,data.n);
    interceptEqtl.sampleFromFC(wAcorr,data.genePheIdxMap,varEps.values,data.neQTLVec);
    // genic effects sampling
    

    if (mcmcType == "EIEO") {
        Vector2d gamma(0,1), piNon(1-piEffNonEqtl.value,piEffNonEqtl.value);
        MatrixXd piGene(2,2);
        piGene << 1-piEffEieo1.value,piEffEieo1.value,1-piEffEieo2.value,piEffEieo2.value;
        snpEffects.sampleFromFCEIEO(data,wbcorr,wAcorr,snpEffectVec,eQTLJointVec,
            snpEffectVecLatent,eQTLJointVecLatent,deltaTrait,gamma,piNon,piGene,
            sigmaSqBetaNonEqtl.value,sigmaSqMats.sigmaSqInvMats,varEps.values,vare.value);
        if (estimatePi) {
            piEffEieo1.sampleFromFC(data.numEqtlOverlap,snpEffects.numNonZerosEqtlVec[0]);
            piEffEieo2.sampleFromFC(data.numEqtlOverlap,snpEffects.numNonZerosEqtlVec[1]);
            piEffNonEqtl.sampleFromFC(data.numNonEqtl,snpEffects.numNonZerosNonEqtl);
        }
        geneEffectVec.sampleSecondary(data,snpEffects.betaTotal,eQTLJointVec,sigmaSqTheta.value,piTheta.value);
        if (geneEffectVec.secondaryAvailable) {
            sigmaSqTheta.sampleFromFC(geneEffectVec.values.squaredNorm(),geneEffectVec.nnGene);
            piTheta.sampleFromFC(geneEffectVec.numGenes,geneEffectVec.nnGene);
        }
        sigmaSqMats.sampleFromFCInvWishartGeneral(snpEffects.ssqBetaEqtl,snpEffects.ssqAlphaEqtl,
            snpEffects.ssqEqtlMat,snpEffects.numNonZerosEqtlVec,snpEffects.numNonZerosEqtlVecAcrossGenesPostIW,false);
    }

    // vare.sampleFromFC(data.y,data.Z,snpEffects.betaTotal, data.n);

    vareMedParameter.value=geneEffectVec.vareMed;
    thetaInterceptParameter.value=geneEffectVec.muRegression;
    if (!data.fixedResidual) vare.sampleFromFC(wbcorr,data.n);
    vareMean.value = vare.value;
    if(data.numKeptGenes != 0) {
        varEpsMean.value = varEps.values.mean();
        // varEps.sampleFromFC(data.genePheVec,data.ZGene,eQTLJointVec,data.neQTLVec, data.genePheIdxMap, data.geneEffectNames,data.gene2cisSnpMap);
        if (!data.fixedResidual) varEps.sampleFromFC(wAcorr,data.neQTLVec);
    }
    for (unsigned g=0;g<data.numKeptGenes;++g) {
        residualMats[g]->values=Vector2d(vare.value,varEps.values[g]).asDiagonal();
    }
    nnzBtw.value=snpEffects.numNonZerosNonEqtl;
    // summary
    nnsTot.getValue(snpEffects.numNonNullSnpTot);  // nnSnp
    if(data.numKeptGenes != 0){
        nnzGen.getValue(snpEffects.numNonZerosEqtlVec(0));
        nnsGen.getValue(snpEffects.numNonNullEqtl);
        nnGene.getValue(geneEffectVec.nnGene);
        if(nnGene.value != 0){
            nnsPG.value = nnsGen.value / nnGene.value; // Average nonZeroEqtl per nonZeroGene nnEqtlPerGene
        } else {
            nnsPG.value = 0;
        }
    } else {
        nnsPG.value = 0;
    }
    if(data.numNonEqtl != 0) {
        sigmaSqBetaNonEqtl.sampleFromFC(snpEffects.ssqNonEqtl, snpEffects.numNonZerosNonEqtl);
    } else {
        nnzBtw.value = 0;
        sigmaSqBetaNonEqtl.value = 0;
    }
    // new code
    varg.compute(snpEffects.ghat);
    vargGeneCis.compute(snpEffects.gwhatMap);
    vargGene.compute(data.numKeptInds,snpEffects.gwhatGwasMap,geneEffectVec.values);
    hsq.compute(varg.value, vare.value);

    if(data.numKeptGenes != 0){
        medHsq.value=vargGene.value/data.varPhenotypic;
        cisHsq.compute(vargGeneCis.values, data.varPhenotypiceQTL);

        sigmaSqAlpha.value = sigmaSqMats.sigmaSqAlphaAll;
        sigmaSqBetaEqtl.value = sigmaSqMats.sigmaSqBetaEqtlPM;
        cisHsqMean.value = cisHsq.values.mean();
    } else {
        medHsq.value = 0;  
        sigmaSqAlpha.value = 0;
        sigmaSqBetaEqtl.value = 0;  
        cisHsqMean.value = 0; 
    }
    ++iter;
}

// Derived predictions are rebuilt from the current coefficients, never accumulated across sweeps.
void BayesCO::SnpEffects::rebuildPredictions(const Data &data, const SnpEffectVec &beta, const EQTLJointVec &alpha) {
    betaTotal.setZero(data.numIncdSnps);
    ghat.setZero(data.numKeptInds);
    for (unsigned j=0;j<data.numIncdSnps;++j) {
        const string &id=data.snpEffectNames[j];
        auto linked=data.gwasSnpID2geneIDMap.find(id);
        if (linked==data.gwasSnpID2geneIDMap.end()) betaTotal[j]=beta[0]->getValue(id);
        else for (const auto &g : linked->second) betaTotal[j]+=beta[data.geneID2IdxMap.at(g)+1]->getValue(id);

    }
    ParallelGibbs::sparseProduct(data.genotype(),betaTotal,ghat);
    values=betaTotal;
    gwhatMap.clear(); gwhatGwasMap.clear();
    vector<VectorXd> localPredictions(data.numKeptGenes),globalPredictions(data.numKeptGenes);
    vector<std::exception_ptr> errors(data.numKeptGenes);
    #pragma omp parallel num_threads(ParallelGibbs::workers(data.numKeptGenes)) if(data.numKeptGenes>1 && !data.genotype().streamed() && !data.geneGenotype().streamed())
    {
    #pragma omp single
    predictionWorkerCount=omp_get_num_threads();
    #pragma omp for schedule(dynamic)
    for (unsigned g=0;g<data.numKeptGenes;++g) {
      try {
        const auto &ids=data.genePheIdxMap.at(data.geneEffectNames[g]);
        VectorXd local=VectorXd::Zero(ids.size()), global=VectorXd::Zero(data.numKeptInds);
        for (const auto &snp : data.gene2cisSnpIDMap.at(g)) {
            const double effect=alpha[g]->getValue(snp);
            if (effect==0.0) continue;
            const int ci=data.cisSnpID2IdxMap.at(snp), gi=data.snpInfoMap.at(snp)->index;
            const auto localColumn=data.geneGenotype().col(ci);
            for (unsigned i=0;i<ids.size();++i) local[i]+=localColumn[ids[i]]*effect;
            ParallelGibbs::axpy(data.genotype().col(gi),effect,global);
        }
        localPredictions[g]=std::move(local);globalPredictions[g]=std::move(global);
      } catch(...) {errors[g]=std::current_exception();}
    }
    }
    for(const auto &error:errors)if(error)std::rethrow_exception(error);
    for(unsigned g=0;g<data.numKeptGenes;++g) {
        gwhatMap.emplace(g,std::move(localPredictions[g]));
        gwhatGwasMap.emplace(g,std::move(globalPredictions[g]));
    }

}

// Independent secondary mixture regression of total cis beta on current alpha.
void BayesCO::GeneEffects::sampleSecondary(const Data &data, const VectorXd &betaTotal,
                                         const vector<ParamSet*> &alpha, double sigmaTheta, double piTheta) {
    std::set<int> activeSet;
    for(unsigned g=0;g<numGenes;++g)
        for(const auto &snp:data.gene2cisSnpIDMap.at(g))
            if(alpha[g]->getValue(snp)!=0)activeSet.insert(data.snpInfoMap.at(snp)->index);
    vector<int> active(activeSet.begin(),activeSet.end());
    deltaGene.setZero(); nnGene=0; ssqGene=0; secondaryAvailable=active.size()>=2;
    if (!secondaryAvailable) { values.setZero(); muRegression=0.0; return; }
    MemoryBudget::require(MemoryBudget::bytes(active.size(),numGenes+6),"secondary theta active-cis workspace");
    MatrixXd H=MatrixXd::Zero(active.size(),numGenes);VectorXd target(active.size());
    std::map<int,unsigned> row;for(unsigned j=0;j<active.size();++j){row.emplace(active[j],j);target[j]=betaTotal[active[j]];}
    for(unsigned g=0;g<numGenes;++g)for(const auto &snp:data.gene2cisSnpIDMap.at(g)){
        const double effect=alpha[g]->getValue(snp);if(effect!=0)H(row.at(data.snpInfoMap.at(snp)->index),g)=effect;
    }
    VectorXd residual=target-H*values;
    muRegression=Normal::sample(residual.mean(),vareMed/active.size());
    residual.array()-=muRegression;
    for (unsigned g=0;g<numGenes;++g) {
        const double old=values[g], xx=H.col(g).squaredNorm();
        const double rhs=(H.col(g).dot(residual)+xx*old)/vareMed;
        const double v=1.0/(xx/vareMed+1.0/sigmaTheta), mean=v*rhs;
        const double odds=0.5*(log(v/sigmaTheta)+mean*rhs)+log(piTheta)-log1p(-piTheta);
        const double prob=odds>=0 ? 1.0/(1.0+exp(-odds)) : exp(odds)/(1.0+exp(odds));
        values[g]=Stat::ranf()<prob ? Normal::sample(mean,v) : 0.0;
        residual.noalias()+=H.col(g)*(old-values[g]);
        if (values[g]!=0.0) {deltaGene[g]=1; ++nnGene;}
    }
    const double residualVariance=(residual.array()-residual.mean()).square().sum()/(active.size()-1.0);
    if (residualVariance>0 && std::isfinite(residualVariance)) vareMed=residualVariance;
    VectorXd prediction=H*values;
    const double genetic=(prediction.array()-prediction.mean()).square().sum()/(active.size()-1.0);
    propMed=genetic/(genetic+vareMed); ssqGene=values.squaredNorm();
}
