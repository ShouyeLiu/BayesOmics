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



# include "ModelBayesRO.hpp"


void BayesRO::SnpEffects::sampleFromFCAIAO(const VectorXd &gamma, VectorXd &wbcorr, vector<VectorXd> &wAcorr, const vector<MatrixDat> &Z, const vector<MatrixDat> ZGene, const map<string, vector<int> > &genePheIdxMap, 
    SnpEffectVec &snpEffectVec, EQTLJointVec &eQTLJointVec, 
        const map<int,string> &gwasSnpIdx2snpIDMap, map<string, int> geneID2IdxMap,const map<string,vector<string> > &gwasSnpID2geneIDMap, const map<string, int> cisSnpID2IdxMap,
        SigmaSqBetaNonEqtl &sigmaSqBetaNonEqtl, SigmaSqMatVec &sigmaSqMatVec,const VectorXd &piEffEqtlVec, const VectorXd &piEffNonEqtlVec, const VectorXd varEps, const double &vare) {
    numNonZeros = 0;
    numNonZerosNonEqtl = 0;
    numNonNullEqtl = 0;
    numNonNullSnpTot = 0;
    numNonNullSnpPerGene = 0;
    nnG = 0;
    // sum of square
    ghat.setZero(wbcorr.size());
    gwhatMap.clear();
    gwhatGwasMap.clear();
    ssqNonEqtl = 0;
    ssqBetaEqtl = 0;
    ssqAlphaEqtl = 0;
    double numNonZerosEqtl = 0;
    double sample = 0.0, oldSample = 0.0,varRes = 0.0,rhs = 0.0,uhat = 0.0,invLhs = 0.0;
    Vector2d sampleVec, oldSampleVec,  varResVec,diagVec, rhsVec;
    unsigned deltaj;
    // Matrix2d invLhsMat;
    // BayesR non-gene part
    VectorXd uhatVec;
    VectorXd invLhsVec;
    // BayesR gene part
    vector<Vector2d > uhatVecVec;
    vector<Matrix2d> invLhsMatVec;
    VectorXd logDeltaVec, probDeltaVec;
    arma::dvec uhatVecArma, sampleVecArma;
    arma::dmat invLhsMatArma;
    vector<int> snpIdxInLD;
    vector<string> geneIDSet;
    int snpIdx,geneIdx,cisSnpIdx;
    string snpID,geneID;
    set<string> nnGeneNameSet;
    nnGeneNameSet.clear();
    std::map<int, std::vector<int>>::iterator ldIter;
    int numTotalSnps = gwasSnpIdx2snpIDMap.size();
    int nGWAS =  Z[0].values.rows();
    unsigned ndists = gamma.size();
    uhatVecVec.resize(ndists);
    invLhsMatVec.resize(ndists);
    logDeltaVec.setZero(ndists);
    probDeltaVec.setZero(ndists);
    // non-gene part
    uhatVec.resize(ndists);
    invLhsVec.resize(ndists);

    int nGenes = geneID2IdxMap.size();
    ssqEqtlMat.resize(nGenes);
    // #pragma omp parallel for schedule(dynamic)
    numNonZerosEqtlVecAcrossGenesPostIW.resize(nGenes);
    for(unsigned j = 0; j < nGenes; ++j){
        ssqEqtlMat[j] = Matrix2d::Zero();
        numNonZerosEqtlVecAcrossGenesPostIW[j] = 0;
    }
    

    for(unsigned snpIdx = 0; snpIdx < numTotalSnps; snpIdx ++ ){
        snpID  = gwasSnpIdx2snpIDMap.at(snpIdx);
        double sZPZ = Z[0].col(snpID).dot(Z[0].col(snpID)); // single value of ZPZ from a snp
        // cout << " snpIdx: " << snpIdx << "snpID: " << snpID << endl;
        if(gwasSnpID2geneIDMap.find(snpID) == gwasSnpID2geneIDMap.end()){
            if(true){
            //////////////////////////////////////////////////////////////////
            ////////////////// run BayesC module ///////////////////////
            //////////////////////////////////////////////////////////////////
            VectorXd logPiEffNonEqtlVec = piEffNonEqtlVec.array().log();
            geneIdx = 0;
            oldSample  =  snpEffectVec[geneIdx]->getValue(snpID); 
            varRes     = vare ;
            rhs        = Z[0].col(snpID).dot(wbcorr);
            rhs        = rhs + sZPZ * oldSample;
            rhs        = rhs/varRes;
            invLhsVec.array() = 1.0/(sZPZ/vare + sigmaSqBetaNonEqtl.sigmaSqNonEqtlInvVec.array());
            invLhsVec[0] = 0;
            uhatVec.array()       = invLhsVec.array() * rhs;
            logDeltaVec.array()  = 0.5*((invLhsVec).array().log() - sigmaSqBetaNonEqtl.sigmaSqBetaNonEqtlLogVec.array() + uhatVec.array() * rhs ) + logPiEffNonEqtlVec.array();
            logDeltaVec[0] = logPiEffNonEqtlVec[0];
            for(unsigned k = 0; k < ndists; k++){
                probDeltaVec[k] = 1.0/(exp(logDeltaVec.array() - logDeltaVec[k] )).sum();
            }
            deltaj = bernoulli.sample(probDeltaVec, Stat::ranf());
            nsnpDistNonEqtl[deltaj] = nsnpDistNonEqtl[deltaj] + 1;
            if (deltaj > 0){
                    sample = uhatVec[deltaj] + Stat::snorm()*sqrt(invLhsVec[deltaj]);
                    snpEffectVec[geneIdx]->setValue(snpID,sample);
                    wbcorr = wbcorr + Z[0].col(snpID) * (oldSample - sample);
                    ssqNonEqtl += sample*sample / gamma[deltaj];
                    numNonZerosNonEqtl ++;
                    numNonZeros++;
                    numNonNullSnpTot++;
            } else {
                    if (oldSample) { 
                        wbcorr = wbcorr + Z[0].col(snpID) * oldSample;
                    }
                    snpEffectVec[geneIdx]->setValue(snpID,0);
                    
            } // end of probDelta1
            betaTotal[snpIdx] = snpEffectVec[geneIdx]->getValue(snpID);
        }

        } else{
            if(true){
            ///////////////////////////////////////////////////////////////////
            //////////////////// run BayesOmics gene module ////////////////////
            ////////////////////////////////////////////////////////////////////
            VectorXd logPiEffEqtlVec = piEffEqtlVec.array().log();
            geneIDSet = gwasSnpID2geneIDMap.at(snpID);
            bool isNonNullSnp = false;
            vector<int> genePheIdxPerGene;
            for(unsigned j = 0; j < geneIDSet.size(); j++) {
                geneID = geneIDSet[j];
                geneIdx  = geneID2IdxMap.at(geneID);
                cisSnpIdx = cisSnpID2IdxMap.at(snpID);
                genePheIdxPerGene = genePheIdxMap.at(geneID);  // individuals with gene info
                if(gwhatMap.find(geneIdx) == gwhatMap.end()) gwhatMap.insert(pair<int,VectorXd>(geneIdx,VectorXd::Zero(genePheIdxPerGene.size())));  // cis-heritability
                if(gwhatGwasMap.find(geneIdx) == gwhatGwasMap.end()) gwhatGwasMap.insert(pair<int,VectorXd>(geneIdx,VectorXd::Zero(nGWAS)));  //med heritabiltiy
                double sZPZGene = (ZGene[0].col(snpID)(genePheIdxPerGene)).dot(ZGene[0].col(snpID)(genePheIdxPerGene));
                Vector2d oldSampleVec = Vector2d( snpEffectVec[geneIdx + 1]->getValue(snpID), eQTLJointVec[geneIdx]->getValue(snpID));
                varResVec    = Vector2d( vare, varEps[geneIdx] );
                rhsVec       = Vector2d(Z[0].col(snpID).dot(wbcorr), ZGene[0].col(snpID)(genePheIdxPerGene).dot(wAcorr[geneIdx]) );
                rhsVec       = rhsVec + Vector2d(sZPZ * oldSampleVec(0), sZPZGene * oldSampleVec(1)  );
                rhsVec       = rhsVec.array() / varResVec.array();
                diagVec      = Vector2d(sZPZ/vare, sZPZGene/varEps[geneIdx]);
                for(unsigned k = 1; k < ndists; k++) {
                    invLhsMatVec[k] = (diagVec.array()).matrix().asDiagonal();
                    invLhsMatVec[k] = (invLhsMatVec[k] + sigmaSqMatVec.sigmaSqInvMatVec[k][geneIdx] ).inverse();
                    uhatVecVec[k]     = invLhsMatVec[k].transpose() * rhsVec;
                    logDeltaVec[k]  = 0.5 * (log(invLhsMatVec[k].determinant()) - sigmaSqMatVec.sigmaSqDetLogVecVec[k][geneIdx] + uhatVecVec[k].dot(rhsVec) ) + logPiEffEqtlVec[k];
                }
                logDeltaVec[0] = logPiEffEqtlVec[0];
                for(unsigned k = 0; k < ndists; k++){
                    probDeltaVec[k] = 1.0/(exp(logDeltaVec.array() - logDeltaVec[k] )).sum();
                }
                deltaj = bernoulli.sample(probDeltaVec,Stat::ranf());
                nsnpDistEqtl[deltaj] = nsnpDistEqtl[deltaj] + 1;
                if (deltaj > 0) {
                        uhatVecArma = arma::dmat(uhatVecVec[deltaj].data(),uhatVecVec[deltaj].rows(),uhatVecVec[deltaj].cols(),false,false);
                        invLhsMatArma = arma::dmat(invLhsMatVec[deltaj].data(),invLhsMatVec[deltaj].rows(),invLhsMatVec[deltaj].cols(),false,false);
                        sampleVecArma = arma::mvnrnd(uhatVecArma, invLhsMatArma, 1);
                        sampleVec = Eigen::Map<Eigen::VectorXd>(sampleVecArma.memptr(),sampleVecArma.n_rows,sampleVecArma.n_cols); 
                        snpEffectVec[geneIdx +1]->setValue(snpID,sampleVec(0));
                        eQTLJointVec[geneIdx]->setValue(snpID,sampleVec(1));
                        wbcorr = wbcorr + Z[0].col(snpID) * (oldSampleVec(0) - sampleVec(0) );
                        wAcorr[geneIdx] = wAcorr[geneIdx] + ZGene[0].col(snpID)(genePheIdxPerGene) *(oldSampleVec(1) - sampleVec(1) );
                        ssqEqtlMat[geneIdx] += ( (sampleVec * sampleVec.transpose()) * 1/ gamma[deltaj]).matrix();
                        ssqBetaEqtl = ssqBetaEqtl + sampleVec(0) * sampleVec(0) / gamma[deltaj];
                        ssqAlphaEqtl = ssqAlphaEqtl + sampleVec(1) * sampleVec(1) / gamma[deltaj];
                        gwhatMap.at(geneIdx) += ZGene[0].col(snpID)(genePheIdxPerGene) * sampleVec(1); // cis-heritability
                        gwhatGwasMap.at(geneIdx) += Z[0].col(snpID) * sampleVec(1); //med-heri
                        numNonZerosEqtl ++;
                        numNonZeros++;
                        isNonNullSnp = true;
                        nnGeneNameSet.insert(geneID);
                } else {
                    if (oldSampleVec[0]){
                            wbcorr = wbcorr + Z[0].col(snpID) * oldSampleVec(0);
                            wAcorr[geneIdx] = wAcorr[geneIdx] + ZGene[0].col(snpID)(genePheIdxPerGene) * oldSampleVec(1);
                        }
                    snpEffectVec[geneIdx +1]->setValue(snpID,0);
                    eQTLJointVec[geneIdx]->setValue(snpID,0);
                }
            }
            if(isNonNullSnp) {
                ++numNonNullEqtl;
                ++numNonNullSnpTot;    
            }
            } // end of SBayesOmics module
            } // end of if statement
        // ghat = ghat + Z[0].col(snpID) * snpEffectMat.row(snpIdx).sum(); // heritability
    } // end of snp loop

    nnG = nnGeneNameSet.size();
    if(nnG != 0){
        numNonNullSnpPerGene = numNonZerosEqtl / nnG; // the 
    }

    numNonZerosEqtlVec = Vector2d(numNonZerosEqtl,numNonZerosEqtl);
    values = betaTotal;
}

void BayesRO::ProbMixComps::sampleFromFC(const VectorXd &snpStore) {
	VectorXd dirx;
	dirx = snpStore + alphaVec;
    values = Dirichlet::sample(ndist, dirx);
    for (unsigned i=0; i<ndist; ++i) {
      (*this)[i]->value=values[i];  
    }
}

void BayesRO::SigmaSqMatVec::setPrior(const VectorXd &gamma, const double &sigmaSqBetaEqtl, const VectorXd &sigmaSqAlphaVec){
    Matrix2d varcov, geneCor;
    geneCor << 1, 0, 0, 1;

    for(int i = 1; i < ndists; i++){
        for(int j = 0; j < numGenes; ++j){
            varcov << sigmaSqBetaEqtl, 0, 0, sigmaSqAlphaVec[j];
            varcovPriorMatVec[i][j] = varcov * 0.5;
            sigmaSqMatVec[i][j] = gamma[i] * varcov;
            sigmaSqInvMatVec[i][j] = sigmaSqMatVec[i][j].inverse();
            sigmaSqDetLogVecVec[i][j] = log( sigmaSqMatVec[i][j].determinant());
            // cout << "i: " << i  << " geneIdx: " << geneNames[j] << " varcovPriorMatVec:" << sigmaSqMatVec[i][j]  <<  " invMats: " << sigmaSqInvMatVec[i][j] << " logMat: " << sigmaSqDetLogVecVec[i][j] << endl;
        } // end of gene loop
    } // end of gamma loop

    sigmaSqAlphaAll = 0;
    scaleAlphaPM = (nua-2)/nua* sigmaSqAlphaVec;
    scaleBetaEqtl = (nub-2)/nub * sigmaSqBetaEqtl;
    scaleAlphaAll = ((nua-2)/nua * sigmaSqAlphaVec.array()).mean();

    sigmaSqBetaEqtlPM = sigmaSqBetaEqtl;
    sigmaSqAlphaAll = sigmaSqAlphaVec.mean();
    sigmaSqAlphaPM = sigmaSqAlphaVec;
}

void BayesRO::SigmaSqMatVec::sampleFromFCInvWishartGeneral(VectorXd &gamma, const double &ssqBetaEqtl, const double &ssqAlphaEqtl,vector<Matrix2d> ssqEqtlMat, 
                            const VectorXd &numNonZerosEqtlVec,const vector<unsigned> numNonZerosEqtlVecAcrossGenesPostIW,const bool messageBool){
    ////////////////////////////////////////////////////////////////
    ///// Step 1. Construct genome-wise gene region proir
    ////////////////////////////////////////////////////////////////
    sigmaSqBetaEqtlPM = InvChiSq::sample(nub + numNonZerosEqtlVec(0), ssqBetaEqtl + nub * scaleBetaEqtl );
    sigmaSqAlphaAll  = InvChiSq::sample(nua + numNonZerosEqtlVec(1), ssqAlphaEqtl + nua * scaleAlphaAll );
    unsigned ndists = gamma.size();

    Matrix2d varcovPriorsChr;
    double dfIW = 0.0;
    varcovPriorsChr << sigmaSqBetaEqtlPM * 0.5 , 0.0, 0.0, sigmaSqAlphaAll * 0.5 ;
    Matrix2d sampleEigen;
    for(int i =0; i < numGenes; ++i){
        if(numNonZerosEqtlVecAcrossGenesPostIW[i] == 0) {continue;}
        MatrixXd effArmEigen = ssqEqtlMat[i] + varcovPriorsChr;
        dfIW = 4.0 + numNonZerosEqtlVecAcrossGenesPostIW[i];
        arma::dmat effArma = arma::dmat(effArmEigen.data(),effArmEigen.rows(),effArmEigen.cols(),false,false);
        arma::dmat psiParam = effArma ; 
        arma::dmat sample = arma::iwishrnd(psiParam, dfIW);
        MatrixXd sampleEigen = Eigen::Map<Eigen::MatrixXd>(sample.memptr(),sample.n_rows, sample.n_cols);

        for (unsigned k = 1; k < ndists; k ++ ) {
            // extended to multiple distributions
            sigmaSqMatVec[k][i] = gamma[k] * sampleEigen;
            sigmaSqInvMatVec[k][i] = sigmaSqMatVec[k][i].inverse();
            sigmaSqDetLogVecVec[k][i] = log( sigmaSqMatVec[k][i].determinant());
        }
    }
}

void BayesRO::SigmaSqMatVec::sampleFromFCIWCorr(VectorXd &gamma, const double &ssqBetaEqtl, const double &ssqAlphaEqtl,vector<Matrix2d> ssqEqtlMat, 
                            const VectorXd &numNonZerosEqtlVec,const vector<unsigned> numNonZerosEqtlVecAcrossGenesPostIW,const bool messageBool) {
    ////////////////////////////////////////////////////////////////
    ///// Step 1. Construct genome-wise gene region proir
    ////////////////////////////////////////////////////////////////
    sigmaSqBetaEqtlPM = InvChiSq::sample(nub + numNonZerosEqtlVec(0), ssqBetaEqtl + nub * scaleBetaEqtl );
    sigmaSqAlphaAll  = InvChiSq::sample(nua + numNonZerosEqtlVec(1), ssqAlphaEqtl + nua * scaleAlphaAll );
    unsigned ndists = gamma.size();

    Matrix2d varcovPriorsChr;
    double dfIW = 0.0;
    varcovPriorsChr << sigmaSqBetaEqtlPM * 0.5 , 0.0, 0.0, sigmaSqAlphaAll * 0.5 ;
    Matrix2d sampleEigen;
    for(int i =0; i < numGenes; ++i){
        // running mean


        if(numNonZerosEqtlVecAcrossGenesPostIW[i] == 0) {continue;}
        MatrixXd effArmEigen = ssqEqtlMat[i] + varcovPriorsChr;
        dfIW = 4.0 + numNonZerosEqtlVecAcrossGenesPostIW[i];
        arma::dmat effArma = arma::dmat(effArmEigen.data(),effArmEigen.rows(),effArmEigen.cols(),false,false);
        arma::dmat psiParam = effArma ; 
        arma::dmat sample = arma::iwishrnd(psiParam, dfIW);
        MatrixXd sampleEigen = Eigen::Map<Eigen::MatrixXd>(sample.memptr(),sample.n_rows, sample.n_cols);

        for (unsigned k = 1; k < ndists; k ++ ) {
            // extended to multiple distributions
            sigmaSqMatVec[k][i] = gamma[k] * sampleEigen;
            sigmaSqInvMatVec[k][i] = sigmaSqMatVec[k][i].inverse();
            sigmaSqDetLogVecVec[k][i] = log( sigmaSqMatVec[k][i].determinant());
        }
    }

}


void BayesRO::SigmaSqBetaNonEqtl::sampleFromFC(const double snpEffSumSq, const unsigned numSnpEff){
    double dfTilde = df + numSnpEff;
    double scaleTilde = snpEffSumSq + df * scale;
    value = InvChiSq::sample(dfTilde, scaleTilde);
}

void BayesRO::setStartVal(Data data){
    if(data.numKeptGenes) sigmaSqMatVec.setPrior(gamma.values,sigmaSqBetaEqtl.value, sigmaSqAlphaVec.values);
}

void BayesRO::sampleUnknowns(){
    /////////////////////////////////////
    ///// debug real simu situation /////
    /////////////////////////////////////
    unsigned outi = 0;

    static int iter = 0;
    unsigned cnt=0;
    intercept.sampleFromFC(wbcorr,vare.value,data.n);
    interceptEqtl.sampleFromFC(wAcorr,data.genePheIdxMap,varEps.values,data.neQTLVec);
    // aiao sampling
    if(mcmcType == "AIAO"){
        do {
            snpEffects.sampleFromFCAIAO(gamma.values, wbcorr,wAcorr,data.ZDat,data.ZGeneDat, data.genePheIdxMap, snpEffectVec, eQTLJointVec,
                                        data.gwasSnpIdx2snpIDMap,data.geneID2IdxMap,data.gwasSnpID2geneIDMap, data.cisSnpID2IdxMap, sigmaSqBetaNonEqtl, sigmaSqMatVec,
                                        piEffEqtlVec.values,piEffNonEqtlVec.values,varEps.values,vare.value);
            if (++cnt == 100) LOGGER.e(0,"Error: Zero SNP effect in the model for 100 cycles of sampling");
        }while (snpEffects.numNonZeros == 0);
        if (estimatePi) {
        if (true) {
            if(data.numKeptGenes != 0){
                piEffEqtlVec.sampleFromFC(snpEffects.nsnpDistEqtl);
            }else{
                piEffEqtlVec.values.setZero(snpEffects.ndists);
            }
            // intergenic region
            if(data.numNonEqtl != 0){
                piEffNonEqtlVec.sampleFromFC(snpEffects.nsnpDistNonEqtl);
            }else{
                // piEffNonEqtl.value = 0;
                piEffNonEqtlVec.values.setZero(snpEffects.ndists);
            }
        }
        } // end of estimatePi
        if(data.numKeptGenes != 0){  
            // geneEffectVec.sampleFromeFC(data,iter,snpEffects.be/taTotal,eQTLJointMat.values,data.gene2gwasSnpMap,data.gene2cisSnpMap,sigmaSqTheta.value,geneEffectVec.vareMed,piTheta.value, geneEffectVec.deltaGene);
            sigmaSqTheta.sampleFromFC(geneEffectVec.values.dot(geneEffectVec.values),geneEffectVec.nnGene);
            piTheta.sampleFromFC(geneEffectVec.numGenes,geneEffectVec.nnGene);
            sigmaSqMatVec.sampleFromFCInvWishartGeneral(gamma.values,snpEffects.ssqBetaEqtl, snpEffects.ssqAlphaEqtl, snpEffects.ssqEqtlMat,snpEffects.numNonZerosEqtlVec,snpEffects.numNonZerosEqtlVecAcrossGenesPostIW,false);
        } // end of gene effect estimation
    }
    vare.sampleFromFC(wbcorr,data.n);
    vareMean.value = vare.value;
    if(data.numKeptGenes != 0) {
        varEpsMean.value = varEps.values.mean();
        // varEps.sampleFromFC(data.genePheVec,data.ZGene,eQTLJointMat.values,data.neQTLVec, data.genePheIdxMap, data.geneEffectNames,data.gene2cisSnpMap);
    }
    if(data.numNonEqtl != 0) {
        sigmaSqBetaNonEqtl.sampleFromFC(snpEffects.ssqNonEqtl, snpEffects.numNonZerosNonEqtl);
    } else {
        nnzBtw.value = 0;
        sigmaSqBetaNonEqtl.value = 0;
    }
    // summary statistics
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
    // new code
    varg.compute(snpEffects.ghat);
    vargGeneCis.compute(snpEffects.gwhatMap);
    vargGene.compute(data.numKeptInds,snpEffects.gwhatGwasMap,geneEffectVec.values);
    hsq.compute(varg.value, vare.value);
    if(data.numKeptGenes != 0){
        medHsq.compute(vargGene.value, varg.value,vare.value);
        cisHsq.compute(vargGeneCis.values, data.varPhenotypiceQTL);
        sigmaSqAlpha.value = sigmaSqMatVec.sigmaSqAlphaAll;
        sigmaSqBetaEqtl.value = sigmaSqMatVec.sigmaSqBetaEqtlPM;
        // cout << "varbEqtl: " << sigmaSqBetaEqtl.value << endl;
        cisHsqMean.value = cisHsq.values.mean();
    } else {
        medHsq.value = 0;  
        sigmaSqAlpha.value = 0;
        sigmaSqBetaEqtl.value = 0;  
        cisHsqMean.value = 0; 
    }
    ++iter;
}