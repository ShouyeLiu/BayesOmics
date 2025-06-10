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


#include "ModelSBayesRO.hpp"

void ApproxBayesRO::SnpEffects::sampleFromFCAIAO(const VectorXd &gamma, vector <VectorXd> &wcorrBlocks, vector<VectorXd> &wAcorr, const vector <MatrixDat> &Qblocks, const vector<MatrixDat> &Qgene, 
                const map<int,vector<int> > &ldblock2gwasSnpMap, SnpEffectVec &snpEffectVec, EQTLJointVec &eQTLJointVec, DeltaVec deltaVecGWAS, 
                const map<int,string> &gwasSnpIdx2snpIDMap, const map<string, int> &geneID2IdxMap, const map<string,vector<string> > &gwasSnpID2geneIDMap, 
                const map<string, int> &cisSnpID2IdxMap, SigmaSqBetaNonEqtl &sigmaSqBetaNonEqtl, SigmaSqMatVec &sigmaSqMatVec, const VectorXd &piEffEqtlVec, 
                const VectorXd &piEffNonEqtlVec, const VectorXd &nGWAS, const vector<VectorDat> &neQTL,const VectorXd &varEps, const double &vare){
    numNonZeros = 0;
    numNonZerosNonEqtl = 0;
    numNonNullEqtl = 0;
    numNonNullSnpTot = 0;
    numNonNullSnpPerGene = 0;
    ssqBetaTotalGenic = 0;
    numNonNullBetaTotGenic = 0;
    nnG = 0;
    // sum of square
    // ghat.setZero(wbcorr.size());
    betaTotal.setZero(gwasSnpIdx2snpIDMap.size());
    gwhatMap.clear();
    gwhatGwasMap.clear();
    ssqNonEqtl = 0;
    ssqBetaEqtl = 0;
    ssqAlphaEqtl = 0;
    double numNonZerosEqtl = 0;
    double sample = 0.0, oldSample = 0.0,varRes = 0.0,rhs = 0.0,uhat = 0.0,invLhs = 0.0;
    // double logDelta0 = 0.0 , logDelta1 = 0.0, probDelta1 = 0.0;
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
    // int nGWAS =  Z[0].values.rows();
    unsigned ndists = gamma.size();
    uhatVecVec.resize(ndists);
    invLhsMatVec.resize(ndists);
    logDeltaVec.setZero(ndists);
    probDeltaVec.setZero(ndists);
    // non-gene part
    uhatVec.resize(ndists);
    invLhsVec.resize(ndists);

    VectorXd wbcorr;

    int nBlocks = ldblock2gwasSnpMap.size();
    int nGenes = geneID2IdxMap.size();
    ssqEqtlMat.resize(nGenes);
    ssqAlphaEqtlPG.resize(nGenes);
    // #pragma omp parallel for schedule(dynamic)
    numNonZerosEqtlVecAcrossGenesPostIW.setZero(nGenes);
    for(unsigned j = 0; j < nGenes; ++j){
        ssqEqtlMat[j] = Matrix2d::Zero();
        ssqAlphaEqtlPG[j] = 0;
    }

    for(unsigned lbs = 0; lbs < nBlocks; lbs++ ){
        snpIdxInLD = ldblock2gwasSnpMap.at(lbs);
        wbcorr     = wcorrBlocks[lbs];
        for(unsigned i = 0; i < snpIdxInLD.size(); i++){
            snpIdx = snpIdxInLD[i];
            // snpIdx = 106;
            snpID  = gwasSnpIdx2snpIDMap.at(snpIdx);
            if(gwasSnpID2geneIDMap.find(snpID) == gwasSnpID2geneIDMap.end()){
            ////////////////////////////////////////////////////////////////////
            //////////////////// run SBayesC gene module ///////////////////////
            ////////////////////////////////////////////////////////////////////
                VectorXd logPiEffNonEqtlVec = piEffNonEqtlVec.array().log();
                geneIdx = 0;
                oldSample  = snpEffectVec[geneIdx]->getValue(snpID);
                varRes     = vare / nGWAS[snpIdx];
                rhs        = Qblocks[lbs].col(snpID).dot(wbcorr);
                rhs        = rhs + oldSample;
                rhs        = rhs/varRes;
                invLhsVec.array() = 1.0/(nGWAS[snpIdx]/vare + sigmaSqBetaNonEqtl.sigmaSqNonEqtlInvVec.array());
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
                    wbcorr = wbcorr + Qblocks[lbs].col(snpID) * (oldSample - sample);
                    ssqNonEqtl += sample*sample / gamma[deltaj];
                    // deltaMat(snpIdx,geneIdx) = deltaj;
                    numNonZerosNonEqtl ++;
                    numNonZeros++;
                    numNonNullSnpTot++;
                } else {
                    if (oldSample) { 
                        wbcorr = wbcorr + Qblocks[lbs].col(snpID) * oldSample;
                    }
                    snpEffectVec[geneIdx]->setValue(snpID,0);
                }
                betaTotal[snpIdx] = snpEffectVec[geneIdx]->getValue(snpID);
            } else{
            ////////////////////////////////////////////////////////////////////
            //////////////////// run SBayesE gene module ///////////////////////
            ////////////////////////////////////////////////////////////////////
                VectorXd logPiEffEqtlVec = piEffEqtlVec.array().log();
                // double logPiGeneComp = logf(1.0 - piEffEqtl);
                geneIDSet = gwasSnpID2geneIDMap.at(snpID);
                bool isNonNullSnp = false;
                for(unsigned j = 0; j < geneIDSet.size(); j++) {
                    geneID = geneIDSet[j];
                    geneIdx  = geneID2IdxMap.at(geneID);
                    cisSnpIdx = cisSnpID2IdxMap.at(snpID);
                    oldSampleVec = Vector2d( snpEffectVec[geneIdx+1]->getValue(snpID) , eQTLJointVec[geneIdx]->getValue(snpID) );
                    varResVec    = Vector2d( vare/nGWAS[snpIdx], varEps[geneIdx]/neQTL[geneIdx].sval(snpID) );
                    rhsVec       = Vector2d( Qblocks[lbs].col(snpID).dot(wbcorr), Qgene[geneIdx].col(snpID).dot(wAcorr[geneIdx]) );
                    rhsVec       = rhsVec + oldSampleVec;
                    rhsVec       = rhsVec.array() / varResVec.array(); //;;;;; ///
                    diagVec      = Vector2d(nGWAS[snpIdx]/vare, neQTL[geneIdx].sval(snpID)/varEps[geneIdx]);
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
                        // VectorXd sampleVecArma = Stat::MultiNormal::sample(1,uhatVec.cast <double> (),invLhsMat.cast <double> ()).cast  <double>  ();
                        snpEffectVec[geneIdx+1]->setValue(snpID,sampleVec(0));
                        eQTLJointVec[geneIdx]->setValue(snpID,sampleVec(1)); 
                        wbcorr = wbcorr + Qblocks[lbs].col(snpID) * (oldSampleVec(0) - sampleVec(0));
                        wAcorr[geneIdx] = wAcorr[geneIdx] + Qgene[geneIdx].col(snpID) *(oldSampleVec(1) - sampleVec(1));
                        ssqBetaEqtl = ssqBetaEqtl + sampleVec(0) * sampleVec(0) / gamma[deltaj];
                        ssqAlphaEqtl = ssqAlphaEqtl + sampleVec(1) * sampleVec(1) / gamma[deltaj];
                        ssqAlphaEqtlPG[geneIdx] += sampleVec(1) * sampleVec(1) / gamma[deltaj];
                        ssqEqtlMat[geneIdx] += ( 1/ gamma[deltaj]* sampleVec * sampleVec.transpose());
                        betaTotal[snpIdx] += snpEffectVec[geneIdx +1]->getValue(snpID);
                        // deltaMat(snpIdx,geneIdx+1) = deltaj;
                        numNonZerosEqtl ++;
                        numNonZerosEqtlVecAcrossGenesPostIW[geneIdx] ++;
                        numNonZeros++;
                        isNonNullSnp = true;
                        nnGeneNameSet.insert(geneID);
                    } else {
                        if (oldSampleVec[0]){
                            wbcorr = wbcorr  + Qblocks[lbs].col(snpID) * oldSampleVec(0);
                            wAcorr[geneIdx] = wAcorr[geneIdx] + Qgene[geneIdx].col(snpID) * oldSampleVec(1);
                        }
                        snpEffectVec[geneIdx+1]->setValue(snpID,0);
                        eQTLJointVec[geneIdx]->setValue(snpID,0); 
                    }
                    // snpIdxIncis ++; 
                }  // end of overlapping genes
                if(isNonNullSnp) {
                    ++numNonNullEqtl;
                    ++numNonNullSnpTot;
                    ++numNonNullBetaTotGenic;
                    ssqBetaTotalGenic += betaTotal[snpIdx] * betaTotal[snpIdx];
                    
                } // end of isNonNullSnp
            } // end of SBayesOmics module
        } // end of loop for snplist with one ld block
        wcorrBlocks[lbs] = wbcorr;
    }  // end of ld blocks loop
    nnG = nnGeneNameSet.size();
    if(nnG != 0){
        // numNonNullSnpPerGene = numNonZerosEqtl/ double(geneID2IdxMap.size());
        numNonNullSnpPerGene = numNonZerosEqtl / nnG; // the 
    }
    numNonZerosEqtlVec = Vector2d(numNonZerosEqtl,numNonZerosEqtl);
    values = betaTotal;
} // end of ApproxBayesRO::SnpEffects::sampleFromFCAIAO

void ApproxBayesRO::SnpEffects::sampleFromFCEIEO(const int iter,const bool diagnose, const string title,
                const VectorXd &gamma, vector <VectorXd> &wcorrBlocks, vector<VectorXd> &wAcorr, const vector <MatrixDat> &Qblocks, const vector<MatrixDat> &Qgene, 
                const map<int,vector<int> > &ldblock2gwasSnpMap, SnpEffectVec &snpEffectVec, EQTLJointVec &eQTLJointVec, SnpEffectVec &snpEffectVecLatent, EQTLJointVec &eQTLJointVecLatent,
                DeltaVec deltaVecGWAS, DeltaVec deltaVecEQTL,
                const map<int,string> &gwasSnpIdx2snpIDMap, const map<string, int> &geneID2IdxMap, const map<string,vector<string> > &gwasSnpID2geneIDMap, 
                const map<string, int> &cisSnpID2IdxMap, const double &sigmaSqBetaNonEqtl, SigmaSqMatVec &sigmaSqMatVec, const VectorXd &piEffEieo1Vec,const VectorXd &piEffEieo2Vec, 
                const VectorXd &piEffNonEqtlVec, const VectorXd &nGWAS, const vector<VectorDat> &neQTL,const VectorXd &varEps, const double &vare){

}

void ApproxBayesRO::ProbMixComps::sampleFromFC(const VectorXd &snpStore) {
	VectorXd dirx;
	dirx = snpStore + alphaVec;
    values = Dirichlet::sample(ndist, dirx);
    for (unsigned i=0; i<ndist; ++i) {
      (*this)[i]->value=values[i];  
    }
}

void ApproxBayesRO::SigmaSqMatVec::setPrior(const VectorXd &gamma, const double &sigmaSqBetaEqtl, const VectorXd &sigmaSqAlphaVec){
    Matrix2d varcov, geneCor;
    geneCor << 1, 0, 0, 1;

    for(int i = 1; i < ndists; i++){
        for(int j = 0; j < numGenes; ++j){
            varcov << sigmaSqBetaEqtl * 0.5 , 0, 0, sigmaSqAlphaVec[j] * 0.5;
            // double simga12 = sigmaSqBetaEqtl * sigmaSqAlphaVec[i];
            // varcov << sigmaSqBetaEqtl, simga12, simga12, sigmaSqAlphaVec[i];
            varcovPriorMatVec[i][j] = varcov * 0.5;
            sigmaSqMatVec[i][j] = gamma[i] * varcov;
            sigmaSqInvMatVec[i][j] = sigmaSqMatVec[i][j].inverse();
            sigmaSqDetLogVecVec[i][j] = log( sigmaSqMatVec[i][j].determinant());
            // cout << "i: " << i  << " geneIdx: " << geneNames[j] << " varcovPriorMatVec:" << sigmaSqMatVec[i][j]  <<  " invMats: " << sigmaSqInvMatVec[i][j] << " logMat: " << sigmaSqDetLogVecVec[i][j] << endl;
        } // end of gene loop
    } // end of gamma loop

    sigmaSqAlphaAll = 0;
    scaleBetaEqtl =  (nub-2)/nub*  sigmaSqBetaEqtl;
    scaleAlphaAll = (nua-2)/nua* ( sigmaSqAlphaVec.array()).mean();

    sigmaSqBetaEqtlPM = sigmaSqBetaEqtl;
    sigmaSqAlphaAll = sigmaSqAlphaVec.mean();
    sigmaSqAlphaPM = sigmaSqAlphaVec;
}

void ApproxBayesRO::SigmaSqMatVec::sampleFromFCInvWishartGeneral(VectorXd &geneEffects,VectorXd &gamma, const double &ssqBetaEqtl, const double &ssqAlphaEqtl,vector<Matrix2d> ssqEqtlMat, 
                            const VectorXd &numNonZerosEqtlVec,const VectorXd &numNonZerosEqtlVecAcrossGenesPostIW,const bool messageBool){
    unsigned ndists = gamma.size();

    Matrix2d varcovPriorsChr;
    double dfIW = 0.0;
    ////////////////////////////////////////////////////////////////
    ///// Step 1. Construct genome-wise gene region proir
    ////////////////////////////////////////////////////////////////
    sigmaSqBetaEqtlPM = InvChiSq::sample(nub + numNonZerosEqtlVec(0), ssqBetaEqtl + nub * scaleBetaEqtl );
    sigmaSqAlphaAll  = InvChiSq::sample(nua + numNonZerosEqtlVec(1), ssqAlphaEqtl + nua * scaleAlphaAll );
    // sigmaSqAlphaAll = sigmaSqAlphaAll * 1/( numGenes * 1.0 );
    // varcovPriorsChr << sigmaSqBetaEqtlPM * 0.5 , 0.0, 0.0, sigmaSqAlphaAll * 0.5 ;
    // double simga12 = sigmaSqAlphaAll * sigmaSqBetaEqtlPM;
    varcovPriorsChr << sigmaSqBetaEqtlPM *0.5 , 0, 0, sigmaSqAlphaAll * 0.5 ;
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

void ApproxBayesRO::SigmaSqMatVec::sampleFromFCIWCorr(int iter, int burnIn,VectorXd &geneEffects,VectorXd &gamma, const double &ssqBetaEqtl, const double &ssqAlphaEqtl,vector<Matrix2d> ssqEqtlMat, 
                            const double & ssqBetaTotalGenic, const unsigned &numNonNullBetaTotGenic,
                            const VectorXd &numNonZerosEqtlVec,const VectorXd &numNonZerosEqtlVecAcrossGenesPostIW,VectorXd ssqAlphaEqtlPG, const bool messageBool) {
    if(iter == 0) LOGGER << "IW distribution with averaged scale matrix prior with correlation is used." << endl;
    unsigned ndists = gamma.size();
    Matrix2d varcovPriorsChr,diagVariance,Dinv,corrSam;
    double dfIW = 0.0;
    ////////////////////////////////////////////////////////////////
    ///// Step 1. Construct genome-wise gene region proir
    ////////////////////////////////////////////////////////////////
    sigmaSqBetaEqtlPM = InvChiSq::sample(nub + numNonZerosEqtlVec(0), ssqBetaEqtl + 2* nub * scaleBetaEqtl );
    // sigmaSqBetaEqtlPM = InvChiSq::sample(nub + numNonNullBetaTotGenic, ssqBetaTotalGenic + nub * scaleBetaEqtl );
    sigmaSqAlphaAll  = InvChiSq::sample(nua + numNonZerosEqtlVec(1), ssqAlphaEqtl + nua * scaleAlphaAll );
    // sigmaSqAlphaAll = sigmaSqAlphaAll * 1/( numGenes * 1.0 );
    // varcovPriorsChr << sigmaSqBetaEqtlPM * 0.5 , 0.0, 0.0, sigmaSqAlphaAll * 0.5 ;
    // double simga12 = sigmaSqAlphaAll * sigmaSqBetaEqtlPM;
    // varcovPriorsChr << sigmaSqBetaEqtlPM *0.5 , 0, 0, sigmaSqAlphaAll * 0.5 ;
    Matrix2d sampleEigen;
    for(int i =0; i < numGenes; ++i){
        if(numNonZerosEqtlVecAcrossGenesPostIW[i] == 0) {continue;}
        // sigmaSqAlphaPM[i]  = InvChiSq::sample(nua + numNonZerosEqtlVecAcrossGenesPostIW[i], ssqAlphaEqtlPG[i] + nua * scaleAlphaAll );
        double simga12 = 0; // sigmaSqAlphaAll * sigmaSqBetaEqtlPM;
        varcovPriorsChr << sigmaSqBetaEqtlPM * 0.5  , simga12, simga12, sigmaSqAlphaPM[i] * 0.5 ;
        // calculate running mean
        if(iter == 1){
            iwScaleMat[i] = ssqEqtlMat[i];
            scaleMatBetaSigmaSquareVec[i] = sigmaSqBetaEqtlPM;
            scaleMatAlphaSigmaSquareVec[i] = sigmaSqAlphaPM[i];
        }
        if(1 < iter && iter < burnIn){
            iwScaleMat[i] += (ssqEqtlMat[i] - iwScaleMat[i] )/iter;
            scaleMatBetaSigmaSquareVec[i] += (sigmaSqBetaEqtlPM - scaleMatBetaSigmaSquareVec[i] )/iter;
            scaleMatAlphaSigmaSquareVec[i] += (sigmaSqAlphaPM[i] - scaleMatAlphaSigmaSquareVec[i])/iter;

            diagVariance << sqrt(scaleMatBetaSigmaSquareVec[i] ), 0,0, sqrt(scaleMatAlphaSigmaSquareVec[i]);
            Dinv = iwScaleMat[i].diagonal().cwiseSqrt().asDiagonal().inverse();
            corrSam = Dinv *  iwScaleMat[i] * Dinv;
            varcovPriorsChr = diagVariance * corrSam *diagVariance;
        }
        if(iter > burnIn) {
            diagVariance << sqrt(scaleMatBetaSigmaSquareVec[i] ), 0,0, sqrt(scaleMatAlphaSigmaSquareVec[i]);
            Dinv = iwScaleMat[i].diagonal().cwiseSqrt().asDiagonal().inverse();
            corrSam = Dinv *  iwScaleMat[i] * Dinv;
            varcovPriorsChr = diagVariance * corrSam *diagVariance;
        }

        MatrixXd effArmEigen = ssqEqtlMat[i] + varcovPriorsChr;
        dfIW = 4.0 + numNonZerosEqtlVecAcrossGenesPostIW[i];
        arma::dmat effArma = arma::dmat(effArmEigen.data(),effArmEigen.rows(),effArmEigen.cols(),false,false);
        arma::dmat psiParam = effArma; 
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


void ApproxBayesRO::SigmaSqMatVec::sampleFromFCInvWishartSeparationStragety(VectorXd &geneEffects,VectorXd &gamma, const double &ssqBetaEqtl, const double &ssqAlphaEqtl,vector<Matrix2d> ssqEqtlMat, 
                            const VectorXd &numNonZerosEqtlVec, const VectorXd &numNonZerosEqtlVecAcrossGenesPostIW,const bool messageBool){
    unsigned ndists = gamma.size();

    Matrix2d diagVariance,varcovPriorsChr,corrPrior;
    Matrix2d sampleEigen,effArmEigen;
    double dfIW = 0.0;
    ////////////////////////////////////////////////////////////////
    ///// Step 1. Construct genome-wise gene region proir
    ////////////////////////////////////////////////////////////////
    sigmaSqBetaEqtlPM = InvChiSq::sample(nub + numNonZerosEqtlVec(0), ssqBetaEqtl + nub * scaleBetaEqtl );
    sigmaSqAlphaAll  = InvChiSq::sample(nua + numNonZerosEqtlVec(1), ssqAlphaEqtl + nua * scaleAlphaAll );
    diagVariance << sqrt(sigmaSqBetaEqtlPM), 0,0, sqrt(sigmaSqAlphaAll);

    varcovPriorsChr << 1 , 0.0, 0.0, 1;
    for(int i =0; i < numGenes; ++i){
        if(numNonZerosEqtlVecAcrossGenesPostIW[i] == 0) {continue;}
        dfIW = 4.0 + numNonZerosEqtlVecAcrossGenesPostIW[i];
        // sample correlation from inverse wishart
        Matrix2d Dinv = ssqEqtlMat[i].diagonal().cwiseSqrt().asDiagonal().inverse();
        Eigen::Matrix2d corrSam = Dinv * ssqEqtlMat[i] * Dinv + varcovPriorsChr*1.0E-5;
        arma::dmat IMat = arma::dmat(corrSam.data(),corrSam.rows(),corrSam.cols(),false,false);

        arma::dmat armaCorrMatPrior = arma::iwishrnd(IMat, 3.0);
        MatrixXd  corrMatPrior = Eigen::Map<Eigen::MatrixXd>(armaCorrMatPrior.memptr(),armaCorrMatPrior.n_rows, armaCorrMatPrior.n_cols);
        // Construct prior based on separation strategy prior 
        varcovPriorsChr = diagVariance * corrMatPrior *diagVariance;
        effArmEigen = ssqEqtlMat[i] + varcovPriorsChr;
        arma::dmat effArma = arma::dmat(effArmEigen.data(),effArmEigen.rows(),effArmEigen.cols(),false,false);
        // arma::dmat psiParam = arma::dmat(varcovPriorsChr.data(), varcovPriorsChr.rows(), varcovPriorsChr.cols(),false, false);
        arma::dmat psiParam = effArma ; 
        arma::dmat sample = arma::iwishrnd(psiParam, dfIW);
        MatrixXd sampleEigen = Eigen::Map<Eigen::MatrixXd>(sample.memptr(),sample.n_rows, sample.n_cols);
        for (unsigned k = 1; k < ndists; k ++ ) {
            // extended to multiple distributions
            sigmaSqMatVec[k][i] = gamma[k] * sampleEigen;
            sigmaSqInvMatVec[k][i] = sigmaSqMatVec[k][i].inverse();
            sigmaSqDetLogVecVec[k][i] = log( sigmaSqMatVec[k][i].determinant());
        }
        // geneEffects[i] = sampleEigen(0,1)/sampleEigen(1,1) ;


    } // end of gene loop

}

void ApproxBayesRO::SigmaSqBetaNonEqtl::sampleFromFC(const double snpEffSumSq, const unsigned numSnpEff){
    double dfTilde = df + numSnpEff;
    double scaleTilde = snpEffSumSq + df * scale;
    value = InvChiSq::sample(dfTilde, scaleTilde);
}

void ApproxBayesRO::setStartVal(Data data){
    if(data.numKeptGenes)
        sigmaSqMatVec.setPrior(gamma.values,sigmaSqBetaEqtl.value, sigmaSqAlphaVec.values);
}

void ApproxBayesRO::sampleUnknowns(){
    /////////////////////////////////////
    ///// debug real simu situation /////
    /////////////////////////////////////
    unsigned outi = 0;
    static int iter = 0;
    unsigned cnt=0;
    // aiao sampling
    if(mcmcType == "AIAO"){
        do {
            snpEffects.sampleFromFCAIAO(gamma.values, wcorrBlocks ,wAcorr,data.QblocksDat,data.QgeneDat, data.ldblock2gwasSnpMap, 
                                        snpEffectVec,eQTLJointVec,deltaVecGWAS,
                                        data.gwasSnpIdx2snpIDMap,data.geneID2IdxMap,data.gwasSnpID2geneIDMap, data.cisSnpID2IdxMap, sigmaSqBetaNonEqtl, sigmaSqMatVec,
                                        piEffEqtlVec.values,piEffNonEqtlVec.values,data.n,data.neQTLVec,varEps.values,vare.value);
            if (++cnt == 100) LOGGER.e(0,"Error: Zero SNP effect in the model for 100 cycles of sampling");
        }while (snpEffects.numNonZeros == 0 && snpEffects.numNonZerosEqtlVec(1) == 0);
        if (true) {
            if(data.numKeptGenes != 0){
                piEffEqtlVec.sampleFromFC(snpEffects.nsnpDistEqtl);
            }else{
                piEffEqtlVec.values.setZero(snpEffects.ndists);
            }
        } // end of estimatePi
    } else if (mcmcType == "EIEO"){
        do {
            snpEffects.sampleFromFCEIEO(iter,diagnose,data.label,
                gamma.values,wcorrBlocks,wAcorr,data.QblocksDat,data.QgeneDat,
                data.ldblock2gwasSnpMap,snpEffectVec,eQTLJointVec,
                snpEffectVecLatent,eQTLJointVecLatent, deltaVecGWAS, deltaVecEQTL,
                data.gwasSnpIdx2snpIDMap,data.geneID2IdxMap,data.gwasSnpID2geneIDMap,
                data.cisSnpID2IdxMap,sigmaSqBetaNonEqtl.value, sigmaSqMatVec,
                piEffEieo1Vec.values,piEffEieo2Vec.values,piEffNonEqtlVec.values, data.n,data.neQTLVec,varEps.values,vare.value);
            if (++cnt == 100) LOGGER.e(0," Zero SNP effect in the model for 100 cycles of sampling");
        } while (snpEffects.numNonZeros == 0 && snpEffects.numNonZerosEqtlVec(1) ==0);
        if(true){
            if(data.numKeptGenes != 0){
                piEffEieo1.sampleFromFC(data.numEqtlOverlap, snpEffects.numNonZerosEqtlVec(0));
                piEffEieo2.sampleFromFC(data.numEqtlOverlap, snpEffects.numNonZerosEqtlVec(1));
            }else{
                piEffEieo1.value = 0;
                piEffEieo2.value = 0;
            }
        } /// end of pi sampling
    } else {
        LOGGER.e(0,"Wrong mcmc type for SBayesRO model, please select AIAO or EIEO.");
    }
    ///////////////////////////////////////////////////
    // general samplings
    ///// residual sampling
    if(true){
        if(true) vare.sampleFromFC(iter,data.varPhenotypic,wcorrBlocks,data.gwasEffectInBlock,data.QblocksDat, data.n, snpEffects.betaTotal, data.ldblock2gwasSnpMap);
        // vareMean.value = vare.valueVec.mean();
        vareMean.value = vare.value;
        // vare.value = 0.5;
        if(data.numKeptGenes != 0) {
            varEps.sampleFromFC(data.varPhenotypiceQTL,wAcorr,data.eQTLEffAcrossGenes,data.QgeneDat,data.neQTLVec,eQTLJointVec);
            varEpsMean.value = varEps.values.mean();
        }
    }
    // gene effec and genetic variance and covariance
    if(true){
        if(data.numKeptGenes != 0){
            geneEffectVec.sampleFromeFC(data,iter,"",snpEffects.betaTotal,eQTLJointVec,data.gene2gwasSnpMap,data.gwasSnpIdx2snpIDMap,sigmaSqTheta.value,geneEffectVec.vareMed,piTheta.value, geneEffectVec.deltaGene);
            sigmaSqTheta.sampleFromFC(geneEffectVec.values.dot(geneEffectVec.values),geneEffectVec.nnGene);
            piTheta.sampleFromFC(geneEffectVec.numGenes,geneEffectVec.nnGene);
            // sampling variance-covariance matrix
            sigmaSqMatVec.sampleFromFCIWCorr(iter,data.burnIn, geneEffectVec.values,gamma.values,snpEffects.ssqBetaEqtl, snpEffects.ssqAlphaEqtl, snpEffects.ssqEqtlMat,snpEffects.ssqBetaTotalGenic,snpEffects.numNonNullBetaTotGenic,snpEffects.numNonZerosEqtlVec,snpEffects.numNonZerosEqtlVecAcrossGenesPostIW,snpEffects.ssqAlphaEqtlPG,false);  
            // sigmaSqMatVec.sampleFromFCInvWishartGeneral(geneEffectVec.values,gamma.values,snpEffects.ssqBetaEqtl, snpEffects.ssqAlphaEqtl, snpEffects.ssqEqtlMat,snpEffects.numNonZerosEqtlVec,snpEffects.numNonZerosEqtlVecAcrossGenesPostIW,false);       
        } 
    }
    /////////////////////////////////
    //  pi for intergenic
    if(data.numNonEqtl != 0){
        piEffNonEqtlVec.sampleFromFC(snpEffects.nsnpDistNonEqtl);
    }else{
        // piEffNonEqtl.value = 0;
        piEffNonEqtlVec.values.setZero(snpEffects.ndists);
    }
    ////////////////////////////////////////////
    /////// summary of various non-zero effects
    if(data.numKeptGenes != 0){
        if(mcmcType == "AIAO"){
            nnzGen.getValue(snpEffects.numNonZerosEqtlVec(0));
            nnsGen.getValue(snpEffects.numNonNullEqtl);
        } else if (mcmcType == "EIEO"){
            nnsGen.getValue(snpEffects.numNonZerosEqtlVec(1));
            nnEqtl.getValue(snpEffects.numNonNullEqtl);
            nsnp00.getValue(snpEffects.numSnpCompVec(0));
            nsnp10.getValue(snpEffects.numSnpCompVec(1));
            nsnp01.getValue(snpEffects.numSnpCompVec(2));
            nsnp11.getValue(snpEffects.numSnpCompVec(3));
        }
    } 
    // general
    nnsTot.getValue(snpEffects.numNonNullSnpTot);  // nnSnp
    if(data.numKeptGenes != 0){
        nnGene.getValue(geneEffectVec.nnGene);
        if(nnGene.value != 0){
            nnsPG.value = nnsGen.value / nnGene.value; // Average nonZeroEqtl per nonZeroGene nnEqtlPerGene
        } else {
            nnsPG.value = 0;
        }
    } else {
        nnsPG.value = 0;
    }

    /////////////////////////////////////////////
    /////// summary for intergenic region
    if(data.numNonEqtl != 0) {
        sigmaSqBetaNonEqtl.sampleFromFC(snpEffects.ssqNonEqtl, snpEffects.numNonZerosNonEqtl);
    } else {
        nnzBtw.value = 0;
        sigmaSqBetaNonEqtl.value = 0;
    }
    ///////////////////////////////////////
    // general summary
    ///////////////////////////////////////
    varg.compute(snpEffects.betaTotal, data.QblocksDat,data.ldblock2gwasSnpMap);
    vargGene.compute(geneEffectVec.values,eQTLJointVec, data.QgeneDat);
    geneEffectVec.values = geneEffectVec.valuesAdjust;
    vargGeneCis.compute(eQTLJointVec,data.QgeneDat);

    hsq.compute(varg.value, data.varPhenotypic);

    if(data.numKeptGenes != 0){
        medHsq.compute(vargGene.value,data.varPhenotypic);
        cisHsq.compute(vargGeneCis.values, data.varPhenotypiceQTL);
        sigmaSqAlpha.value = sigmaSqMatVec.sigmaSqAlphaAll;
        sigmaSqBetaEqtl.value = sigmaSqMatVec.sigmaSqBetaEqtlPM;
        cisHsqMean.value = cisHsq.values.mean();
    } else {
        medHsq.value = 0;  
        sigmaSqAlpha.value = 0;
        sigmaSqBetaEqtl.value = 0;  
        cisHsqMean.value = 0; 
    }
    ++iter;
}

