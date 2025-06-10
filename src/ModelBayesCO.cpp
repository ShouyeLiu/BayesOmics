
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


#include "ModelBayesCO.hpp"

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

// AIAO model
void BayesCO::SnpEffects::sampleFromFCAIAO(VectorXd &wbcorr, vector<VectorXd> &wAcorr,const vector<MatrixDat> &Z,const vector<MatrixDat> ZGene, const map<string, vector<int> > &genePheIdxMap, 
    SnpEffectVec &snpEffectVec, EQTLJointVec &eQTLJointVec, MatrixXd &deltaMat, const map<int,string> &gwasSnpIdx2snpIDMap, map<string, int> geneID2IdxMap,
    const map<string,vector<string> > &gwasSnpID2geneIDMap, const map<string, int> cisSnpID2IdxMap, const double sigmaSqBetaNonEqtl,
    SigmaSqMat &sigmaSqMats, const double &piEffEqtl, const double &piEffNonEqtl,const VectorXd varEps, const double &vare) {
    numNonZeros = 0;
    numNonZerosNonEqtl = 0;
    numNonNullEqtl = 0;
    numNonNullSnpTot = 0;
    numNonNullSnpPerGene = 0;
    nnG = 0;
    // sum of square
    ssqNonEqtl = 0;
    ssqBetaEqtl = 0;
    ssqAlphaEqtl = 0;
    
    ghat.setZero(wbcorr.size());
    gwhatMap.clear();
    gwhatGwasMap.clear();
    double numNonZerosEqtl = 0;
    double sample = 0.0, oldSample = 0.0,varRes = 0.0,rhs = 0.0,uhat = 0.0,invLhs = 0.0;
    double logDelta0 = 0.0 , logDelta1 = 0.0, probDelta1 = 0.0;
    Vector2d sampleVec, oldSampleVec,  varResVec,diagVec, rhsVec, uhatVec;
    Matrix2d invLhsMat;
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
        if(gwasSnpID2geneIDMap.find(snpID) == gwasSnpID2geneIDMap.end()){
            if(true){
            ////////////////////////////////////////////////////////////////////
            //////////////////// run BayesC module ///////////////////////
            ////////////////////////////////////////////////////////////////////
            double logPi = logf(piEffNonEqtl);
            double logPiComp = logf(1.0 - piEffNonEqtl);
            geneIdx = 0;
            oldSample  = snpEffectVec[geneIdx]->getValue(snpID);
            varRes     = vare ;
            rhs        = Z[0].col(snpID).dot(wbcorr);
            rhs        = rhs + sZPZ * oldSample;
            rhs        = rhs/varRes;
            invLhs     = 1.0/(sZPZ/vare + 1.0/sigmaSqBetaNonEqtl);
            uhat       = invLhs * rhs;
            logDelta0  = logPiComp;
            logDelta1  = 0.5*(logf(invLhs) - logf(sigmaSqBetaNonEqtl) + uhat*rhs) + logPi;
            probDelta1 = 1.0/(1.0 + expf(logDelta0-logDelta1));
                if (Stat::ranf() < probDelta1){
                    sample = uhat + Stat::snorm()*sqrt(invLhs);
                    snpEffectVec[geneIdx]->setValue(snpID,sample);
                    wbcorr = wbcorr + Z[0].col(snpID) * (oldSample - sample);
                    ssqNonEqtl += sample*sample;
                    deltaMat(snpIdx,geneIdx) = 1;
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
                ghat = ghat + Z[0].col(snpID) * betaTotal[snpIdx]; 
            }
        } else{
            if(true){
            ///////////////////////////////////////////////////////////////////
            //////////////////// run BayesCO gene module ////////////////////
            ////////////////////////////////////////////////////////////////////
            double logPiGene = logf(piEffEqtl);
            double logPiGeneComp = logf(1.0 - piEffEqtl);
            geneIDSet = gwasSnpID2geneIDMap.at(snpID);
            bool isNonNullSnp = false;
            vector<int> genePheIdxPerGene;
            // cout << "nGWAS: " << nGWAS << endl;
            for(unsigned j = 0; j < geneIDSet.size(); j++) {
                geneID = geneIDSet[j];
                geneIdx  = geneID2IdxMap.at(geneID);
                cisSnpIdx = cisSnpID2IdxMap.at(snpID);
                // cout << "geneID :" << geneID << " snpID: " << snpID << endl;
                genePheIdxPerGene = genePheIdxMap.at(geneID);  // individuals with gene info
                if(gwhatMap.find(geneIdx) == gwhatMap.end()) gwhatMap.insert(pair<int,VectorXd>(geneIdx,VectorXd::Zero(genePheIdxPerGene.size())));  // cis-heritability
                // cout << "size: " << genePheIdxPerGene.size() << endl;
                if(gwhatGwasMap.find(geneIdx) == gwhatGwasMap.end()) gwhatGwasMap.insert(pair<int,VectorXd>(geneIdx,VectorXd::Zero(nGWAS)));  //med heritabiltiy
                double sZPZGene = (ZGene[0].col(snpID)(genePheIdxPerGene)).dot(ZGene[0].col(snpID)(genePheIdxPerGene));
                // cout << " snpID: " << snpID << std::setprecision (15) <<  ZGene[0].col(snpID)(genePheIdxPerGene) << endl;
                Vector2d oldSampleVec = Vector2d( snpEffectVec[geneIdx + 1]->getValue(snpID), eQTLJointVec[geneIdx]->getValue(snpID) );
                varResVec    = Vector2d( vare, varEps[geneIdx] );
                rhsVec       = Vector2d(Z[0].col(snpID).dot(wbcorr), ZGene[0].col(snpID)(genePheIdxPerGene).dot(wAcorr[geneIdx]) );
                rhsVec       = rhsVec + Vector2d(sZPZ * oldSampleVec(0), sZPZGene * oldSampleVec(1)  );
                rhsVec       = rhsVec.array() / varResVec.array();
                diagVec      = Vector2d(sZPZ/vare, sZPZGene/varEps[geneIdx]);
                invLhsMat    =  (diagVec.array()).matrix().asDiagonal();
                invLhsMat    = (invLhsMat + sigmaSqMats.sigmaSqInvMats[geneIdx] ).inverse();
                uhatVec      = invLhsMat.transpose() * rhsVec;
                logDelta0    = logPiGeneComp;
                logDelta1    = 0.5 * (logf(invLhsMat.determinant()) - sigmaSqMats.sigmaSqDetLogVec[geneIdx] + uhatVec.dot(rhsVec) ) + logPiGene;
                probDelta1   = 1.0/(1.0 + expf(logDelta0-logDelta1));
                if (Stat::ranf() < probDelta1) {
                    Vector2d sampleVec;
                    uhatVecArma = arma::dmat(uhatVec.data(),uhatVec.rows(),uhatVec.cols(),false,false);
                    arma::dmat invLhsMatArma = arma::dmat(invLhsMat.data(),invLhsMat.rows(),invLhsMat.cols(),false,false);
                    arma::dvec sampleVecArma = arma::mvnrnd(uhatVecArma, invLhsMatArma, 1); 
                    sampleVec = Eigen::Map<Eigen::VectorXd>(sampleVecArma.memptr(),sampleVecArma.n_rows,sampleVecArma.n_cols); 
                    snpEffectVec[geneIdx +1]->setValue(snpID,sampleVec(0) );
                    eQTLJointVec[geneIdx]->setValue(snpID,sampleVec(1));
                    wbcorr = wbcorr + Z[0].col(snpID) * (oldSampleVec(0) - sampleVec(0) );
                    wAcorr[geneIdx] = wAcorr[geneIdx] + ZGene[0].col(snpID)(genePheIdxPerGene) *(oldSampleVec(1) - sampleVec(1) );

                    gwhatMap.at(geneIdx) += ZGene[0].col(snpID)(genePheIdxPerGene) * eQTLJointVec[geneIdx]->getValue(snpID); // cis-heritability
                    gwhatGwasMap.at(geneIdx) += Z[0].col(snpID) * eQTLJointVec[geneIdx]->getValue(snpID); //med-heri
                    ssqEqtlMat[geneIdx] += ( sampleVec * sampleVec.transpose());
                    ssqBetaEqtl = ssqBetaEqtl + sampleVec(0) * sampleVec(0) ;
                    ssqAlphaEqtl = ssqAlphaEqtl + sampleVec(1) * sampleVec(1);
                    deltaMat(snpIdx,geneIdx+1) = 1;
                    numNonZerosEqtl ++;
                    numNonZeros++;
                    numNonZerosEqtlVecAcrossGenesPostIW[geneIdx] ++;
                    isNonNullSnp = true;
                    nnGeneNameSet.insert(geneID);
                } else {
                    if (oldSampleVec[0]){
                            wbcorr = wbcorr + Z[0].col(snpID) * oldSampleVec(0);
                        }
                    if(oldSampleVec[1]){
                            wAcorr[geneIdx] = wAcorr[geneIdx] + ZGene[0].col(snpID)(genePheIdxPerGene) * oldSampleVec(1);
                        }
                    snpEffectVec[geneIdx +1]->setValue(snpID,0 );
                    eQTLJointVec[geneIdx]->setValue(snpID, 0 );
                }
            }
            betaTotal[snpIdx] += snpEffectVec[geneIdx +1]->getValue(snpID);
            ghat = ghat + Z[0].col(snpID) * betaTotal[snpIdx]; 
            if(isNonNullSnp) {
                ++numNonNullEqtl;
                ++numNonNullSnpTot;    
            }
            } // end of SBayesCO module

        } // end of if statement
    } // end of snp loop

    nnG = nnGeneNameSet.size();
    if(nnG != 0){
        numNonNullSnpPerGene = numNonZerosEqtl / nnG; // the 
    }

    numNonZerosEqtlVec = Vector2d(numNonZerosEqtl,numNonZerosEqtl);
    values = betaTotal;
} // end of sampleFromFCAIAO


// EIEO model
void BayesCO::SnpEffects::sampleFromFCEIEO(VectorXd &wbcorr, vector<VectorXd> &wAcorr,const vector<MatrixDat> &Z,const vector<MatrixDat> ZGene, const map<string, vector<int> > &genePheIdxMap, 
    SnpEffectVec &snpEffectVec, EQTLJointVec &eQTLJointVec, SnpEffectVec &snpEffectVecLatent, EQTLJointVec &eQTLJointVecLatent,
    MatrixXd &deltaMat, const map<int,string> &gwasSnpIdx2snpIDMap, map<string, int> geneID2IdxMap, DeltaVec deltaVecGWAS,DeltaVec deltaVecEQTL,
    const map<string,vector<string> > &gwasSnpID2geneIDMap, const map<string, int> cisSnpID2IdxMap, const double sigmaSqBetaNonEqtl,
    SigmaSqMat &sigmaSqMats, const double &piEffEieo1, const double &piEffEieo2, const double &piEffNonEqtl,const VectorXd varEps, const double &vare) {

    numNonZeros = 0;
    numNonZerosEqtlVec = Vector2d::Zero();
    numSnpCompVec = Vector4d::Zero();
    numNonZerosNonEqtl = 0;
    numNonNullEqtl = 0;
    numNonNullSnpTot = 0;
    numNonNullSnpPerGene = 0;
    ssqBetaTotalGenic = 0;
    numNonNullBetaTotGenic = 0;
    nnG = 0;
    // sum of square
    ssqNonEqtl = 0;
    ssqBetaEqtl = 0;
    ssqAlphaEqtl = 0;
    betaTotal.setZero(gwasSnpIdx2snpIDMap.size());
    betaTotalLatent.setZero(gwasSnpIdx2snpIDMap.size());

    ghat.setZero(wbcorr.size());
    gwhatMap.clear();
    gwhatGwasMap.clear();
    double numNonZerosEqtl = 0;
    double sample = 0.0, oldSample = 0.0,varRes = 0.0,rhs = 0.0,uhat = 0.0,invLhs = 0.0;
    double logDelta0 = 0.0 , logDelta1 = 0.0, probDelta1 = 0.0;
    Vector2d sampleVec, oldSampleVec,  varResVec,diagVec, rhsVec, uhatVec;
    Matrix2d invLhsMat;
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

    int nGenes = geneID2IdxMap.size();
    ssqEqtlMat.resize(nGenes);
    // #pragma omp parallel for schedule(dynamic)
    numNonZerosEqtlVecAcrossGenesPostIW.resize(nGenes);
    for(unsigned j = 0; j < nGenes; ++j){
    ssqEqtlMat[j] = Matrix2d::Zero();
    numNonZerosEqtlVecAcrossGenesPostIW[j] = 0;
    }
    // numNonZerosEqtlVecAcrossGenesPostIW.assign(nGenes,0);

    for(unsigned snpIdx = 0; snpIdx < numTotalSnps; snpIdx ++ ){
        snpID  = gwasSnpIdx2snpIDMap.at(snpIdx);
        double sZPZ = Z[0].col(snpID).dot(Z[0].col(snpID)); // single value of ZPZ from a snp
        if(gwasSnpID2geneIDMap.find(snpID) == gwasSnpID2geneIDMap.end()){
            if(true){
                ////////////////////////////////////////////////////////////////////
                //////////////////// run BayesC module ///////////////////////
                ////////////////////////////////////////////////////////////////////
                double logPi = logf(piEffNonEqtl);
                double logPiComp = logf(1.0 - piEffNonEqtl);
                geneIdx = 0;
                oldSample  = snpEffectVec[geneIdx]->getValue(snpID); 
                varRes     = vare ;
                rhs        = Z[0].col(snpID).dot(wbcorr);
                rhs        = rhs + sZPZ * oldSample;
                rhs        = rhs/varRes;
                invLhs     = 1.0/(sZPZ/vare + 1.0/sigmaSqBetaNonEqtl);
                uhat       = invLhs * rhs;
                logDelta0  = logPiComp;
                logDelta1  = 0.5*(logf(invLhs) - logf(sigmaSqBetaNonEqtl) + uhat*rhs) + logPi;
                probDelta1 = 1.0/(1.0 + expf(logDelta0-logDelta1));
                if (Stat::ranf() < probDelta1){
                sample = uhat + Stat::snorm()*sqrt(invLhs);
                snpEffectVec[geneIdx]->setValue(snpID,sample);
                wbcorr = wbcorr + Z[0].col(snpID) * (oldSample - sample);
                ssqNonEqtl += sample*sample;
                deltaMat(snpIdx,geneIdx) = 1;
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
                //////////////////// run BayesCO gene module ////////////////////
                ////////////////////////////////////////////////////////////////////
                Vector2d logPiGeneVec = Vector2d( log(piEffEieo1), log(piEffEieo2));
                Vector2d logPiGeneCompVec = Vector2d(log(1- piEffEieo1), log(1-piEffEieo2));
                bool isNonNullSnp = false;
                bool isNonNullSnpTrait = false;
                bool isNonNullSnpGene = false;
                geneIDSet = gwasSnpID2geneIDMap.at(snpID);
                vector<int> genePheIdxPerGene;
                for(unsigned j = 0; j < geneIDSet.size(); j++) {
                    string geneID = geneIDSet[j];
                    int geneIdx  = geneID2IdxMap.at(geneID);
                    int cisSnpIdx = cisSnpID2IdxMap.at(snpID);
                    genePheIdxPerGene = genePheIdxMap.at(geneID);  // individuals with gene info
                    if(gwhatMap.find(geneIdx) == gwhatMap.end()) gwhatMap.insert(pair<int,VectorXd>(geneIdx,VectorXd::Zero(genePheIdxPerGene.size())));  // cis-heritability
                    if(gwhatGwasMap.find(geneIdx) == gwhatGwasMap.end()) gwhatGwasMap.insert(pair<int,VectorXd>(geneIdx,VectorXd::Zero(nGWAS)));  //med heritabiltiy
                    double sZPZGene = (ZGene[0].col(snpID)(genePheIdxPerGene)).dot(ZGene[0].col(snpID)(genePheIdxPerGene));
                    Vector2d sampleLatentVec = Vector2d(snpEffectVecLatent[geneIdx +1]->getValue(snpID), eQTLJointVecLatent[geneIdx]->getValue(snpID));
                    Vector2d oldSampleVec = Vector2d( snpEffectVec[geneIdx + 1]->getValue(snpID), eQTLJointVec[geneIdx]->getValue(snpID));
                    Vector2d newSampleVec = oldSampleVec;
                    Vector2d varResVec    = Vector2d( vare, varEps[geneIdx]);
                    Vector2d rhsVec       =  Vector2d(Z[0].col(snpID).dot(wbcorr), ZGene[0].col(snpID)(genePheIdxPerGene).dot(wAcorr[geneIdx]) );
                    //////////////////////////////////////////////
                    rhsVec       = rhsVec + Vector2d(sZPZ * oldSampleVec(0), sZPZGene * oldSampleVec(1)  );
                    rhsVec       = rhsVec.array() / varResVec.array();
                    diagVec      = Vector2d(sZPZ/vare, sZPZGene/varEps[geneIdx]);
                    //// here use eieo sampler I
                    for(int traitIdx = 0; traitIdx < 2; traitIdx ++){
                        int altTraitIdx = 1- traitIdx;
                        double Ginv11 = sigmaSqMats.sigmaSqInvMats[geneIdx](traitIdx,traitIdx);
                        double Ginv12 = sigmaSqMats.sigmaSqInvMats[geneIdx](traitIdx,altTraitIdx);
                        double C11 = 0;
                        if(traitIdx ==0){
                            C11 = Ginv11 + diagVec(0);
                        } else {
                            C11 = Ginv11 + diagVec(1);
                        }
                        double C12 = Ginv12; // assuming residual covariance = 0
                        /// When delta_{jk} = 0
                        double invLhs0 = 1/ Ginv11;
                        double rhs0 = - Ginv12 * sampleLatentVec[altTraitIdx];
                        double uhat0 = rhs0 * invLhs0;

                        // when delta_{jk} = 1
                        double invLhs1 = 1/C11;
                        double rhs1 = rhsVec[traitIdx] - C12 * sampleLatentVec[altTraitIdx];
                        double uhat1 = rhs1 * invLhs1;

                        /// Sample delta_{jk}
                        double logDelta0 = - 0.5 * (log(Ginv11) - uhat0 * uhat0 * Ginv11) + logPiGeneCompVec[traitIdx];
                        double logDelta1 = - 0.5 * (log(C11) - uhat1 * uhat1 * C11) + logPiGeneVec[traitIdx];
                        double probDelta1   = 1.0/(1.0 + exp(logDelta0-logDelta1));
                        // Sample marker effects
                        if (Stat::ranf() < probDelta1) {
                        // if (Stat::ranf() < 1) {
                            sampleLatentVec(traitIdx) = uhat1 + Stat::snorm()*sqrt(invLhs1);
                            // sampleLatentVec(traitIdx) = uhat1 ;
                            newSampleVec(traitIdx) = sampleLatentVec(traitIdx);
                            if(traitIdx == 0) {
                                isNonNullSnpTrait = true;
                                wbcorr = wbcorr + Z[0].col(snpID) * (oldSampleVec(traitIdx) - newSampleVec(traitIdx));
                                numNonZerosEqtlVec(0) = numNonZerosEqtlVec(0) + 1;
                                // deltaMatGWAS(snpIdx,geneIdx+1) = 1;
                                deltaVecGWAS[geneIdx +1]->setValue(snpID,1.0);
                                numNonZerosGenicGwasVec[geneIdx]++;
                                ssqBetaEqtlPG[geneIdx] += sampleLatentVec(traitIdx) * sampleLatentVec(traitIdx);
                            } else {
                                isNonNullSnpGene = true;
                                wAcorr[geneIdx] = wAcorr[geneIdx] + ZGene[0].col(snpID)(genePheIdxPerGene) * (oldSampleVec(traitIdx) - newSampleVec(traitIdx));
                                numNonZerosEqtlVec(1) = numNonZerosEqtlVec(1) + 1;
                                // deltaMateQTL(snpIdx,geneIdx+1) = 1;      
                                deltaVecEQTL[geneIdx]->setValue(snpID,1.0);
                                numNonZerosGenicEqtlVec[geneIdx]++;
                                // ssqAlphaEqtlPG[geneIdx] += sampleLatentVec(traitIdx) * sampleLatentVec(traitIdx);
                            }   
                        } else {
                            sampleLatentVec(traitIdx) = uhat0 + Stat::snorm()*sqrt(invLhs0);
                            newSampleVec(traitIdx)  = 0;
                            if(oldSampleVec(traitIdx)){
                                if(traitIdx == 0){
                                    wbcorr = wbcorr + Z[0].col(snpID) * (oldSampleVec(traitIdx));
                                } else {
                                    wAcorr[geneIdx] = wAcorr[geneIdx] + ZGene[0].col(snpID)(genePheIdxPerGene) * (oldSampleVec(traitIdx));
                                }
                            } // deltaTrait
                        } // end if statement when sampling marker effects
                        // update marker effects
                        if(traitIdx == 0){
                            snpEffectVec[geneIdx +1]->setValue(snpID,newSampleVec(traitIdx));
                            snpEffectVecLatent[geneIdx+1]->setValue(snpID,sampleLatentVec(traitIdx));
                            // ssqBetaEqtl = ssqBetaEqtl + newSampleVec(traitIdx) * newSampleVec(traitIdx);
                            ssqBetaEqtl += sampleLatentVec(traitIdx) * sampleLatentVec(traitIdx);
                            betaTotal[snpIdx] += snpEffectVec[geneIdx +1]->getValue(snpID);
                            betaTotalLatent[snpIdx] += snpEffectVecLatent[geneIdx +1]->getValue(snpID);
                        } else {
                            eQTLJointVec[geneIdx]->setValue(snpID,newSampleVec(traitIdx));
                            eQTLJointVecLatent[geneIdx]->setValue(snpID, sampleLatentVec(traitIdx));
                            ssqAlphaEqtl = ssqAlphaEqtl + newSampleVec(traitIdx) * newSampleVec(traitIdx);
                            ssqAlphaEqtlPG[geneIdx] += sampleLatentVec(traitIdx) * sampleLatentVec(traitIdx);
                        }
                    } /// end of eieo sampler I
                    Vector2d sampleLatentPairVec = Vector2d(snpEffectVecLatent[geneIdx + 1]->getValue(snpID),eQTLJointVecLatent[geneIdx]->getValue(snpID));
                    ssqEqtlMat[geneIdx] += ( sampleLatentPairVec * sampleLatentPairVec.transpose());
                    numNonZerosEqtlVecAcrossGenesPostIW[geneIdx] = numNonZerosEqtlVecAcrossGenesPostIW[geneIdx] + 1;
                }  // end of overlapping genes
                ////////////////////
                if(isNonNullSnpTrait) {
                    numNonNullSnpTot = numNonNullSnpTot + 1; 
                    numNonZeros++;   
                    ++numNonNullBetaTotGenic;  
                    ssqBetaTotalGenic += betaTotal[snpIdx] * betaTotal[snpIdx];             
                } 
                if(isNonNullSnpGene){
                    numNonNullEqtl = numNonNullEqtl + 1;
                }
                // snp role
                if(!isNonNullSnpTrait && !isNonNullSnpGene){
                    numSnpCompVec(0) = numSnpCompVec(0) + 1;  // snp00
                } else if(isNonNullSnpTrait && !isNonNullSnpGene){
                    numSnpCompVec(1) = numSnpCompVec(1) + 1;  // snp10
                } else if(!isNonNullSnpTrait && isNonNullSnpGene){
                    numSnpCompVec(2) = numSnpCompVec(2) + 1;  // snp01
                } else {
                    numSnpCompVec(3) = numSnpCompVec(3) + 1;  // snp11
                }
            } // end of SBayesCO module
        } // end of if statement
    } // end of snp loop

    nnG = nnGeneNameSet.size();
    if(nnG != 0){
        numNonNullSnpPerGene = numNonZerosEqtl / nnG; // the 
    }
    numNonZerosEqtlVec = Vector2d(numNonZerosEqtl,numNonZerosEqtl);
    values = betaTotal;
} // end of sampleFromFCAIAO




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
        if(numNonZerosEqtlVecAcrossGenesPostIW[i] == 0) {continue;}
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
    H = eQTLJointMat(eQTLCommonCis,all);
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

void BayesCO::setStartVal(Data data){
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
    // aiao sampling
    if(mcmcType == "AIAO"){
        do {
            snpEffects.sampleFromFCAIAO(wbcorr,wAcorr,data.ZDat,data.ZGeneDat, data.genePheIdxMap, snpEffectVec, eQTLJointVec,
            deltaMat[0]->values,data.gwasSnpIdx2snpIDMap,data.geneID2IdxMap,data.gwasSnpID2geneIDMap, data.cisSnpID2IdxMap, 
            sigmaSqBetaNonEqtl.value, sigmaSqMats,piEffEqtl.value,piEffNonEqtl.value,varEps.values,vare.value);
            if (++cnt == 100) LOGGER.e(0,"Error: Zero SNP effect in the model for 100 cycles of sampling");
        }while (snpEffects.numNonZeros == 0);
        if (estimatePi) {
        // if (true) {
            if(data.numKeptGenes != 0){
                piEffEqtl.sampleFromFC(data.numEqtlOverlap, snpEffects.numNonZerosEqtlVec(0));
            }else{
                piEffEqtl.value = 0;
            }
            // intergenic region
            if(data.numNonEqtl != 0){
                piEffNonEqtl.sampleFromFC(data.numNonEqtl, snpEffects.numNonZerosNonEqtl);
            }else{
                piEffNonEqtl.value = 0;
            }
        } // end of estimatePi
        if(data.numKeptGenes != 0){  
            // geneEffectVec.sampleFromeFC(data,iter,snpEffects.betaTotal,eQTLJointMat.values,data.gene2gwasSnpMap,data.gene2cisSnpMap,sigmaSqTheta.value,geneEffectVec.vareMed,piTheta.value, geneEffectVec.deltaGene);
            sigmaSqTheta.sampleFromFC(geneEffectVec.values.dot(geneEffectVec.values),geneEffectVec.nnGene);
            piTheta.sampleFromFC(geneEffectVec.numGenes,geneEffectVec.nnGene);
            sigmaSqMats.sampleFromFCInvWishartGeneral(snpEffects.ssqBetaEqtl, snpEffects.ssqAlphaEqtl, snpEffects.ssqEqtlMat,snpEffects.numNonZerosEqtlVec,snpEffects.numNonZerosEqtlVecAcrossGenesPostIW,false);
        } // end of gene effect estimation
    }

    // vare.sampleFromFC(data.y,data.Z,snpEffects.betaTotal, data.n);
    vare.sampleFromFC(wbcorr,data.n);
    vareMean.value = vare.value;
    if(data.numKeptGenes != 0) {
        varEpsMean.value = varEps.values.mean();
        // varEps.sampleFromFC(data.genePheVec,data.ZGene,eQTLJointVec,data.neQTLVec, data.genePheIdxMap, data.geneEffectNames,data.gene2cisSnpMap);
        varEps.sampleFromFC(wAcorr,data.neQTLVec);
    }
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
    deltaMat[0]->values.setZero();
    if(data.numKeptGenes != 0){
        medHsq.compute(vargGene.value, varg.value,vare.value);
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