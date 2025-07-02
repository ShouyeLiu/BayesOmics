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


#include "ModelSBayesCO.hpp"
#define STATS_ENABLE_ARMA_WRAPPERS
#define STATS_ENABLE_EIGEN_WRAPPERS
#define STATS_GO_INLINE


void ApproxBayesCO::SnpEffects::sampleFromFCAIAO(Data data,const vector<MatrixXd> &QblocksMat,const int iter,const bool diagnose, const string title, 
                vector <VectorXd> &wcorrBlocks, vector<VectorXd> &wAcorr, vector<VectorXd> &wbcorrGene,
                 const vector <MatrixDat> &Qblocks, const vector<MatrixDat> &Qgene, 
                const map<int,vector<int> > &ldblock2gwasSnpMap,
                SnpEffectVec &snpEffectVec, EQTLJointVec &eQTLJointVec,
                const map<int,string> &gwasSnpIdx2snpIDMap, const map<string, int> &geneID2IdxMap, const map<string,vector<string> > &gwasSnpID2geneIDMap, 
                const map<string, int> &cisSnpID2IdxMap, const double &sigmaSqBetaNonEqtl,  SigmaSqMat &sigmaSqMats,SigmaSqMatResidual &sigmaSqMatRes, const double &piEffEqtl, const double &piEffNonEqtl,
                const VectorXd &nGWAS, const vector<VectorDat> &neQTL,const VectorXd &varEps, const VectorXd &vareVec){
    numNonZeros = 0;
    numNonZerosNonEqtl = 0;
    numNonNullEqtl = 0;
    numNonZerosEqtlVec = Vector2d::Zero();
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
    Vector2d meanVec;
    meanVec.setZero();
    set<string> nnGeneNameSet;
    nnGeneNameSet.clear();
    std::map<int, std::vector<int>>::iterator ldIter;
    int numTotalSnps = gwasSnpIdx2snpIDMap.size();

    int nBlocks = ldblock2gwasSnpMap.size();
    vargInt = 0;
    vargGenic = 0;
    whatBlocksInt.resize(nBlocks);
    whatBlocksGen.resize(nBlocks);
    whatBlocks.resize(nBlocks);
    ssqBlocks.resize(nBlocks);
    for (unsigned i=0; i<nBlocks; ++i) {
        whatBlocksInt[i].setZero(wcorrBlocks[i].size());
        whatBlocksGen[i].setZero(wcorrBlocks[i].size());
        whatBlocks[i].setZero(wcorrBlocks[i].size());
        ssqBlocks[i] = 0;
    }
    int nGenes = geneID2IdxMap.size();
    ssqEqtlMat.resize(nGenes);
    numNonZerosGenicGwasVec.setZero(nGenes);
    numNonZerosGenicEqtlVec.setZero(nGenes);
    numNonZerosEqtlVecAcrossGenesPostIW.setZero(nGenes);
    ssqAlphaEqtlPG.resize(nGenes);
    ssqBetaEqtlPG.resize(nGenes);
    for(unsigned j = 0; j < nGenes; ++j){
        ssqEqtlMat[j] = Matrix2d::Zero();
        ssqAlphaEqtlPG[j] = 0;
        ssqBetaEqtlPG[j] = 0;
    }

    /// here we need to output everything to debug
    std::ofstream file1, file2;
    if(diagnose ){
        string outPath = title;
        string outFile;
        // // gwas snplist
        outFile = "snpeffect-cpp-" + to_string(iter) + ".txt";
        file1.open((outPath + outFile).c_str()); // this one 
        if(iter ==0){
            for(unsigned lbs = 0; lbs < nBlocks; lbs++ ){
            outFile = "scocpp.lblocks." + to_string(lbs) + "." + to_string(iter);
            file2.open((outPath + outFile).c_str());
             vector<int> snpIdxInLD = ldblock2gwasSnpMap.at(lbs);
            for(unsigned i = 0; i < snpIdxInLD.size(); i++){
                int snpIdx = snpIdxInLD[i];
                string snpID  = gwasSnpIdx2snpIDMap.at(snpIdx);
                file2 << snpID << " ";
            }
            file2 << endl;
            file2 << Qblocks[lbs].values ;
            file2 << endl;
            file2.close();
            // if(iter ==0){
                outFile = "scocpp.wbcorr." + to_string(lbs) + "." + to_string(iter);
                file2.open((outPath + outFile).c_str());
                file2 << wcorrBlocks[lbs] ;
                file2 << endl;
                file2.close();
                // }
            }
        }
    }

    for(unsigned lbs = 0; lbs < nBlocks; lbs++ ){
        vector<int> snpIdxInLD = ldblock2gwasSnpMap.at(lbs);
        // vare  = vareBlocks[lbs];
        double vare =  vareVec[lbs];
        Ref<const MatrixXd> Qlbs = QblocksMat[lbs];
        Ref<VectorXd> wbcorr      = wcorrBlocks[lbs];
        LDBlockInfo *blockInfo = data.keptLdBlockInfoVec[lbs];
        unsigned blockStart = blockInfo->startSnpIdx;
        unsigned blockEnd   = blockInfo->endSnpIdx;

        // for(unsigned i = 0; i < snpIdxInLD.size(); i++){
        for(unsigned i = blockStart,snpii = 0; i <= blockEnd; i++){
            // int snpIdx = snpIdxInLD[i];
            int snpIdx = i;
            string snpID  = gwasSnpIdx2snpIDMap.at(snpIdx);
            Ref<const VectorXd> Qlbsi = Qlbs.col(i - blockStart);
            if(gwasSnpID2geneIDMap.find(snpID) == gwasSnpID2geneIDMap.end()){
            ////////////////////////////////////////////////////////////////////
            //////////////////// run SBayesC gene module ///////////////////////
            ////////////////////////////////////////////////////////////////////
                double logPi = logf(piEffNonEqtl);
                double logPiComp = logf(1.0 - piEffNonEqtl);
                int geneIdx = 0;
                double oldSample  = snpEffectVec[geneIdx]->getValue(snpID); //1
                double varRes     = vare / nGWAS[snpIdx]; // 2
                // double rhs        = Qblocks[lbs].col(snpID).dot(wbcorr); //3
                double rhs        = Qlbsi.dot(wbcorr); //3
                rhs        = rhs + oldSample;  //4 
                rhs        = rhs/varRes; // 5
                double invLhs     = 1.0/(nGWAS[snpIdx]/vare + 1.0/sigmaSqBetaNonEqtl ); //6
                double uhat       = invLhs * rhs; // 7
                double logDelta0  = logPiComp; // 8
                double logDelta1  = 0.5*(log(invLhs) - log(sigmaSqBetaNonEqtl) + uhat*rhs) + logPi; // 9
                double probDelta1 = 1.0/(1.0 + exp(logDelta0-logDelta1)); // 10
                if (Stat::ranf() < probDelta1){
                // if (Stat::ranf() < 1){
                    double sample = uhat + Stat::snorm()*sqrt(invLhs);
                    snpEffectVec[geneIdx]->setValue(snpID,sample);
                    wbcorr = wbcorr + Qlbsi * (oldSample - sample);
                    ssqNonEqtl += sample*sample;
                    numNonZerosNonEqtl ++;
                    numNonZeros++;
                    numNonNullSnpTot++;
                    ssqBlocks[lbs] += sample*sample; // zhili
                    whatBlocksInt[lbs] += Qlbsi* sample; // zhili
                    whatBlocks[lbs] += Qlbsi* sample;
                } else {
                    if (oldSample) { 
                        wbcorr = wbcorr + Qlbsi * oldSample;
                    }
                    snpEffectVec[geneIdx]->setValue(snpID,0);
                }
                betaTotal[snpIdx] = snpEffectVec[geneIdx]->getValue(snpID);
                if(diagnose){
                    file1 << lbs << "\t" << snpIdx << "\t" << snpID << "\t" << logDelta0 << "\t" << logDelta1 << "\t"
                    << snpEffectVec[geneIdx]->getValue(snpID) << "\t" << ssqNonEqtl << "\t" << sigmaSqBetaNonEqtl << "\t"
                    << uhat << "\t" << invLhs << "\t" << rhs << "\t" << vare << "\t" << piEffNonEqtl << "\t"
                    << wbcorr.squaredNorm() << "\t" << numNonZeros << endl;
                }

            } else{
            ////////////////////////////////////////////////////////////////////
            //////////////////// run SBayesE gene module ///////////////////////
            ////////////////////////////////////////////////////////////////////
                double logPiGene = logf(piEffEqtl);
                double logPiGeneComp = logf(1.0 - piEffEqtl);
                vector<string> geneIDSet = gwasSnpID2geneIDMap.at(snpID);
                bool isNonNullSnp = false;
                for(unsigned j = 0; j < geneIDSet.size(); j++) {
                    string geneID = geneIDSet[j];
                    int geneIdx  = geneID2IdxMap.at(geneID);
                    int cisSnpIdx = cisSnpID2IdxMap.at(snpID);
                    Vector2d oldSampleVec = Vector2d( snpEffectVec[geneIdx + 1]->getValue(snpID), eQTLJointVec[geneIdx]->getValue(snpID));
                    Vector2d varResVec    = Vector2d( vare/nGWAS[snpIdx], varEps[geneIdx]/neQTL[geneIdx].sval(snpID) );
                    Vector2d rhsVec       = Vector2d( Qlbsi.dot(wbcorr), Qgene[geneIdx].col(snpID).dot(wAcorr[geneIdx]) );
                    rhsVec       = rhsVec + oldSampleVec;
                    rhsVec       = rhsVec.array() / varResVec.array(); //
                    Matrix2d invLhsMat    =  (1.0/varResVec.array()).matrix().asDiagonal();
                    invLhsMat    = (invLhsMat + sigmaSqMats.sigmaSqInvMats[geneIdx]).inverse();
                    Vector2d uhatVec      = invLhsMat.transpose() * rhsVec;
                    double logDelta0    = logPiGeneComp;
                    double logDelta1    = 0.5 * (log(invLhsMat.determinant()) - sigmaSqMats.sigmaSqDetLogVec[geneIdx] + uhatVec.dot(rhsVec) ) + logPiGene;
                    double probDelta1   = 1.0/(1.0 + exp(logDelta0-logDelta1));
                    if (Stat::ranf() < probDelta1) {
                    // if (Stat::ranf() < 1) {
                        arma::dvec uhatVecArma = arma::dmat(uhatVec.data(),uhatVec.rows(),uhatVec.cols(),false,false);
                        arma::dmat invLhsMatArma = arma::dmat(invLhsMat.data(),invLhsMat.rows(),invLhsMat.cols(),false,false);
                        arma::dvec sampleVecArma = arma::mvnrnd(uhatVecArma, invLhsMatArma, 1); 
                        VectorXd sampleVec = Eigen::Map<Eigen::VectorXd>(sampleVecArma.memptr(),sampleVecArma.n_rows,sampleVecArma.n_cols); 
                        snpEffectVec[geneIdx +1]->setValue(snpID,sampleVec(0));
                        eQTLJointVec[geneIdx]->setValue(snpID,sampleVec(1));
                        wbcorr = wbcorr + Qlbsi * (oldSampleVec(0) - sampleVec(0));
                        wAcorr[geneIdx] += Qgene[geneIdx].col(snpID) *(oldSampleVec(1) - sampleVec(1));
                        wbcorrGene[geneIdx] += Qgene[geneIdx].col(snpID) *(oldSampleVec(0) - sampleVec(0));
                        numNonZerosEqtlVec(0) ++;
                        numNonZerosEqtlVec(1)++;
                        numNonZerosEqtlVecAcrossGenesPostIW[geneIdx] ++;
                        ssqEqtlMat[geneIdx] += ( sampleVec * sampleVec.transpose()); // data covariance without centered
                        ssqBetaEqtl = ssqBetaEqtl + sampleVec(0) * sampleVec(0) ;
                        ssqAlphaEqtl = ssqAlphaEqtl + sampleVec(1) * sampleVec(1);
                        ssqBetaEqtlPG[geneIdx] += sampleVec(0) * sampleVec(0);
                        ssqAlphaEqtlPG[geneIdx] += sampleVec(1) * sampleVec(1);
                        // betaTotal[snpIdx] += snpEffectVec[geneIdx +1]->getValue(snpID);
                        betaTotal[snpIdx] += sampleVec(0);
                        numNonZerosGenicGwasVec[geneIdx]++;
                        numNonZerosGenicEqtlVec[geneIdx]++;
                        ssqBlocks[lbs] += sampleVec(0)*sampleVec(0); // zhili
                        whatBlocks[lbs] += Qlbsi* sampleVec(0); // zhili
                        whatBlocksGen[lbs]  += Qlbsi* sampleVec(0); // zhili
                        // numNonZeros++;
                        isNonNullSnp = true;
                        nnGeneNameSet.insert(geneID);
                    } else {
                        if (oldSampleVec[0]){
                            wbcorr = wbcorr  + Qlbsi * oldSampleVec(0);
                            wAcorr[geneIdx] = wAcorr[geneIdx] + Qgene[geneIdx].col(snpID) * oldSampleVec(1);
                            wbcorrGene[geneIdx] += Qgene[geneIdx].col(snpID) *(oldSampleVec(0));
                        }
                        snpEffectVec[geneIdx +1]->setValue(snpID,0);
                        eQTLJointVec[geneIdx]->setValue(snpID,0);
                    }
                    if(diagnose){
                        file1 << snpIdx << "\t" << snpID << "\t" << geneIdx << "\t" << geneID << "\t" << logDelta0 << "\t" << logDelta1 << "\t"
                        << snpEffectVec[geneIdx + 1]->getValue(snpID) << "\t" << eQTLJointVec[geneIdx]->getValue(snpID) << "\t" 
                        << ssqBetaEqtl << "\t" << ssqAlphaEqtl << "\t" << ssqEqtlMat[geneIdx].determinant() << "\t" 
                        << wAcorr[geneIdx].squaredNorm() << "\t" << wbcorr.squaredNorm() << "\t" << numNonZeros << endl;
                    }
                }  // end of overlapping genes
                if(isNonNullSnp) {
                    ++numNonNullEqtl;
                    ++numNonNullSnpTot;
                    ++numNonNullBetaTotGenic;
                    numNonZeros++;
                    ssqBetaTotalGenic += betaTotal[snpIdx] * betaTotal[snpIdx];
                    
                } // end of isNonNullSnp
            } // end of SBayesOmics module
            // ghat = ghat + Z[0].col(snpID) * snpEffectMat.row(snpIdx).sum(); // heritability
        } // end of loop for snplist with one ld block
        wcorrBlocks[lbs] = wbcorr;
        vargInt += whatBlocksInt[lbs].squaredNorm();
        vargGenic += whatBlocksGen[lbs].squaredNorm();
    }  // end of ld blocks loop
    file1.close();
    nnG = nnGeneNameSet.size();
    if(nnG != 0){
        // numNonNullSnpPerGene = numNonZerosEqtl/ double(geneID2IdxMap.size());
        numNonNullSnpPerGene = numNonZerosEqtlVec(0) / nnG; // the 
    }
    values = betaTotal;
} // end of ApproxBayesRO::SnpEffects::sampleFromFCAIAO

void ApproxBayesCO::SnpEffects::sampleFromFCEIEO(Data data,const vector<MatrixXd> &QblocksMat,const int iter,const bool diagnose, const string title, vector <VectorXd> &wcorrBlocks, 
                vector<VectorXd> &wAcorr, const vector <MatrixDat> &Qblocks, const vector<MatrixDat> &Qgene, const map<int,vector<int> > &ldblock2gwasSnpMap, 
                SnpEffectVec &snpEffectVec, EQTLJointVec &eQTLJointVec, SnpEffectVec &snpEffectVecLatent, EQTLJointVec &eQTLJointVecLatent,
                const map<int,string> &gwasSnpIdx2snpIDMap, const map<string, int> &geneID2IdxMap, const map<string,vector<string> > &gwasSnpID2geneIDMap, 
                const map<string, int> &cisSnpID2IdxMap, const double &sigmaSqBetaNonEqtl,  SigmaSqMat &sigmaSqMats, const double &piEffEieo1, const double &piEffEieo2, const double &piEffNonEqtl,
                const VectorXd &nGWAS, const vector<VectorDat> &neQTL,const VectorXd &varEps, const VectorXd &vareVec) {
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
    // Matrix2d invLhsMat;
    // BayesR non-gene part

    vector<int> snpIdxInLD;
    vector<string> geneIDSet;
    // int snpIdx,geneIdx,cisSnpIdx;
    // string snpID,geneID;
    set<string> nnGeneNameSet;
    nnGeneNameSet.clear();
    std::map<int, std::vector<int>>::iterator ldIter;
    int numTotalSnps = gwasSnpIdx2snpIDMap.size();

    int nBlocks = ldblock2gwasSnpMap.size();
        vargInt = 0;
    vargGenic = 0;
    whatBlocksInt.resize(nBlocks);
    whatBlocksGen.resize(nBlocks);
    whatBlocks.resize(nBlocks);
    ssqBlocks.resize(nBlocks);
    for (unsigned i=0; i<nBlocks; ++i) {
        whatBlocksInt[i].setZero(wcorrBlocks[i].size());
        whatBlocksGen[i].setZero(wcorrBlocks[i].size());
        whatBlocks[i].setZero(wcorrBlocks[i].size());
        ssqBlocks[i] = 0;
    }

    int nGenes = geneID2IdxMap.size();
    ssqEqtlMat.resize(nGenes);

    numNonZerosEqtlVecAcrossGenesPostIW.setZero(nGenes);
    numNonZerosGenicGwasVec.setZero(nGenes);
    numNonZerosGenicEqtlVec.setZero(nGenes);
    ssqAlphaEqtlPG.resize(nGenes);
    ssqBetaEqtlPG.resize(nGenes);
    for(unsigned j = 0; j < nGenes; ++j){
        ssqEqtlMat[j] = Matrix2d::Zero();
        ssqAlphaEqtlPG[j] = 0;
        ssqBetaEqtlPG[j] = 0;
    }
    /// here we need to output everything to debug
    std::ofstream file1, file2;
    if(diagnose){
        string outPath = title;
        string outFile;
        // // gwas snplist
        outFile = "snpeffect-cpp-" + to_string(iter) + ".txt";
        file1.open((outPath + outFile).c_str()); // this one 
        for(unsigned lbs = 0; lbs < nBlocks; lbs++ ){
            outFile = "scocpp.qblocks." + to_string(lbs);
            file2.open((outPath + outFile).c_str());
            snpIdxInLD = ldblock2gwasSnpMap.at(lbs);
            for(unsigned i = 0; i < snpIdxInLD.size(); i++){
                int snpIdx = snpIdxInLD[i];
                string snpID  = gwasSnpIdx2snpIDMap.at(snpIdx);
                file2 << snpID << " ";
            }
            file2 << endl;
            file2 << Qblocks[lbs].values ;
            file2 << endl;
            file2.close();
            if(iter ==0){
                outFile = "scocpp.wbcorr." + to_string(lbs);
                file2.open((outPath + outFile).c_str());
                file2 << wcorrBlocks[lbs] ;
                file2 << endl;
                file2.close();
            }
        }
    }
    
    for(unsigned lbs = 0; lbs < nBlocks; lbs++ ){
        // vector<int> snpIdxInLD = ldblock2gwasSnpMap.at(lbs);
        // vare  = vareBlocks[lbs];
        double vare =  vareVec[lbs];
        Ref<const MatrixXd> Qlbs = QblocksMat[lbs];
        Ref<VectorXd> wbcorr      = wcorrBlocks[lbs];
        LDBlockInfo *blockInfo = data.keptLdBlockInfoVec[lbs];
        unsigned blockStart = blockInfo->startSnpIdx;
        unsigned blockEnd   = blockInfo->endSnpIdx;
        unsigned blockSize = blockInfo->numSnpInBlock;


        // shuffling the SNP index for faster convergence
        vector<int> snpIndexVec = Gadget::shuffle_index(blockStart, blockEnd);

        //for(unsigned i = blockStart; i <= blockEnd; i++){
        for (unsigned t = 0; t < blockSize; t++) {
            unsigned i = snpIndexVec[t];
        // for(unsigned i = 0; i < snpIdxInLD.size(); i++){
        // for(unsigned i = blockStart,snpii = 0; i <= blockEnd; i++){
            // int snpIdx = snpIdxInLD[i];
            int snpIdx = i;
            SnpInfo *snp = blockInfo->snpInfoVec[i-blockStart];
            if (snp->skip) {
                // valuesPtr[i] = 0.0;
                continue;
            }
            if (badSnps[i]) {
                // valuesPtr[i] = 0.0;
                continue;
            }
            string snpID  = gwasSnpIdx2snpIDMap.at(snpIdx);
            Ref<const VectorXd> Qlbsi = Qlbs.col(i - blockStart);

            if(gwasSnpID2geneIDMap.find(snpID) == gwasSnpID2geneIDMap.end()){
            ////////////////////////////////////////////////////////////////////
            //////////////////// run SBayesC gene module ///////////////////////
            ////////////////////////////////////////////////////////////////////
                double logPi = logf(piEffNonEqtl);
                double logPiComp = logf(1.0 - piEffNonEqtl);
                int geneIdx = 0;
                double oldSample  = snpEffectVec[geneIdx]->getValue(snpID); //1
                double varRes     = vare / nGWAS[snpIdx]; // 2
                // double rhs        = Qblocks[lbs].col(snpID).dot(wbcorr); //3
                double rhs        = Qlbsi.dot(wbcorr); //3
                rhs        = rhs + oldSample;  //4 
                rhs        = rhs/varRes; // 5
                double invLhs     = 1.0/(nGWAS[snpIdx]/vare + 1.0/sigmaSqBetaNonEqtl ); //6
                double uhat       = invLhs * rhs; // 7
                double logDelta0  = logPiComp; // 8
                double logDelta1  = 0.5*(log(invLhs) - log(sigmaSqBetaNonEqtl) + uhat*rhs) + logPi; // 9
                double probDelta1 = 1.0/(1.0 + exp(logDelta0-logDelta1)); // 10
                if (Stat::ranf() < probDelta1){
                // if (Stat::ranf() < 1){
                    double sample = uhat + Stat::snorm()*sqrt(invLhs);
                    snpEffectVec[geneIdx]->setValue(snpID,sample);
                    wbcorr = wbcorr + Qlbsi * (oldSample - sample);
                    ssqNonEqtl += sample*sample;
                    numNonZerosNonEqtl ++;
                    numNonZeros++;
                    numNonNullSnpTot++;
                    ssqBlocks[lbs] += sample*sample; // zhili
                    whatBlocksInt[lbs] += Qlbsi* sample; // zhili
                    whatBlocks[lbs] += Qlbsi* sample; // zhili
                } else {
                    if (oldSample) { 
                        wbcorr = wbcorr + Qlbsi * oldSample;
                    }
                    snpEffectVec[geneIdx]->setValue(snpID,0);
                }
                betaTotal[snpIdx] = snpEffectVec[geneIdx]->getValue(snpID);
                if(diagnose){
                    file1 << lbs << "\t" << snpIdx << "\t" << snpID << "\t" << logDelta0 << "\t" << logDelta1 << "\t"
                    << snpEffectVec[geneIdx]->getValue(snpID) << "\t" << ssqNonEqtl << "\t" << sigmaSqBetaNonEqtl << "\t"
                    << uhat << "\t" << invLhs << "\t" << rhs << "\t" << vare << "\t" << piEffNonEqtl << "\t"
                    << wbcorr.squaredNorm() << "\t" << numNonZeros << endl;
                }
            } else{
            ////////////////////////////////////////////////////////////////////
            //////////////////// run SBayesE gene module ///////////////////////
            ////////////////////////////////////////////////////////////////////
                Vector2d logPiGeneVec = Vector2d( log(piEffEieo1), log(piEffEieo2));
                Vector2d logPiGeneCompVec = Vector2d(log(1- piEffEieo1), log(1-piEffEieo2));
                vector<string> geneIDSet = gwasSnpID2geneIDMap.at(snpID);
                bool isNonNullSnp = false;
                bool isNonNullSnpTrait = false;
                bool isNonNullSnpGene = false;
                for(unsigned j = 0; j < geneIDSet.size(); j++) {
                    string geneID = geneIDSet[j];
                    int geneIdx  = geneID2IdxMap.at(geneID);
                    int cisSnpIdx = cisSnpID2IdxMap.at(snpID);
                    // Vector2d sampleLatentVec = Vector2d(snpEffectMatLatent(snpIdx,geneIdx +1), eQTLJointMatLatent(cisSnpIdx,geneIdx));
                    // Vector2d oldSampleVec = Vector2d( snpEffectMat(snpIdx,geneIdx+1), eQTLJointMat(cisSnpIdx,geneIdx) );
                    Vector2d sampleLatentVec = Vector2d(snpEffectVecLatent[geneIdx +1]->getValue(snpID), eQTLJointVecLatent[geneIdx]->getValue(snpID));
                    Vector2d oldSampleVec = Vector2d( snpEffectVec[geneIdx + 1]->getValue(snpID), eQTLJointVec[geneIdx]->getValue(snpID));
                    Vector2d newSampleVec = oldSampleVec;
                    Vector2d varResVec    = Vector2d( vare/nGWAS[snpIdx], varEps[geneIdx]/neQTL[geneIdx].sval(snpID) );
                    // Vector2d varResVec    = Vector2d( vare, varEps[geneIdx] );
                    Vector2d rhsVec       = Vector2d( Qlbsi.dot(wbcorr), Qgene[geneIdx].col(snpID).dot(wAcorr[geneIdx]) );
                    //////////////////////////////////////////////
                    rhsVec       = rhsVec + oldSampleVec;
                    rhsVec       = rhsVec.array() / varResVec.array(); //;;
                    //// here use eieo sampler I
                    for(int traitIdx = 0; traitIdx < 2; traitIdx ++){
                        int altTraitIdx = 1- traitIdx;
                        double Ginv11 = sigmaSqMats.sigmaSqInvMats[geneIdx](traitIdx,traitIdx);
                        double Ginv12 = sigmaSqMats.sigmaSqInvMats[geneIdx](traitIdx,altTraitIdx);
                        double C11 = 0;
                        if(traitIdx ==0){
                            C11 = Ginv11 + nGWAS[snpIdx]/vare;
                        } else {
                            C11 = Ginv11 + neQTL[geneIdx].sval(snpID)/varEps[geneIdx];
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
                                wbcorr = wbcorr + Qlbsi * (oldSampleVec(traitIdx) - newSampleVec(traitIdx));
                                numNonZerosEqtlVec(0) = numNonZerosEqtlVec(0) + 1;
                                numNonZerosGenicGwasVec[geneIdx]++;
                                ssqBetaEqtlPG[geneIdx] += sampleLatentVec(traitIdx) * sampleLatentVec(traitIdx);
                                whatBlocksGen[lbs]  += Qlbsi* newSampleVec(traitIdx); // zhili
                                ssqBlocks[lbs] += newSampleVec(traitIdx)*newSampleVec(traitIdx); // zhili
                                whatBlocks[lbs] += Qlbsi* newSampleVec(traitIdx); // zhili
                            } else {
                                isNonNullSnpGene = true;
                                wAcorr[geneIdx] = wAcorr[geneIdx] + Qgene[geneIdx].col(snpID) * (oldSampleVec(traitIdx) - newSampleVec(traitIdx));
                                numNonZerosEqtlVec(1) = numNonZerosEqtlVec(1) + 1;
                                numNonZerosGenicEqtlVec[geneIdx]++;
                                // ssqAlphaEqtlPG[geneIdx] += sampleLatentVec(traitIdx) * sampleLatentVec(traitIdx);
                            }   
                        } else {
                            sampleLatentVec(traitIdx) = uhat0 + Stat::snorm()*sqrt(invLhs0);
                            newSampleVec(traitIdx)  = 0;
                            if(oldSampleVec(traitIdx)){
                                if(traitIdx == 0){
                                    wbcorr = wbcorr + Qlbsi * (oldSampleVec(traitIdx));
                                } else {
                                    wAcorr[geneIdx] = wAcorr[geneIdx] + Qgene[geneIdx].col(snpID) * (oldSampleVec(traitIdx));
                                }
                            } 
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
                        // add diagnose 
                        if(diagnose){
                            if(traitIdx == 0){
                                file1 << snpIdx << "\t" << snpID << "\t" << geneIdx << "\t" << geneID << "\t" << logDelta0 << "\t" << logDelta1 << "\t" 
                                << rhs0 << "\t" << invLhs0  << "\t" << rhs1 << "\t" << invLhs1 << "\t";
                            }else {
                                file1  << logDelta0 << "\t" << logDelta1 << "\t" 
                                    << rhs0 << "\t" << invLhs0  << "\t" 
                                    << rhs1 << "\t" << invLhs1 << "\t";
                            }
                        }
                    } /// end of eieo sampler I
                    Vector2d sampleLatentPairVec = Vector2d(snpEffectVecLatent[geneIdx + 1]->getValue(snpID),eQTLJointVecLatent[geneIdx]->getValue(snpID));
                    ssqEqtlMat[geneIdx] += ( sampleLatentPairVec * sampleLatentPairVec.transpose());
                    numNonZerosEqtlVecAcrossGenesPostIW[geneIdx] = numNonZerosEqtlVecAcrossGenesPostIW[geneIdx] + 1;
                    /////////////////////////////////////////////////////////
                    if(diagnose){
                        file1 << snpEffectVec[geneIdx + 1]->getValue(snpID) << "\t" << eQTLJointVec[geneIdx]->getValue(snpID)  << "\t" 
                        << snpEffectVecLatent[geneIdx+1]->getValue(snpID) << "\t" << eQTLJointVecLatent[geneIdx]->getValue(snpID) << "\t"
                        << ssqBetaEqtl << "\t" << ssqAlphaEqtl << "\t" << ssqEqtlMat[geneIdx].determinant() << "\t" 
                        << wAcorr[geneIdx].squaredNorm() << "\t" << wbcorr.squaredNorm() << endl;
                    }
                    // check snp00 
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
            } // end of SBayesOmics module
        } // end of loop for snplist with one ld block
        wcorrBlocks[lbs] = wbcorr;
        vargInt += whatBlocksInt[lbs].squaredNorm();
        vargGenic += whatBlocksGen[lbs].squaredNorm();
    }  // end of ld blocks loop
    // file1.close();
    nnG = nnGeneNameSet.size();
    if(nnG != 0){
        // numNonNullSnpPerGene = numNonZerosEqtl/ double(geneID2IdxMap.size());
        numNonNullSnpPerGene = numNonZerosEqtlVec(1)  / nnG; // the 
    }
    // numNonZeros = numNonZerosEqtlVec(0);
    values = betaTotal;
} // end of ApproxBayesRO::SnpEffects::sampleFromFCAIAO

void ApproxBayesCO::GenotypicVar::compute(VectorXd &betaTotal, const vector <MatrixDat> &Qblocks,const map<int,vector<int> > &ldblock2gwasSnpMap){
        vector<int> snpIdxInLD;
        value = 0;
        int numBlocks = ldblock2gwasSnpMap.size();
        for(unsigned lbs = 0; lbs < numBlocks; lbs++ ){
            vector<int> snpIdxInLD = ldblock2gwasSnpMap.at(lbs);
            VectorXd Qbeta =  Qblocks[lbs].values * betaTotal(snpIdxInLD); 
            value = value +  Qbeta.dot(Qbeta);
        }
}

void ApproxBayesCO::GenotypicVar::compute(int niter,VectorXd &betaTotal, const vector<MatrixXd> &Qblocks,const map<int,vector<int> > &ldblock2gwasSnpMap){
        vector<int> snpIdxInLD;
        value = 0;
        // cout << "iter " << niter << " "; 
        int numBlocks = ldblock2gwasSnpMap.size();
        for(unsigned lbs = 0; lbs < numBlocks; lbs++ ){
            vector<int> snpIdxInLD = ldblock2gwasSnpMap.at(lbs);
            VectorXd Qbeta =  Qblocks[lbs] * betaTotal(snpIdxInLD); 
            // cout << "lbs " << lbs << " " << Qbeta.dot(Qbeta) << " beta " << betaTotal(snpIdxInLD).sum()  << " lbs: " << Qblocks[lbs].sum() << endl;
            // for (unsigned k = 0; k < snpIdxInLD.size();k++ ) cout << " " << snpIdxInLD[k] << " "; cout << endl;
            value = value +  Qbeta.dot(Qbeta);
        }
        // cout << " vg: " << value << endl;
}

void ApproxBayesCO::GenotypicVarGene::compute(VectorXd &geneEffects, EQTLJointVec &eQTLJointVec,const vector<MatrixDat> &Qgene){
        value = 0;
        for(unsigned i = 0; i < numGenes; i++ ){
            VectorXd QgeneATheta = Qgene[i].values * eQTLJointVec[i]->values * geneEffects[i];
            value = value + QgeneATheta.dot(QgeneATheta);
        }                           
}

void ApproxBayesCO::GenotypicVarGeneCis::computeHsqBeta(SnpEffectVec &snpEffectVec,const vector<MatrixDat> &Qgene){
    for(unsigned i = 0; i < numGenes; i++ ){
        VectorXd Qbeta = Qgene[i].values * snpEffectVec[i+1]->values;
        valueBetaVec[i] =  Qbeta.dot(Qbeta);
    }  
}

void ApproxBayesCO::GenotypicVarGeneCis::compute( EQTLJointVec &eQTLJointVec,const vector<MatrixDat> &Qgene){
        for(unsigned i = 0; i < numGenes; i++ ){
            VectorXd QgeneA = Qgene[i].values * eQTLJointVec[i]->values;
            values[i] = QgeneA.dot(QgeneA) ;
        }                        
}

void ApproxBayesCO::ResidualVar::sampleFromFC(int iter,const double &varPhenotypic, const vector <VectorXd> &wcorrBlocks, const vector<VectorXd> &gwasEffectInBlock, 
                                              const vector <MatrixDat> &Qblocks, const VectorXd &nGWAS, const VectorXd &betaTotal,
                                              const map<int,vector<int> > &ldblock2gwasSnpMap){
    if(iter ==0) LOGGER << "Use adjusted SSE to sample residuals." << endl;
    int nobs = 0;
    double sse = 0.0;
    double sample;
    double sseSub = 0.0;
    double sseTotal = 0;
    double sseSubAcrossBlocks = 0;
    double sseTotalAcrossBlocks = 0;

    for(int i = 0; i < numBlocks; i++){
        // extract cis-snp index for given gene
        vector<int> snpIdxInBlockSet = ldblock2gwasSnpMap.at(i);
        // calculate sse
        sseSub = betaTotal(snpIdxInBlockSet).dot(gwasEffectInBlock[i] + Qblocks[i].values.transpose() * wcorrBlocks[i]);
        sseTotal = (varPhenotypic - sseSub) * (double) nGWAS.mean();
        sseSubAcrossBlocks += sseSub;
        if(sseTotal > 0){
            double dfTilde = df + (double) nGWAS.mean();
            double scaleTilde = sseTotal + df*scale;
            double sample = InvChiSq::sample(dfTilde, scaleTilde);
            valueVec[i] = sample;
            // vareAcrossBlock[i] = sample;
        } else {
            valueVec[i] = vary;
        }
    }
    // sseTotalAcrossBlocks = (varPhenotypic - sseSubAcrossBlocks) * (double) nGWAS.mean();
    // if(sseTotalAcrossBlocks > 0){
    //     double dfTilde = df + (double) nGWAS.mean();
    //     double scaleTilde = sseTotalAcrossBlocks + df*scale;
    //     double sample = InvChiSq::sample(dfTilde, scaleTilde);
    //     valueVec.setConstant(sample);
    // }
}

void ApproxBayesCO::ResidualVar::sampleFromFC(int iter,vector<VectorXd> &wcorrBlocks, const vector<VectorXd> &whatBlocks, VectorXd &ssqBlocks, const VectorXd &nGWASblocks, const VectorXd &numEigenvalBlock) {
    if(iter ==0) LOGGER << "Use SBayesRC's constraits to sample residuals." << endl;
    for (unsigned i=0; i<numBlocks; ++i) {
        float sse = wcorrBlocks[i].squaredNorm() * nGWASblocks[i];
        float dfTilde = df + numEigenvalBlock[i];
        float scaleTilde = sse + df*scale;
        float sample = InvChiSq::sample(dfTilde, scaleTilde);
        /// calculate vargBlock
        if (ssqBlocks[i]/(whatBlocks[i].squaredNorm()) > threshold & sample/vary > 0.9) {
            valueVec[i] = sample;
        } else {
            valueVec[i] = vary;
        }
        //cout << "vare " << i << " " << nGWASblocks[i] << " " << sse << " " << dfTilde << " " << values[i] << " " << ssqBlocks[i] << " " << vargBlocks[i] << endl;
    }
    mean = valueVec.mean();
}


void ApproxBayesCO::ResidualVareEQTL::sampleFromFC(const VectorXd &varPhenotypiceQTL,const vector<VectorXd> &wAcorr, const vector<VectorXd> &eQTLEffAcrossGenes, const vector<MatrixDat> &Qgene, const vector<VectorDat> &neQTL, EQTLJointVec &eQTLJointVec){
    int nobs = 0;
    double sse = 0.0,sseSub = 0.0;
    for (int i = 0; i < numGenes; i++){
        nobs = neQTL[i].values.mean();
        // extract cis-snp index for given gene
        // calculate sse
        sseSub = eQTLJointVec[i]->values.dot(eQTLEffAcrossGenes[i] + Qgene[i].values.transpose() * wAcorr[i]);
        // LOGGER << "eQTLJointVec[i]->values.dot(eQTLEffAcrossGenes[i]  " << eQTLEffAcrossGenes[i].sum()  << endl;
        // LOGGER << "Qgene[i].values.transpose() * wAcorr[i] " << (Qgene[i].values.transpose() * wAcorr[i]).sum() << endl;
        sse = (varPhenotypiceQTL[i] - sseSub) * (double) neQTL[i].values.mean();
        // construct parameters for inverse-chi-square distribution
        if(sse > 0){
            double dfTilde = df  + (double) nobs;
            double scaleTilde = sse + df * scales[i];
            values[i] = Stat::InvChiSq::sample(dfTilde, scaleTilde);
        } else {
            values[i] = varPhenotypiceQTL[i];
        }
    }
}

void ApproxBayesCO::SigmaSqBetaNonEqtl::sampleFromFC(const double &snpEffSumSq, const unsigned &numSnpEff){
    double dfTilde = df + numSnpEff;
    double scaleTilde = snpEffSumSq + df * scale;
    value = InvChiSq::sample(dfTilde, scaleTilde);
}
void ApproxBayesCO::SigmaSqBetaNonEqtl::sampleFromFC(const double &snpEffSumSq, const unsigned &numSnpEff,const double varg,const double piEffNonEqtl,const unsigned iter, const unsigned burnIn){
                // here we set running mean
    // double scaleIteri = 0;
    // if(iter > 0 && iter < burnIn){
    //     // if (noscale){
    //         // } else {
    //         scaleIteri = 0.5 * varg / (numNonEqtl*(piEffNonEqtl));
    //             // }
    //         scaleRunningMean += (scaleIteri - scaleRunningMean)/iter; 
    //         // scale = scaleRunningMean;
    //     }
    // if(iter >= burnIn) {scale = scaleRunningMean;}

    double dfTilde = df + numSnpEff;
    double scaleTilde = snpEffSumSq + df * scale;
    value = InvChiSq::sample(dfTilde, scaleTilde);
}

void ApproxBayesCO::SigmaSqBetaNonEqtl::sampleFromFC(const double &snpEffSumSq, const unsigned &numSnpEff,const unsigned hsqBetaInt){
    double dfTilde = df + numSnpEff;
    double scaleTilde = snpEffSumSq + df * hsqBetaInt;
    value = InvChiSq::sample(dfTilde, scaleTilde);
}

void ApproxBayesCO::SigmaSqAlphaVec::sampleFromFC(VectorXd ssqAlphaEqtlPG, const VectorXd &numNonZerosEqtlVecAcrossGenesPostIW,const VectorXd varg,const double piEffEqtl,const unsigned iter, const unsigned burnIn){
    // construct snpEffectSumSq and numSnpEff by using eQTLJointMat 
    double numSnpEff,snpEffSumSq,dfTilde, scaleTilde;
    for(int i =0; i < numGenes; ++i){
        // if(numNonZerosEqtlVecAcrossGenesPostIW[i] == 0) {continue;}
        // running mean for scale
        // double scaleIteri = 0;
        // if(iter > 0 && iter <= burnIn){
        // // if (noscale){
        //     // } else {
        //         scaleIteri = 0.5 * varg[i] / (numEqtlPG[i] *(piEffEqtl));
        //         // if (scaleIteri > 1e-4){
        //             scaleRunMeanVec[i] += (scaleIteri - scaleRunMeanVec[i])/iter; 
        //             scaleVec[i] = scaleRunMeanVec[i];
        //         // }
        //         // }
        // }
        // if(iter > burnIn) {scaleVec[i] = scaleRunMeanVec[i];}
        numSnpEff = numNonZerosEqtlVecAcrossGenesPostIW[i];
        snpEffSumSq = ssqAlphaEqtlPG[i];
        dfTilde = df + numSnpEff;
        scaleTilde = snpEffSumSq + df * scaleVec[i];
        values[i] = InvChiSq::sample(dfTilde, scaleTilde);
    }
}

void ApproxBayesCO::SigmaSqBetaEqtl::sampleFromFC(const double &snpEffSumSq, const unsigned &numSnpEff,const double varg,const double piEffEqtl,const unsigned iter, const unsigned burnIn){
                // here we set running mean
    // double scaleIteri = 0;
    // if(iter > 0 && iter <= burnIn){
    //     // if (noscale){
    //         // } else {
    //         scaleIteri = 0.5 * varg  / (numEqtl * piEffEqtl);
    //             // }
    //         scaleRunningMean += (scaleIteri - scaleRunningMean)/iter; 
    //         // scale = scaleRunningMean;
    //     }
    // if(iter > burnIn) {scale = scaleRunningMean;}

    double dfTilde = df + numSnpEff;
    double scaleTilde = snpEffSumSq + df * scale;
    value = InvChiSq::sample(dfTilde, scaleTilde);
}

void ApproxBayesCO::SigmaSqBetaEqtl::sampleFromFC(VectorXd ssqAlphaEqtlPG, const VectorXd &numNonZerosEqtlVecAcrossGenesPostIW){

}
void ApproxBayesCO::SigmaSqAlphaVec::sampleFromPrior(){
    for(unsigned i = 0; i < geneNames.size();++i){
        values[i] = InvChiSq::sample((*this)[i]->df, (*this)[i]->scale);
    }
}

void ApproxBayesCO::SigmaSqMatResidual::setPrior(const double &gwasVare, const VectorXd &varEps){
    Matrix2d varcov, geneCor;
    geneCor << 1, 0, 0, 1;
    for(int i =0; i  < numGenes; ++i){
        // varcov << gwasVare, 0, 0, varEps[i];
        varcov << 1, 0, 0, 1;
        varcovPriors[i] = varcov ;
        sigmaSqMats[i] = varcov;
        sigmaSqInvMats[i] = sigmaSqMats[i].inverse();
    }
}
void ApproxBayesCO::SigmaSqMatResidual::sampleFromFC(const vector<VectorXd> &wbcorrGene,const vector<VectorXd> &wAcorrGene,const VectorXd &varEps, const VectorXd &vareVec, const map<int,vector<int> > &ldblock2gwasSnpMap,const map<string, int> &geneID2IdxMap,const map<int,string> &gwasSnpIdx2snpIDMap,const map<string,vector<string> > &gwasSnpID2geneIDMap) {
    double dfIW = 0.0;
    MatrixXd  sampleEigen, effArmEigen;
    Matrix2d varcovPriorsSnp;
    int nBlocks = ldblock2gwasSnpMap.size();
    for(unsigned lbs = 0; lbs < nBlocks; lbs++ ){
        vector<int> snpIdxInLD = ldblock2gwasSnpMap.at(lbs);
        double vare =  vareVec[lbs];
        for(unsigned i = 0; i < snpIdxInLD.size(); i++){
            int snpIdx = snpIdxInLD[i];
            string snpID  = gwasSnpIdx2snpIDMap.at(snpIdx);
            if(gwasSnpID2geneIDMap.find(snpID) != gwasSnpID2geneIDMap.end()){
                vector<string> geneIDSet = gwasSnpID2geneIDMap.at(snpID);
                for(unsigned j = 0; j < geneIDSet.size(); j++) {
                    string geneID = geneIDSet[j];
                    int geneIdx  = geneID2IdxMap.at(geneID);
                    int nnz = wbcorrGene[geneIdx].size();
                    Eigen::MatrixXd wcorrRes(nnz, 2);
                    varcovPriorsSnp << vare,0,0, varEps[geneIdx];
                    wcorrRes.col(0) = wbcorrGene[geneIdx];
                    wcorrRes.col(1) = wAcorrGene[geneIdx];
                    effArmEigen = wcorrRes.transpose()* wcorrRes;
                    effArmEigen += varcovPriorsSnp ;
                    //////// Reconstruct the matrix with modified eigenvalues
                    dfIW = 4 + nnz;
                    arma::dmat effArma = arma::dmat(effArmEigen.data(),effArmEigen.rows(),effArmEigen.cols(),false,false);
                    arma::dmat psiParam = effArma ; 
                    arma::dmat sample = arma::iwishrnd(psiParam, dfIW);
                    sampleEigen = Eigen::Map<Eigen::MatrixXd>(sample.memptr(),sample.n_rows, sample.n_cols);
                    /////////////////////////////////////////////////////////////////////
                    sigmaSqMats[snpIdx] = sampleEigen;
                    sigmaSqInvMats[snpIdx] = sigmaSqMats[snpIdx].inverse();
                }
            }  // end of overlapping genes
        } // end of loop for snplist with one ld block
    }  // end of ld blocks loop
}


void ApproxBayesCO::SigmaSqMat::setPrior(const double &sigmaSqBetaEqtl, const VectorXd &sigmaSqAlphaVec){
    Matrix2d varcov, geneCor;
    geneCor << 1, 0, 0, 1;
    for(int i =0; i  < numGenes; ++i){
        // double simga12 = sigmaSqBetaEqtl * sigmaSqAlphaVec[i];
        // double sigma12 = 0;
        varcov << sigmaSqBetaEqtl, 0, 0, sigmaSqAlphaVec[i];
        varcovPriors[i] = varcov;
        sigmaSqMats[i] = varcov;
        sigmaSqInvMats[i] = sigmaSqMats[i].inverse();
        sigmaSqDetLogVec[i] = log( sigmaSqMats[i].determinant());
    }
    
    sigmaSqAlphaAll = 0;
    scaleAlphaPM = (nua-2)/nua* sigmaSqAlphaVec;
    scaleBetaEqtl = (nub-2)/nub * sigmaSqBetaEqtl;
    scaleAlphaAll = ((nua-2)/nua * sigmaSqAlphaVec.array()).mean();

    sigmaSqBetaEqtlPM = sigmaSqBetaEqtl;
    sigmaSqAlphaAll = sigmaSqAlphaVec.mean();
    sigmaSqAlphaPM = sigmaSqAlphaVec;

}

void ApproxBayesCO::SigmaSqMat::sampleFromFCIWIndSMatPriorAIAO(int iter, int burnIn, VectorXd &geneEffects,const double &ssqBetaEqtl,const double &ssqAlphaEqtl, 
         const double & ssqBetaTotalGenic, const unsigned &numNonNullBetaTotGenic,
         const  VectorXd ssqBetaEqtlPG,const  VectorXd ssqAlphaEqtlPG, vector<Matrix2d> ssqEqtlMat,
         const VectorXd &numNonZerosEqtlVec,const VectorXd &numNonZerosEqtlVecAcrossGenesPostIW,
         const double hsqGenic,const VectorXd cisHsqMean ,const bool messageBool){
    if(iter == 0) LOGGER << "IW distribution with independent prior is used." << endl;
    Matrix2d diagVariance,corrPrior,corrMatPrior;
    Matrix2d effArmEigen,sampleEigen;
    sigmaSqBetaEqtlPM = InvChiSq::sample(nub + numNonNullBetaTotGenic, ssqBetaTotalGenic + nub * scaleBetaEqtl );
    sigmaSqAlphaAll  = InvChiSq::sample(nua + numNonZerosEqtlVec(1), ssqAlphaEqtl + nua * scaleAlphaAll );
    // sigmaSqAlphaAll = sigmaSqAlphaPM.mean();
    double scaleAlphaGene ;
    double scaleBetaGene;
    Matrix2d varcovPriorsChr;
    varcovPriorsChr << sigmaSqBetaEqtlPM  , 0.0, 0.0, sigmaSqAlphaAll ;
    double dfIW = 0.0;
    for(int i =0; i < numGenes; ++i){
        if(numNonZerosEqtlVecAcrossGenesPostIW[i] == 0) {continue;}
        if (cisHsqMean[i] != 0){
            // scaleBetaEqtl = hsqGenic / (double) numNonZerosEqtlVec(0);
            // scaleAlphaPM[i] = cisHsqMean[i] / (double) numNonZerosEqtlVecAcrossGenesPostIW[i];
            // sigmaSqAlphaPM[i] = cisHsqMean[i] / (double) numNonZerosEqtlVecAcrossGenesPostIW[i];
        }
        varcovPriorsChr << sigmaSqBetaEqtlPM  , 0.0, 0.0, sigmaSqAlphaPM[i];
        MatrixXd effArmEigen = ssqEqtlMat[i] + varcovPriorsChr;
        dfIW = 4 + numNonZerosEqtlVecAcrossGenesPostIW[i];
        arma::dmat effArma = arma::dmat(effArmEigen.data(),effArmEigen.rows(),effArmEigen.cols(),false,false);
        arma::dmat psiParam = effArma ; 
        arma::dmat sample = arma::iwishrnd(psiParam, dfIW);
        sampleEigen = Eigen::Map<Eigen::MatrixXd>(sample.memptr(),sample.n_rows, sample.n_cols);
        /////////////////////////////////////////////////////////////////////
        sigmaSqMats[i] = sampleEigen;
        sigmaSqInvMats[i] = sigmaSqMats[i].inverse();
        sigmaSqDetLogVec[i] = logf( sigmaSqMats[i].determinant());
        // geneEffects[i] = sampleEigen(0,1)/sampleEigen(1,1) ;
    }
}

void ApproxBayesCO::SigmaSqMat::sampleFromFCIWIndSMatPriorEIEO(int iter, int burnIn, const double &sigmaSqBetaEqtlPM,const VectorXd &sigmaSqAlphaPM,const vector<Matrix2d> & ssqEqtlMat, double &geneticCov, 
         VectorXd &geneEffects, const VectorXd &numEqtlPG,
         const bool messageBool){
    
    if(iter == 0) LOGGER << "IW distribution with independent prior is used." << endl;
    Matrix2d diagVariance,corrPrior,corrMatPrior;
    // if (hsqGenic != 0 && cisHsqMean != 0){
    //     scaleBetaEqtl = hsqGenic / (double) numNonZerosEqtlVec(0);
    //     scaleAlphaAll = cisHsqMean / (double) numNonZerosEqtlVec(1);
    // }
    // sigmaSqBetaEqtlPM = InvChiSq::sample(nub + numNonZerosEqtlVec(0), ssqBetaEqtl + nub * scaleBetaEqtl );
    // sigmaSqAlphaAll  = InvChiSq::sample(nua + numNonZerosEqtlVec(1), ssqAlphaEqtl + nua * scaleAlphaAll );
    // sigmaSqBetaEqtlPM = InvChiSq::sample(nub + numNonNullBetaTotGenic, ssqBetaTotalGenic + nub * scaleBetaEqtl );
    // sigmaSqAlphaAll  = InvChiSq::sample(nua + numNonZerosEqtlVec(1), ssqAlphaEqtl + nua * scaleAlphaAll );
    // sigmaSqAlphaAll = sigmaSqAlphaPM.mean();
    double scaleAlphaGene ;
    double scaleBetaGene;
    
    Matrix2d varcovPriorsChr,sampleEigen,effArmEigen;
    float dfPrior;
    // varcovPriorsChr << sigmaSqBetaEqtlPM *0.5 , 0.0, 0.0, sigmaSqAlphaAll *0.5;
    double dfIW = 0.0;
    for(int i =0; i < numGenes; ++i){
        varcovPriorsChr << sigmaSqBetaEqtlPM , 0.0, 0.0, sigmaSqAlphaPM[i];
        if(true){
            if(iter ==0 && i == 0){LOGGER << "use 10% prior " << endl;}
            // dfPrior = 0.1*numEqtlPG[i];
            // effArmEigen = ssqEqtlMat[i] + varcovPriorsChr * (dfPrior );
            // dfIW = numEqtlPG[i] + dfPrior;
            dfPrior = 4;
            effArmEigen = ssqEqtlMat[i] + varcovPriorsChr * (dfPrior -3 );
            dfIW = numEqtlPG[i] + dfPrior;

        } else {
            dfPrior = 100000;
            effArmEigen = ssqEqtlMat[i] + varcovPriorsChr * (dfPrior );
            dfIW = numEqtlPG[i] + dfPrior;
        }
        arma::dmat effArma = arma::dmat(effArmEigen.data(),effArmEigen.rows(),effArmEigen.cols(),false,false);
        arma::dmat psiParam = effArma ; 
        arma::dmat sample = arma::iwishrnd(psiParam, dfIW);
        sampleEigen = Eigen::Map<Eigen::MatrixXd>(sample.memptr(),sample.n_rows, sample.n_cols);
        // cout << "iter " << iter << " sample\n " << sampleEigen << " \nprior\n " <<  varcovPriorsChr << endl;
        // int tmp;
        // cin >> tmp;
        
        /////////////////////////////////////////////////////////////////////
        // if(i ==0 && iter == 0) LOGGER << " use independent prior" << endl;
        // varcovPriorsChr << sigmaSqBetaEqtlPM , 0.0, 0.0, sigmaSqAlphaPM[i];
        // sampleEigen = varcovPriorsChr;
        // sampleEigen(0,1) = sampleEigen(1,0) = 0;
        /////////////////////////////////////////
        sigmaSqMats[i] = sampleEigen;

        sigmaSqInvMats[i] = sigmaSqMats[i].inverse();
        sigmaSqDetLogVec[i] = logf( sigmaSqMats[i].determinant());
        geneEffects[i] = sampleEigen(0,1)/sampleEigen(1,1) ;
        geneticCov += sampleEigen(0,1)/(sqrt(sampleEigen(0,0))*sqrt(sampleEigen(1,1))); 
    }
}

void ApproxBayesCO::SigmaSqMat::sampleFromFCIWSMatCorrPriorAIAO(int iter, int burnIn, double sigmaSqBetaEqtlPM,VectorXd sigmaSqAlphaPM, VectorXd &geneEffects,const double &ssqBetaEqtl,const double &ssqAlphaEqtl, 
         const double & ssqBetaTotalGenic, const unsigned &numNonNullBetaTotGenic,
         const  VectorXd ssqBetaEqtlPG,const  VectorXd ssqAlphaEqtlPG, vector<Matrix2d> ssqEqtlMat,
         const VectorXd &numNonZerosEqtlVec , const VectorXd &numNonZerosEqtlVecAcrossGenesPostIW,
         const double hsqGenic,const double cisHsqMean ,const bool messageBool) {
    // this prior will lead to convergence isssue for varbEqtl. 
    if(iter == 0) LOGGER << "IW distribution with averaged scale matrix prior with correlation is used." << endl;
    Matrix2d effArmEigen, varcovPriorsChr, sampleEigen;
    double dfIW = 0.0;
    for(int i =0; i < numGenes; ++i){
        if(numNonZerosEqtlVecAcrossGenesPostIW[i] == 0) {continue;}
        varcovPriorsChr << sigmaSqBetaEqtlPM , 0.0, 0.0, sigmaSqAlphaPM[i];
        effArmEigen = ssqEqtlMat[i] + varcovPriorsChr;
        dfIW = 4 + numNonZerosEqtlVecAcrossGenesPostIW[i];
        arma::dmat effArma = arma::dmat(effArmEigen.data(),effArmEigen.rows(),effArmEigen.cols(),false,false);
        arma::dmat psiParam = effArma ; 
        arma::dmat sample = arma::iwishrnd(psiParam, dfIW);
        MatrixXd   sampleEigen = Eigen::Map<Eigen::MatrixXd>(sample.memptr(),sample.n_rows, sample.n_cols);
        /////////////////////////////////////////////////////////////////////
        sigmaSqMats[i] = sampleEigen;
        sigmaSqInvMats[i] = sigmaSqMats[i].inverse();
        sigmaSqDetLogVec[i] = logf( sigmaSqMats[i].determinant());
    }
}

void ApproxBayesCO::GeneEffects::sampleFromeFC(Data data, int iter ,const string title, VectorXd &betaTotal, MatrixXd &eQTLJointMat,
                            const map<int,vector<int> > &gene2gwasSnpMap,const map<int,vector<int>> &gene2cisSnpMap, 
                            double &sigmaSqTheta, double &vareMed, const double &piTheta, VectorXd &deltaGene ){
    double sample, oldSample,varRes,rhs,uhat,invLhs,muHat;
    double logDelta0, logDelta1, probDelta1;
    double sigmaSqThetaInv, sigmaSqThetaLog;
    double muBeta;
    int eQTLNumAcrossGenes =0;
    vector<int> eQTLCommon; // gwas-snp
    vector<int> eQTLCommonCis; // cis-snp
    eQTLCommon.clear();
    eQTLCommonCis.clear();
    deltaGene.setZero();
    VectorXd betaCorr, betaHat;
    MatrixXd H, HPH; 
    set<int> eQTLUniqSet; 
    set<int>::iterator iterSet;
    vector<int> geneIdxVec;
    geneIdxVec.clear();
    bool geneHasSnps;
    std::ofstream file1, file2;
    string outPath = data.label;
    string outFile;
    for(int i =0; i < numGenes; ++i){
        geneHasSnps = false;
        vector<int> snpIdxInGeneSet = gene2gwasSnpMap.at(i);
        vector<int> gene2cisSnpIdxPerGene = gene2cisSnpMap.at(i);
        // string outFile = "geneEffect_alpha_gene_"+ to_string(i) + "_" + to_string(iter) + ".txt";
        // file1.open((outPath + outFile).c_str());
        // file1 << eQTLJointMat(gene2cisSnpIdxPerGene,i) << endl;
        // file1.close();
        for(int j = 0; j < snpIdxInGeneSet.size(); j++){
            int gwasIdx = snpIdxInGeneSet[j];
            int eqtlIdx = gene2cisSnpIdxPerGene[j];
            if(eQTLUniqSet.find(gwasIdx) != eQTLUniqSet.end()){
                continue;
            }
            if (betaTotal[gwasIdx] != 0 && eQTLJointMat(eqtlIdx,i) != 0 ){
                // if(iter == 97){
                //     cout << "betaTotal[gwasIdx]  " << betaTotal[gwasIdx] << " eQTLJointMat(eqtlIdx,i) " << eQTLJointMat(eqtlIdx,i) << endl;
                // }
                eQTLCommon.push_back(gwasIdx);
                eQTLCommonCis.push_back(eqtlIdx);
                eQTLUniqSet.insert(gwasIdx); 
                geneHasSnps = true;
            }
        }
        if(geneHasSnps) geneIdxVec.push_back(i);
    }
    nnGene = 0;
    ssqGene = 0;
    eQTLNumAcrossGenes = eQTLCommon.size();
    sigmaSqThetaInv = 1/sigmaSqTheta;
    sigmaSqThetaLog = logf(sigmaSqTheta);
    betaCorr = betaTotal(eQTLCommon);
    H = eQTLJointMat(eQTLCommonCis,geneIdxVec);
    HPH = (H.transpose() * H).diagonal();
    // // gwas snplist
    // std::ofstream file1, file2;
    // outFile = "geneEffect_beta_" + to_string(iter) + ".txt";
    // file1.open((outPath + outFile).c_str());
    // file1 << betaCorr << endl;
    // file1.close();
    // outFile = "geneEffect_H_" + to_string(iter) + ".txt";
    // file1.open((outPath + outFile).c_str());
    // file1 << H << endl;
    // file1.close();
    // if(eQTLNumAcrossGenes < 9){
    //     cout << "here" << endl;
    // }
    if(eQTLNumAcrossGenes >= 2){
        muBeta = betaCorr.mean();
        betaCorr.array() = betaCorr.array() - muBeta - (H * values(geneIdxVec)).array();
        /// Start mcmc process
        for(unsigned iteri =0; iteri < 1; iteri++){
            betaCorr.array() += muBeta;
            rhs = betaCorr.sum();
            invLhs = 1/(double) eQTLNumAcrossGenes;
            muHat = invLhs*rhs;
            muBeta = Normal::sample(muHat,invLhs*vareMed);
            betaCorr.array() -= muBeta;
            betaHat.setZero(eQTLNumAcrossGenes); // betaHat = 0;
            for (unsigned i = 0; i < geneIdxVec.size(); i++) {
                int geneIdx = geneIdxVec[i];
                double betaSEScale = Gadget::calcMean(data.scalingGWASFactorVec);
                double alphaSEScale = Gadget::calcMean(data.scalingeQTLFactorVecVec[geneIdx]);
                double logPi = logf(piTheta);
                double logPiComp = logf(1-piTheta);
                oldSample = values[geneIdx];
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
                    sample = uhat + Stat::snorm()*sqrtf(invLhs);
                    values[geneIdx] = sample;
                    valuesAdjust[geneIdx] = sample  * alphaSEScale/ betaSEScale;
                    ssqGene += sample * sample;
                    betaCorr = betaCorr + H.col(i)*(oldSample - values[geneIdx]);
                    betaHat = betaHat + H.col(i) * values[geneIdx];
                    deltaGene[geneIdx] = 1;
                    nnGene = nnGene + 1;
                } else {
                    if (oldSample) {
                        betaCorr = betaCorr + H.col(i)*oldSample;
                    }
                    values[geneIdx] = 0;
                }
            } // End of gene loop
        } // End of mcmc sampling loop
        // double scaleGene = 1/2 * Gadget::calcVariance(betaHat);
        // double sse = (betaTotal(eQTLCommon) - H * values).dot(betaTotal(eQTLCommon) - H * values);
        // double dfTilde = 4 + nnGene;
        // double scaleTilde = sse + 4 * scaleGene;
        // vareMed = InvChiSq::sample(dfTilde, scaleTilde); 
        vareMed = Gadget::calcVariance(betaCorr);
        propMed = Gadget::calcVariance(betaHat)/(Gadget::calcVariance(betaHat) + vareMed);
    } else {
        for (unsigned i = 0; i < numGenes; i++) {
            values[i] = 0; // here we assign gene effect as zero since it's not available.
            valuesAdjust[i]  = 0;
        } 
    }
    // if(eQTLNumAcrossGenes < 9){
    //     cout << values << endl;
    // }
    // cout << "gene effect: " << values << endl;
    // outFile = "geneEffect_theta_" + to_string(iter) + ".txt";
    // file1.open((outPath + outFile).c_str());
    // file1 << values << endl;
    // file1.close();
}

void ApproxBayesCO::GeneEffects::sampleFromeFC(Data data, int iter ,const string title, VectorXd &betaTotal, EQTLJointVec &eQTLJointVec,
                            const map<int,vector<int> > &gene2gwasSnpMap,const map<int,string> &gwasSnpIdx2snpIDMap,
                            double &sigmaSqTheta, double &vareMed, const double &piTheta, VectorXd &deltaGene ){
    double sample, oldSample,varRes,rhs,uhat,invLhs,muHat;
    double logDelta0, logDelta1, probDelta1;
    double sigmaSqThetaInv, sigmaSqThetaLog;
    double muBeta;
    int eQTLNumAcrossGenes =0;
    map<int,vector<int>> gene2gwasSnpIdxMapLocal; // gwas-snp
    
    vector<int > eQTLCommonCis; // cis-snp
    map<int,int> gwasSnpIdx2HidxMap;
    gene2gwasSnpIdxMapLocal.clear();
    eQTLCommonCis.clear();
    deltaGene.setZero();

    VectorXd betaCorr, betaHat;
    MatrixXd H, HPH; 
    set<int> eQTLUniqSet; 
    vector<int> geneIdxVec;
    geneIdxVec.clear();
    bool geneHasSnps;
    //  debug
    std::ofstream file1, file2;
    string outPath = data.label;
    string outFile;
    gwasSnpIdx2HidxMap.clear();
    int Hidx = 0;
    geneEffectScaleFactorMean.setZero(numGenes);
    for(int i =0; i < numGenes; ++i){
        geneHasSnps = false;
        //////////////////////////////////////////////////////////////////
        // restrict SNPs in the genic region and find nonzero snp sets
        //////////////////////////////////////////////////////////////////
        vector<int> snpIdxInGeneSet = gene2gwasSnpMap.at(i);
        vector<int> eQTLCommon;
        // string outFile = "geneEffect_alpha_gene_"+ to_string(i) + "_" + to_string(iter) + ".txt";
        // file1.open((outPath + outFile).c_str());
        // file1 << eQTLJointMat(gene2cisSnpIdxPerGene,i) << endl;
        // file1.close();
        /// loop all eqtl within the gene to find nonzero eqtls
        for(int j = 0; j < eQTLJointVec[i]->values.size(); j++){
            int gwasIdx = snpIdxInGeneSet[j];
            string gwasID = gwasSnpIdx2snpIDMap.at(gwasIdx);
            if(eQTLUniqSet.find(gwasIdx) != eQTLUniqSet.end()){
                continue;
            }
            if (betaTotal[gwasIdx] != 0 && eQTLJointVec[i]->values(j) != 0 ){
                eQTLCommon.push_back(gwasIdx);
                eQTLUniqSet.insert(gwasIdx); 
                gwasSnpIdx2HidxMap.insert(pair<int,int>(gwasIdx,Hidx));
                geneHasSnps = true;
                Hidx++;
                // add scale factor
                double singleScale = data.scalingeQTLFactorVecVec[i][j]/ data.scalingGWASFactorVec[gwasIdx];
                geneEffectScaleFactorMean[i] = geneEffectScaleFactorMean[i] + singleScale;
            }
        }
        // if nonzero eqtl exists, add scale factor
        if(geneHasSnps){
            geneIdxVec.push_back(i);
            gene2gwasSnpIdxMapLocal.insert(pair<int,vector<int>>(i,eQTLCommon));
            geneEffectScaleFactorMean[i] = geneEffectScaleFactorMean[i]/eQTLCommon.size();
        }
    }
    //////////////////////////////////////////////////////////////////
    /// construct betaCorr vector and H matrix based on nonzero vector
    /////////////////////////////////////////////////////////////////
    H.setZero(eQTLUniqSet.size(),geneIdxVec.size());
    betaCorr.setZero(eQTLUniqSet.size());
    for (unsigned j = 0,k=0; j < geneIdxVec.size(); j++) {
        int geneIdx = geneIdxVec[j];
        vector<int> snpIdxInGeneSet = gene2gwasSnpMap.at(j);
        vector<int> eQTLCommon = gene2gwasSnpIdxMapLocal.at(geneIdx);
        double betaSingleScaleFactor, alphaSingleScaleFactor;
        for(unsigned i = 0; i < eQTLCommon.size(); i++){
            int gwasIdx = eQTLCommon[i];
            string gwasID = gwasSnpIdx2snpIDMap.at(gwasIdx);
            
            Hidx = gwasSnpIdx2HidxMap.at(gwasIdx);
            betaCorr(Hidx) = betaTotal[gwasIdx];
            H(Hidx,k) = eQTLJointVec[geneIdx]->getValue(gwasID);

            
            // calculate scale factor
            // double betaSEScale = Gadget::calcMean(data.scalingGWASFactorVec);

                // double alphaSEScale = Gadget::calcMean(data.scalingeQTLFactorVecVec[geneIdx]);
        }
        k++;
    }

    nnGene = 0;
    ssqGene = 0;
    eQTLNumAcrossGenes = eQTLUniqSet.size();
    sigmaSqThetaInv = 1/sigmaSqTheta;
    sigmaSqThetaLog = logf(sigmaSqTheta);

    // betaCorr = betaTotal(eQTLUniqSet);
    // H = eQTLJointMat(eQTLCommonCis,geneIdxVec);
    HPH = (H.transpose() * H).diagonal();
    // // gwas snplist
    // std::ofstream file1, file2;
    // outFile = "geneEffect_beta_" + to_string(iter) + ".txt";
    // file1.open((outPath + outFile).c_str());
    // file1 << betaCorr << endl;
    // file1.close();
    // outFile = "geneEffect_H_" + to_string(iter) + ".txt";
    // file1.open((outPath + outFile).c_str());
    // file1 << H << endl;
    // file1.close();

    // if(eQTLNumAcrossGenes < 9){
    //     cout << "here" << endl;
    // }
    if(eQTLNumAcrossGenes >= 2){
        muBeta = betaCorr.mean();
        betaCorr.array() = betaCorr.array() - muBeta - (H * values(geneIdxVec)).array();
        /// Start mcmc process
        for(unsigned iteri =0; iteri < 1; iteri++){
            betaCorr.array() += muBeta;
            rhs = betaCorr.sum();
            invLhs = 1/(double) eQTLNumAcrossGenes;
            muHat = invLhs*rhs;
            muBeta = Normal::sample(muHat,invLhs*vareMed);
            betaCorr.array() -= muBeta;
            betaHat.setZero(eQTLNumAcrossGenes); // betaHat = 0;
            for (unsigned i = 0; i < geneIdxVec.size(); i++) {
                int geneIdx = geneIdxVec[i];

                double logPi = logf(piTheta);
                double logPiComp = logf(1-piTheta);
                oldSample = values[geneIdx];
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
                    sample = uhat + Stat::snorm()*sqrtf(invLhs);
                    values[geneIdx] = sample;
                    valuesAdjust[geneIdx] = sample * geneEffectScaleFactorMean[geneIdx] ;
                    ssqGene += sample * sample;
                    betaCorr = betaCorr + H.col(i)*(oldSample - values[geneIdx]);
                    betaHat = betaHat + H.col(i) * values[geneIdx];
                    deltaGene[geneIdx] = 1;
                    nnGene = nnGene + 1;
                } else {
                    if (oldSample) {
                        betaCorr = betaCorr + H.col(i)*oldSample;
                    }
                    values[geneIdx] = 0;
                    valuesAdjust[geneIdx] = 0;
                }
            } // End of gene loop
        } // End of mcmc sampling loop
        // double scaleGene = 1/2 * Gadget::calcVariance(betaHat);
        // double sse = (betaTotal(eQTLCommon) - H * values).dot(betaTotal(eQTLCommon) - H * values);

        // double dfTilde = 4 + nnGene;
        // double scaleTilde = sse + 4 * scaleGene;
        // vareMed = InvChiSq::sample(dfTilde, scaleTilde); 
        vareMed = Gadget::calcVariance(betaCorr);
        propMed = Gadget::calcVariance(betaHat)/(Gadget::calcVariance(betaHat) + vareMed);
    } else {
        for (unsigned i = 0; i < numGenes; i++) {
            values[i] = 0; // here we assign gene effect as zero since it's not available.
            valuesAdjust[i] = 0;
        } 
    }
    // if(eQTLNumAcrossGenes < 9){
    //     cout << values << endl;
    // }
    // cout << "gene effect: " << values << endl;
    // outFile = "geneEffect_theta_" + to_string(iter) + ".txt";
    // file1.open((outPath + outFile).c_str());
    // file1 << values << endl;
    // file1.close();
}

void ApproxBayesCO::Rounding::computeRcorr(vector<MatrixXd> &LDBlockUs,vector<VectorXd> &LDBlockLambdas, const vector<ChromInfo*> &chromInfoVec, const VectorXd &snpEffects, VectorXd &rcorr){
    if (count++ % 100) return;
    VectorXd rcorrOld = rcorr;
    for (unsigned chr=0; chr<chromInfoVec.size(); ++chr) {
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        for (unsigned i=chrStart; i<=chrEnd; ++i) {
        }
    }
}

void ApproxBayesCO::SigmaSqTheta::sampleFromFC(const double snpEffSumSq, const unsigned numSnpEff){
    double dfTilde = df + numSnpEff;
    double scaleTilde = snpEffSumSq + df * scale;
    value = InvChiSq::sample(dfTilde, scaleTilde);
}

void ApproxBayesCO::setStartVal(void){
    if(data.numKeptGenes)
        sigmaSqMats.setPrior(sigmaSqBetaEqtl.value, sigmaSqAlphaVec.values);
        sigmaSqMatRes.setPrior(vare.value,varEps.values);
}

void ApproxBayesCO::sampleUnknowns() {
    static int iter = 0;
    unsigned cnt=0;
    if(mcmcType == "AIAO"){
        do {
            ////// Step 1. Sampling effect pair
            snpEffects.sampleFromFCAIAO(data,data.Qblocks,iter,diagnose,data.label,wcorrBlocks,wAcorr,wbcorrGene,data.QblocksDat,data.QgeneDat,data.ldblock2gwasSnpMap,
            snpEffectVec, eQTLJointVec,
            data.gwasSnpIdx2snpIDMap,data.geneID2IdxMap,data.gwasSnpID2geneIDMap,data.cisSnpID2IdxMap,sigmaSqBetaNonEqtl.value, sigmaSqMats,sigmaSqMatRes,
            piEffEqtl.value,piEffNonEqtl.value,data.n,data.neQTLVec,varEps.values,vare.valueVec);
            if (++cnt == 100) LOGGER.e(0," Zero SNP effect in the model for 100 cycles of sampling");
        } while (snpEffects.numNonZeros == 0 && snpEffects.numNonZerosEqtlVec(1) == 0);
        if(true){
            ////// Step 2.1 Sampling pi for effect pair
            // genic region
            if(iter == 0) LOGGER << "Sample piEffEqtl for genic region." << endl;
            if(data.numKeptGenes != 0){
                piEffEqtl.sampleFromFC(data.numEqtlOverlap, snpEffects.numNonZerosEqtlVec(0));
            }else{
                piEffEqtl.value = 0;
            }
        } // end of estimatePi
        if(data.numKeptGenes != 0){
            if(iter == 0) LOGGER << "Sample gene effect for genic region." << endl;
            sigmaSqTheta.sampleFromFC(geneEffectVec.values.dot(geneEffectVec.values),geneEffectVec.nnGene);
            piTheta.sampleFromFC(geneEffectVec.numGenes,geneEffectVec.nnGene);
            /////// Step 5. Sampling gene effect theta k
            // geneEffectVec.sampleFromeFC(data,iter,data.label,snpEffects.betaTotal,eQTLJointMat.values,data.gene2gwasSnpMap,data.gene2cisSnpMap,sigmaSqTheta.value,geneEffectVec.vareMed,piTheta.value, geneEffectVec.deltaGene);
            geneEffectVec.sampleFromeFC(data,iter,"",snpEffects.betaTotal,eQTLJointVec,data.gene2gwasSnpMap,data.gwasSnpIdx2snpIDMap,sigmaSqTheta.value,geneEffectVec.vareMed,piTheta.value, geneEffectVec.deltaGene);

        } // end of gene effect estimation
        if(true){
            /////// Step 3 Sampling variance covariance matrix 
            if(data.numKeptGenes != 0){
                if(iter == 0) LOGGER << "Sample genetic variance-covariance matrix for genic region." << endl;
                // sigmaSqBetaEqtl.value = sigmaSqBetaNonEqtl.value;
                // sigmaSqBetaEqtl.sampleFromFC(snpEffects.ssqBetaTotalGenic,snpEffects.numNonNullBetaTotGenic,snpEffects.vargGenic,piEffEqtl.value,iter,data.burnIn);
                // sigmaSqBetaEqtl.sampleFromFC(snpEffects.ssqBetaEqtl,snpEffects.numNonZerosEqtlVec(0),snpEffects.vargGenic,piEffEqtl.value,iter,data.burnIn);
                // sigmaSqAlphaVec.sampleFromFC(snpEffects.ssqAlphaEqtlPG,snpEffects.numNonZerosEqtlVecAcrossGenesPostIW,vargGeneCis.values,piEffEqtl.value,iter,data.burnIn);
                 // sigmaSqMats.sampleFromFCIWIndSMatPriorAIAO(iter,data.burnIn, geneEffectVec.values,snpEffects.ssqBetaEqtl, snpEffects.ssqAlphaEqtl,snpEffects.ssqBetaTotalGenic,snpEffects.numNonNullBetaTotGenic,snpEffects.ssqBetaEqtlPG,snpEffects.ssqAlphaEqtlPG, snpEffects.ssqEqtlMat,snpEffects.numNonZerosEqtlVec,snpEffects.numNonZerosEqtlVecAcrossGenesPostIW,cisHsq.valueGenicVec.sum(),cisHsq.values,message);
                sigmaSqMats.sampleFromFCIWSMatCorrPriorAIAO(iter,data.burnIn,sigmaSqBetaEqtl.value,sigmaSqAlphaVec.values, geneEffectVec.values,snpEffects.ssqBetaEqtl, snpEffects.ssqAlphaEqtl,snpEffects.ssqBetaTotalGenic,snpEffects.numNonNullBetaTotGenic,snpEffects.ssqBetaEqtlPG,snpEffects.ssqAlphaEqtlPG, snpEffects.ssqEqtlMat,snpEffects.numNonZerosEqtlVec,snpEffects.numNonZerosEqtlVecAcrossGenesPostIW,cisHsq.valueGenicVec.sum(),cisHsq.values.sum(),message);
            } // end of gene effect estimation
        }
    } else if (mcmcType == "EIEO"){
        ////// Step 1. Sampling effect pair
        //do {
            snpEffects.sampleFromFCEIEO(data,data.Qblocks,iter,diagnose,data.label,wcorrBlocks,wAcorr,data.QblocksDat,data.QgeneDat,data.ldblock2gwasSnpMap,
            snpEffectVec,eQTLJointVec,snpEffectVecLatent,eQTLJointVecLatent,
            data.gwasSnpIdx2snpIDMap,data.geneID2IdxMap,data.gwasSnpID2geneIDMap,data.cisSnpID2IdxMap,sigmaSqBetaNonEqtl.value, sigmaSqMats,
            piEffEieo1.value,piEffEieo2.value,piEffNonEqtl.value,data.n,data.neQTLVec,varEps.values,vare.valueVec);
            //if (++cnt == 100) LOGGER.e(0," Zero SNP effect in the model for 100 cycles of sampling");
        //} while (snpEffects.numNonZeros == 0 && snpEffects.numNonZerosEqtlVec(1) ==0);
        /////// Step 2.1 Sampling Pi for effect pair 
        if(true){
            if(data.numKeptGenes != 0){
                piEffEieo1.sampleFromFC(data.numEqtlOverlap, snpEffects.numNonZerosEqtlVec(0));
                piEffEieo2.sampleFromFC(data.numEqtlOverlap, snpEffects.numNonZerosEqtlVec(1));
            }else{
                piEffEieo1.value = 0;
                piEffEieo2.value = 0;
            }
        } /// end of pi sampling
        if(true){
            if(data.numKeptGenes != 0){
                ////// Step 5. Sampling gene effect 
                sigmaSqTheta.sampleFromFC(geneEffectVec.values.dot(geneEffectVec.values),geneEffectVec.nnGene);
                piTheta.sampleFromFC(geneEffectVec.numGenes,geneEffectVec.nnGene);
                // if(eieoLatent){
                //     geneEffectVec.sampleFromeFC(data,iter,"",snpEffects.betaTotalLatent,eQTLJointMatLatent.values,data.gene2gwasSnpMap,data.gene2cisSnpMap,sigmaSqTheta.value,geneEffectVec.vareMed,piTheta.value, geneEffectVec.deltaGene);
                // } else{
                //     geneEffectVec.sampleFromeFC(data,iter,"",snpEffects.betaTotal,eQTLJointMat.values,data.gene2gwasSnpMap,data.gene2cisSnpMap,sigmaSqTheta.value,geneEffectVec.vareMed,piTheta.value, geneEffectVec.deltaGene);
                // }
                if(eieoLatent){
                    geneEffectVec.sampleFromeFC(data,iter,"",snpEffects.betaTotalLatent,eQTLJointVecLatent,data.gene2gwasSnpMap,data.gwasSnpIdx2snpIDMap,sigmaSqTheta.value,geneEffectVec.vareMed,piTheta.value, geneEffectVec.deltaGene);
                } else{
                geneEffectVec.sampleFromeFC(data,iter,"",snpEffects.betaTotal,eQTLJointVec,data.gene2gwasSnpMap,data.gwasSnpIdx2snpIDMap,sigmaSqTheta.value,geneEffectVec.vareMed,piTheta.value, geneEffectVec.deltaGene);
                }
            } // end of gene effect estimation
        }
        // gene effec and genetic variance and covariance
        if(true){
            ////// Step 3. Sampling variance-covariance matrix 
            if(data.numKeptGenes != 0){
                // test it use independent prior;
                // sigmaSqBetaEqtl.sampleFromFC(snpEffects.ssqBetaEqtl,snpEffects.numNonZerosEqtlVec(0),snpEffects.vargGenic,piEffEieo1.value,iter,data.burnIn);
                // sigmaSqAlphaVec.sampleFromFC(snpEffects.ssqAlphaEqtlPG,snpEffects.numNonZerosGenicEqtlVec,vargGeneCis.values,piEffEieo2.value,iter,data.burnIn);
                // sigmaSqBetaEqtl.sampleFromFC(snpEffects.ssqBetaEqtl,data.numEqtl,snpEffects.vargGenic,piEffEieo1.value,iter,data.burnIn);
                // sigmaSqAlphaVec.sampleFromFC(snpEffects.ssqAlphaEqtlPG,sigmaSqAlphaVec.numEqtlPG,vargGeneCis.values,piEffEieo2.value,iter,data.burnIn);
                // sigmaSqMats.sampleFromFCIWSMatCorrPrior(iter,data.burnIn, geneEffectVec.values,snpEffects.ssqBetaEqtl, snpEffects.ssqAlphaEqtl,snpEffects.ssqBetaTotalGenic,snpEffects.numNonNullBetaTotGenic, snpEffects.ssqBetaEqtlPG,snpEffects.ssqAlphaEqtlPG, snpEffects.ssqEqtlMat,snpEffects.numNonZerosEqtlVec,snpEffects.numNonZerosEqtlVecAcrossGenesPostIW,cisHsq.valueGenicVec.sum(),cisHsq.values.sum(),message);
                geneticCorr.value = 0;
                if(iter == 0) LOGGER << "Now we fixed sigmaSq alpha and beta " << endl;
                sigmaSqMats.sampleFromFCIWIndSMatPriorEIEO(iter,data.burnIn,sigmaSqBetaEqtl.value,sigmaSqAlphaVec.values,snpEffects.ssqEqtlMat,geneticCorr.value,
                  geneEffectVec.values,sigmaSqAlphaVec.numEqtlPG,message);
                geneticCorr.value = geneticCorr.value/data.numGenes;

            } // end of gene effect estimation
        }
    } else {
        LOGGER.e(0,"Wrong mcmc type for SBayesCO model, please select AIAO or EIEO.");
    }
   
    ///////////////////////////////////////////////////
    ///// residual sampling
    ////// Step 4. Sampling residual variance.
    if(true){
        bool isSamVare = false;
        if(iter == 0) LOGGER << "Sample GWAS residual for whole-genome region." << endl;
        if(!sampleVareBool){
            vare.sampleFromFC(iter,data.varPhenotypic,wcorrBlocks,data.gwasEffectInBlock,data.QblocksDat, data.n, snpEffects.betaTotal, data.ldblock2gwasSnpMap);
        } else{
            vare.sampleFromFC(iter,wcorrBlocks, snpEffects.whatBlocks, snpEffects.ssqBlocks, data.nGWASblock, data.numEigenvalBlock);
        }

        if(data.numKeptGenes != 0) {
            if(iter == 0) LOGGER << "Sample xQTL residuals for genic region." << endl;
            // varEps.sampleFromFC(data.varPhenotypiceQTL,wAcorr,data.eQTLEffAcrossGenes,data.QgeneDat,data.neQTLVec,eQTLJointMat.values,data.gene2cisSnpMap);
            varEps.sampleFromFC(data.varPhenotypiceQTL,wAcorr,data.eQTLEffAcrossGenes,data.QgeneDat,data.neQTLVec,eQTLJointVec);
            // sigmaSqMatRes.sampleFromFC(wbcorrGene,wAcorr,varEps.values,vare.valueVec,data.ldblock2gwasSnpMap,data.geneID2IdxMap,data.gwasSnpIdx2snpIDMap,data.gwasSnpID2geneIDMap);
        }
    }
    vareMean.value = vare.valueVec.mean();
    if(data.numKeptGenes) varEpsMean.value = varEps.values.mean();
    //////////////////////////////////////////
    // pi for intergenic 
    // intergenic region
    ///// Step 2.2
    if(true){
        if(data.numNonEqtl != 0){
        piEffNonEqtl.sampleFromFC(data.numNonEqtl, snpEffects.numNonZerosNonEqtl);
        // piEffNonEqtl.sampleFromFC(data.numIncdSnps, snpEffects.numNonZerosNonEqtl + snpEffects.numNonNullBetaTotGenic);
        }else{
            piEffNonEqtl.value = 0;
        }
    }
    ////////////////////////////////////////////
    /////// make convergence of mcmc chain
    nBadSnps.compute_eigen(snpEffects.badSnps, snpEffects.values, snpEffects.posteriorMean, data.b, wcorrBlocks, data.Qblocks, data.keptLdBlockInfoVec, iter);
    ////////////////////////////////////////////
    /////// summary of various non-zero effects
    if(data.numKeptGenes != 0){
        if(mcmcType == "AIAO"){
            nnzGen.getValue(snpEffects.numNonZerosEqtlVec(0));
            nnsGen.getValue(snpEffects.numNonNullEqtl);
            nnEqtlOverlap.getValue(data.numEqtlOverlap);
        } else if (mcmcType == "EIEO"){
            nnEqtlOverlap.getValue(data.numEqtlOverlap);
            nnzGen.getValue(snpEffects.numNonZerosEqtlVec(0));
            nnsGen.getValue(snpEffects.numNonZerosEqtlVec(1));
            nnEqtl.getValue(snpEffects.numNonNullEqtl);
            nsnp00.getValue(snpEffects.numSnpCompVec(0));
            nsnp10.getValue(snpEffects.numSnpCompVec(1));
            nsnp01.getValue(snpEffects.numSnpCompVec(2));
            nsnp11.getValue(snpEffects.numSnpCompVec(3));
        }
    } 
    /// general
    nnsTot.getValue(snpEffects.numNonZeros);
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
    if(true){
        if(data.numNonEqtl != 0) {
            nnzBtw.getValue(snpEffects.numNonZerosNonEqtl);
            sigmaSqBetaNonEqtl.sampleFromFC(snpEffects.ssqNonEqtl, snpEffects.numNonZerosNonEqtl ,snpEffects.vargInt,piEffNonEqtl.value,iter,data.burnIn);
        } else {
        nnzBtw.value = 0;
        sigmaSqBetaNonEqtl.value = 0;
        }
    }
    // general summary
    ///////////////////////////////////////
    varg.compute(iter,snpEffects.betaTotal, data.Qblocks,data.ldblock2gwasSnpMap);
    vargGene.compute(geneEffectVec.values,eQTLJointVec, data.QgeneDat);
    geneEffectVec.values = geneEffectVec.valuesAdjust;
    vargGeneCis.compute(eQTLJointVec,data.QgeneDat);
    vargGeneCis.computeHsqBeta(snpEffectVec,data.QgeneDat);
    // varg.value = snpEffects.vargGenic + snpEffects.vargInt;
    hsq.value = varg.value / data.varPhenotypic;
    cisHsq.valueGenicVec = (vargGeneCis.valueBetaVec / data.varPhenotypic);
    
    if(data.numKeptGenes != 0){
        medHsq.compute(vargGene.value,data.varPhenotypic);
        cisHsq.compute(vargGeneCis.values, data.varPhenotypiceQTL);
        gwasCisHsq.compute(vargGeneCis.valueBetaVec,data.varPhenotypic);
        genicGwasEnrich.compute(gwasCisHsq.values,snpEffects.numNonZerosGenicGwasVec);
        genicEqtlEnrich.compute(cisHsq.values,snpEffects.numNonZerosGenicEqtlVec);
        sigmaSqAlpha.value = sigmaSqAlphaVec.values.mean();
        sigmaSqBetaEqtl.value = sigmaSqBetaEqtl.value;
        cisHsqMean.value = cisHsq.values.mean();
    } else {
        medHsq.value = 0;  
        sigmaSqAlpha.value = 0;
        sigmaSqBetaEqtl.value = 0;  
        cisHsqMean.value = 0; 
    }
    if(diagnose){
        string outPath = data.label;
        string outFile;
        // // gwas snplist
        std::ofstream file1, file2;
        outFile = "sum-all-parameter-cpp-" + to_string(iter) + ".txt";
        file1.open((outPath + outFile).c_str());
        file1 << iter << "\t" << nnsTot.value << "\t" << nnsGen.value << "\t" << nnGene.value << "\t" << nnsPG.value << "\t" 
            << sigmaSqBetaNonEqtl.value << "\t" << sigmaSqBetaEqtl.value << "\t" << sigmaSqAlpha.value << "\t" << hsq.value << "\t"
            << medHsq.value << "\t" << cisHsqMean.value << "\t" << vareMean.value  << "\t" 
            << varEpsMean.value
            <<  endl;
        file1.close();
    }
    ++iter;
}
