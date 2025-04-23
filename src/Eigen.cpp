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


#include "Data.hpp"
#include <sys/stat.h>
#include <dirent.h>

///////////////////////////////////////////////////////////////////////////////////////
////////    Step 1. perform eigen-decomposition for ld blocks                   ///////
///////////////////////////////////////////////////////////////////////////////////////

void Data::readLDBlockInfoFile(const string &ldBlockInfoFile){
    // Read bim file: recombination rate is defined between SNP i and SNP i-1
    ifstream in(ldBlockInfoFile.c_str());
    if (!in) LOGGER.e(0,"can not open the file [" + ldBlockInfoFile + "] to read.");
    LOGGER << "Reading ld block info from file [" + ldBlockInfoFile + "]." << endl;
    ldBlockInfoVec.clear();
    ldBlockInfoMap.clear();
    string header;
    string id;
    int  chr,start, stop;
    int idx = 0;
    getline(in, header);
    while (in >> id >> chr >> start >> stop) {
        LDBlockInfo *ld = new LDBlockInfo(idx++, id, chr);
        ld->startPos    = start;
        ld->endPos     = stop;
        ldBlockInfoVec.push_back(ld);
        if (ldBlockInfoMap.insert(pair<string, LDBlockInfo*>(id, ld)).second == false) {
            LOGGER.e(0,"Duplicate LD block ID found: \"" + id + "\".");
        }
    }
    in.close();
    numLDBlocks = (unsigned) ldBlockInfoVec.size();
    LOGGER << numLDBlocks << " GWAS LD Blocks to be included from [" + ldBlockInfoFile + "]." << endl;
}

void Data::eigenDecomposition( const MatrixXf &X, const float &prop, VectorXf &eigenValAdjusted, MatrixXf &eigenVecAdjusted, float &sumPosEigVal){
    // VectorXf cumsumNonNeg; // cumulative sums of non-negative values
    sumPosEigVal = 0.0;

    SelfAdjointEigenSolver<MatrixXf> eigensolver(X);
    VectorXf eigenVal = eigensolver.eigenvalues();
    MatrixXf eigenVec = eigensolver.eigenvectors();
    int revIdx = eigenVal.size();
    VectorXf cumsumNonNeg(revIdx);
    cumsumNonNeg.setZero();
    revIdx = revIdx -1;
    if(eigenVal(revIdx) < 0) cout << "Error, all eigenvector are negative" << endl;
    cumsumNonNeg(revIdx) = eigenVal(revIdx);
    sumPosEigVal = eigenVal(revIdx);
    revIdx = revIdx -1;
    
    while( eigenVal(revIdx) > 1e-10 ){
        sumPosEigVal = sumPosEigVal + eigenVal(revIdx);
        cumsumNonNeg(revIdx) = eigenVal(revIdx) + cumsumNonNeg(revIdx + 1);
        revIdx =revIdx - 1;
        if(revIdx < 0) break;
    }
    // cout << "revIdx: " << revIdx << endl;
    // cout << "size: " << eigenVal.size()  << " eigenVal: " << eigenVal << endl;
    // cout << "cumsumNoNeg: " << cumsumNonNeg << endl;
    cumsumNonNeg = cumsumNonNeg/sumPosEigVal;
    bool haveValue = false;
    // cout << "cumsumNonNeg: " << cumsumNonNeg << endl;
    // cout << "revIdx : " << revIdx << endl;
    // cout << "revIdx: "  << revIdx  << endl;
    for (revIdx = revIdx + 1; revIdx < eigenVal.size(); revIdx ++ ){
        // cout << "revIdx: " << revIdx << " cumsumNonNeg: " << cumsumNonNeg(revIdx) << endl;
        if(prop >= cumsumNonNeg(revIdx) ){
            revIdx = revIdx -1;
            haveValue = true;
            break;
        }
    }
    // cout << "revIdx : " << revIdx << endl;
    if(!haveValue) revIdx = eigenVal.size() - 1;
    // cout << "cumsumNonNeg: " << cumsumNonNeg.size() << endl;
    // cout << "cumsumNoNeg: " << cumsumNonNeg << endl;
    // cout << "revIdx: "  << revIdx  << endl;

    eigenVecAdjusted = eigenVec.rightCols(eigenVal.size() - revIdx);
    eigenValAdjusted = eigenVal.tail(eigenVal.size() - revIdx);
    // cout << "eigenValAdjusted size: " << eigenValAdjusted.size() << " eigenValue eventually: " << eigenValAdjusted << endl;
    //eigenvalueNum = eigenVal.size() - revIdx;
    // cout << endl;
}

MatrixXd Data::generateLDmatrixPerBlock(const string &bedFile, const vector<string> &snplists){
    int numSnpInRange = snplists.size();
    IndInfo *indi = NULL;
    SnpInfo *snpj = NULL;
    SnpInfo *snpk = NULL;

    snpj = snpInfoMap.at(snplists[0]); // start
    snpk = snpInfoMap.at(snplists[numSnpInRange - 1]); // end;
    unsigned start = snpj->index;
    unsigned end = snpk->index;
    
    if (numIncdSnps == 0) LOGGER.e(0,"No SNP is retained for analysis.");
    if (numKeptInds == 0) LOGGER.e(0,"No individual is retained for analysis.");
    // if (start >= numIncdSnps) LOGGER.e(0,"Specified a SNP range of " + snpRange + " but " + to_string(static_cast<long long>(numIncdSnps)) + " SNPs are included.");

    //////////////////////////////////////////////////////
    // Step 1. read in the genotypes of SNPs in the given range
    //////////////////////////////////////////////////////
    
    const int bedToGeno[4] = {2, -9, 1, 0};
    unsigned size = (numInds+3)>>2;
    
    MatrixXd ZP(numSnpInRange, numKeptInds);  // SNP x Ind
    VectorXd Dtmp;
    Dtmp.setZero(numSnpInRange);

    if (numKeptInds < 2) LOGGER.e(0, " Cannot calculate LD matrix with number of individuals < 2.");
    
    FILE *in1 = fopen(bedFile.c_str(), "rb");
    if (!in1) LOGGER.e(0, " can not open the file [" + bedFile + "] to read.");
    // cout << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
    char header[3];
    fread(header, sizeof(header), 1, in1);
    if (!in1 || header[0] != 0x6c || header[1] != 0x1b || header[2] != 0x01) {
        cerr << "Error: Incorrect first three bytes of bed file: " << bedFile << endl;
        exit(1);
    }

    int genoValue;
    unsigned i, j, k;
    unsigned incj, inck; // index of included SNP
    unsigned long long skipj = 0;
    unsigned nmiss;
    float mean;
    
    set<int> chromInRange;
    
    for (j = 0, incj = 0; j < numSnps; j++) {
        snpj = snpInfoVec[j];
        
        if (snpj->index < start || !snpj->included) {
            skipj += size;
            continue;
        }

        // check if snp exist in snplist in case gap situation
        auto it = std::find(snplists.begin(), snplists.end(), snpj->rsID);
        if (it == snplists.end()) {
            skipj += size;
            continue;
        } 
        
        if (skipj) fseek(in1, skipj, SEEK_CUR);
        skipj = 0;
        
        char *bedLineIn = new char[size];
        fread(bedLineIn, sizeof(char), size, in1);
        
        chromInRange.insert(snpj->chrom);
        
        mean = 0.0;
        nmiss = 0;
        
        for (i = 0; i < numInds; i++) {
            indi = indInfoVec[i];
            if (!indi->kept) continue;
            genoValue = bedToGeno[(bedLineIn[i>>2]>>((i&3)<<1))&3];
            ZP(incj, indi->index) = genoValue;
            if (genoValue == -9) ++nmiss;
            else mean += genoValue;
        }
        delete[] bedLineIn;
        
        // fill missing values with the mean
        snpj->sampleSize = numKeptInds-nmiss;
        mean /= float(snpj->sampleSize);
        if (nmiss) {
            for (i=0; i<numKeptInds; ++i) {
                if (ZP(incj, i) == -9) ZP(incj, i) = mean;
            }
        }
        
        // compute allele frequency
        snpj->af = 0.5f*mean;
        snp2pq[incj] = snpj->twopq = 2.0f*snpj->af*(1.0f-snpj->af);
        
        if (snp2pq[incj]==0) LOGGER.e(0, " " + snpj->rsID + " is a fixed SNP (MAF=0)!");
        
        Dtmp[incj] = Gadget::calcVariance(ZP.row(incj))*numKeptInds;
        ZP.row(incj) = (ZP.row(incj).array() - ZP.row(incj).mean())/sqrt(Dtmp[incj]);
        if (++incj == numSnpInRange) break;
    }
    
    fclose(in1);
    
    ZPZdiag = ZP.rowwise().squaredNorm();
    
    //////////////////////////////////////////////////////
    // Step 2. read in the bed file again to compute Z'Z
    //////////////////////////////////////////////////////
    
    MatrixXd denseZPZ;
    denseZPZ.setZero(numSnpInRange, numSnpInRange);
    VectorXd Zk(numKeptInds);
    Dtmp.setZero(numSnpInRange);
    
    FILE *in2 = fopen(bedFile.c_str(), "rb");
    fseek(in2, 3, SEEK_SET);
    unsigned long long skipk = 0;
    
    set<int>::iterator setend = chromInRange.end();

    if (numSkeletonSnps) {
        for (k = 0, inck = 0; k < numSnps; k++) {
            snpk = snpInfoVec[k];
            
            // if (!snpk->included) {
            if (snpk->index < start || !snpk->included) {
                skipk += size;
                continue;
            }

            // check if snp exist in snplist in case gap situation
            auto it = std::find(snplists.begin(), snplists.end(), snpk->rsID);
            if (it == snplists.end()) {
                skipk += size;
                continue;
            } 

            if (chromInRange.find(snpk->chrom) == setend && !snpk->skeleton) {
                skipk += size;
                ++inck;       // ensure the index is correct
                continue;
            }
            
            if (skipk) fseek(in2, skipk, SEEK_CUR);
            skipk = 0;
            
            char *bedLineIn = new char[size];
            fread(bedLineIn, sizeof(char), size, in2);
            
            mean = 0.0;
            nmiss = 0;
            
            for (i = 0; i < numInds; i++) {
                indi = indInfoVec[i];
                if (!indi->kept) continue;
                genoValue = bedToGeno[(bedLineIn[i>>2]>>((i&3)<<1))&3];
                Zk[indi->index] = genoValue;
                if (genoValue == -9) ++nmiss;   // missing genotype
                else mean += genoValue;
            }
            delete[] bedLineIn;
            
            // fill missing values with the mean
            snpk->sampleSize = numKeptInds-nmiss;
            mean /= float(snpk->sampleSize);
            if (nmiss) {
                for (i=0; i<numKeptInds; ++i) {
                    if (Zk[i] == -9) Zk[i] = mean;
                }
            }
            
            // compute allele frequency
            snpk->af = 0.5f*mean;
            snp2pq[inck] = snpk->twopq = 2.0f*snpk->af*(1.0f-snpk->af);
            
            if (snp2pq[inck]==0) LOGGER.e(0, " " + snpk->rsID + " is a fixed SNP (MAF=0)!");
            Dtmp[inck] = Gadget::calcVariance(Zk.row(inck))*numKeptInds;   // calculate  variances fro each snp;
            Zk = (Zk.array() - Zk.mean())/sqrt(Dtmp[inck]); // center and standardize genotype
            denseZPZ.col(inck) = ZP * Zk;

            //++inck;
            if (++inck == numSnpInRange) break;
        }
    }
    else {
        for (k = 0, inck = 0; k < numSnps; k++) {
            snpk = snpInfoVec[k];
            
            // if (!snpk->included) {
            if (snpk->index < start || !snpk->included) {
                skipk += size;
                continue;
            }

            // check if snp exist in snplist in case gap situation
            auto it = std::find(snplists.begin(), snplists.end(), snpk->rsID);
            if (it == snplists.end()) {
                skipk += size;
                continue;
            } 

            if (skipk) fseek(in2, skipk, SEEK_CUR);
            skipk = 0;
            
            char *bedLineIn = new char[size];
            fread(bedLineIn, sizeof(char), size, in2);
            
            mean = 0.0;
            nmiss = 0;
            
            for (i = 0; i < numInds; i++) {
                indi = indInfoVec[i];
                if (!indi->kept) continue;
                genoValue = bedToGeno[(bedLineIn[i>>2]>>((i&3)<<1))&3];
                Zk[indi->index] = genoValue;
                if (genoValue == -9) ++nmiss;   // missing genotype
                else mean += genoValue;
            }
            delete[] bedLineIn;
            
            // fill missing values with the mean
            snpk->sampleSize = numKeptInds-nmiss;
            mean /= float(snpk->sampleSize);
            if (nmiss) {
                for (i=0; i<numKeptInds; ++i) {
                    if (Zk[i] == -9) Zk[i] = mean;
                }
            }
            
            // compute allele frequency
            snpk->af = 0.5f*mean;
            snp2pq[inck] = snpk->twopq = 2.0f*snpk->af*(1.0f-snpk->af);
            
            if (snp2pq[inck]==0) LOGGER.e(0, " " + snpk->rsID + " is a fixed SNP (MAF=0)!");

            Dtmp[inck] = Gadget::calcVariance(Zk)*numKeptInds;

            Zk = (Zk.array() - Zk.mean())/sqrt(Dtmp[inck]);
            
            denseZPZ.col(inck) = ZP * Zk;

            if (++inck == numSnpInRange) break;
        }
    }

    fclose(in2);

    for (k = 0, inck = 0; k < numSnps; k++) {
        snpk = snpInfoVec[k];
        // if (!snpk->included) {
        if (snpk->index < start || !snpk->included) {
            skipk += size;
            continue;
        }
        if (skipk) fseek(in2, skipk, SEEK_CUR);
       skipk = 0;
    //    if(inck < 6){
    //     cout << snpk->rsID << " ";
    //    }
       inck ++; 
    }
    cout << endl;
    // LOGGER << "corr: " << denseZPZ.block(0, 0, 5, 5) << endl;
    return denseZPZ;
}

void Data::getEigenDataFromFullLDM(const string &filename, const float eigenCutoff){
    
    string outfilename = filename + ".eigen.bin";
    FILE *out3 = fopen(outfilename.c_str(), "wb");

    for (unsigned blk=0; blk < numKeptLDBlocks; blk++){
        
        MatrixXf eigenVec;
        VectorXf eigenVal;
        float sumPosEigVal;
        // cout << "rval: " << rval << endl;
        // cout << "rval cols: " << rval.cols() << " rval rows: " << rval.rows() << endl;
        eigenDecomposition(ZPZmat.cast<float>(), eigenCutoff, eigenVal, eigenVec, sumPosEigVal);
        // cout << "rval: " << rval.row(0) << endl;
        // cout << " Generate and save SVD of LD matrix from LD block " << i << "\r" << flush;
        // save svd matrix
        
        int32_t numEigenValue = eigenVal.size();
        int32_t numSnpInBlock = keptLdBlockInfoVec[blk]->numSnpInBlock;  // TMP

        cout << " Generate Eigen decomposition result for LD block " << blk << ", number of SNPs " << numSnpInBlock << ", number of selected eigenvalues " << numEigenValue << endl;
        
        // save summary
        // 1, the number of SNPs in the block
        fwrite(&numSnpInBlock, sizeof(int32_t), 1, out3);
        // 2, the number of eigenvalues at with the given cutoff
        fwrite(&numEigenValue, sizeof(int32_t), 1, out3);
        // 3. sum of all the positive eigenvalues
        fwrite(&sumPosEigVal, sizeof(float), 1, out3);
        // 4. eigenvalue cutoff based on the proportion of variance explained in LD
        fwrite(&eigenCutoff, sizeof(float), 1, out3);
        // 5. the selected eigenvalues
        fwrite(eigenVal.data(), sizeof(float), numEigenValue, out3);
        // 6. the selected eigenvector;
        uint64_t nElements = (uint64_t) numSnpInBlock * (uint64_t) numEigenValue;
        fwrite(eigenVec.data(), sizeof(float), nElements, out3);
        
    }
    
    fclose(out3);

    //outputEigenDataForLDM(filename);
    
}

void Data::mapSnpsToBlocks(int ldBlockRegionWind){
    int i,j;
    vector<locus_bp> snpVec;
    SnpInfo *snp;

    map<int, string>  chrEndSnp;
    for (i = 1; i < numIncdSnps; i++) {
        snp = incdSnpInfoVec[i];
        if(incdSnpInfoVec[i]->chrom != incdSnpInfoVec[i-1]->chrom){
            chrEndSnp.insert(pair<int, string>(incdSnpInfoVec[i - 1]->chrom,incdSnpInfoVec[i - 1]->rsID ));
        }
    }
    chrEndSnp.insert(pair<int, string>(incdSnpInfoVec[numIncdSnps - 1]->chrom,incdSnpInfoVec[numIncdSnps - 1]->rsID ));
    /////////////////////////////////////////
    // Step 2. Map snps to blocks
    /////////////////////////////////////////
    vector<string> block2snp_1(numLDBlocks), block2snp_2(numLDBlocks);
    map<string,int> keptLdBlock2AllLdBlcokMap;
    vector<locus_bp>::iterator iter;
    map<int, string>::iterator chrIter;
    LDBlockInfo *ldblock;
    for (i = 0; i < numIncdSnps ; i++) {
        snp = incdSnpInfoVec[i];
        snpVec.push_back(locus_bp(snp->rsID, snp->chrom, snp->physPos ));
    }

#pragma omp parallel for private(iter, chrIter)
    for (i = 0; i < numLDBlocks; i++) {
        // find lowest snp_name in the block
        ldblock = ldBlockInfoVec[i];

        iter = find_if(snpVec.begin(), snpVec.end(), locus_bp( ldblock->ID ,ldblock->chrom, ldblock->startPos - ldBlockRegionWind));
        if (iter != snpVec.end()) block2snp_1[i] = iter->locusName;
        else block2snp_1[i] = "NA";
    }
#pragma omp parallel for private(iter, chrIter)
    for (i = 0; i < numLDBlocks; i++) { 
        ldblock = ldBlockInfoVec[i];               
        if (block2snp_1[i] == "NA") {
            block2snp_2[i] = "NA";
            continue;
        }
        iter = find_if(snpVec.begin(), snpVec.end(), locus_bp(ldblock->ID, ldblock->chrom, ldblock->endPos + ldBlockRegionWind));
        if (iter != snpVec.end()){
            if (iter->bp ==  ldblock->endPos + ldBlockRegionWind){
                block2snp_2[i] = iter->locusName;
            }else {
                if(iter!=snpVec.begin()){
                    iter--;
                    block2snp_2[i] = iter->locusName;
                }
                else block2snp_2[i] = "NA";
            }
        }
        else {
            chrIter = chrEndSnp.find(ldblock->chrom);
            if (chrIter == chrEndSnp.end()) block2snp_2[i] = "NA";
            else block2snp_2[i] = chrIter->second;
        }
    }
    int mapped = 0;
    for (i = 0; i < numLDBlocks; i++) {
        ldblock = ldBlockInfoVec[i];
        if (block2snp_1[i] != "NA" && block2snp_2[i] != "NA") 
        {
            mapped++;
            // ldblock->kept = true;
            keptLdBlock2AllLdBlcokMap.insert(pair<string, int>(ldblock->ID, i));
        } else {
            ldblock->kept = false;
        }
    }
    if (mapped < 1) LOGGER.e(0, "No SNP can be mapped to the provided LD block list. Please check the input data regarding chromosome and bp.");
    else LOGGER << mapped << " GWAS LD block(s) are retained." << endl;

    keptLdBlockInfoVec = makeKeptLDBlockInfoVec(ldBlockInfoVec);
    numKeptLDBlocks = (unsigned) keptLdBlockInfoVec.size();

    string lastSNPIDInPreLDblocks;
    map<string, int>::iterator iter1, iter2;
    map<string, int> snpNameMap;
    
    for (i = 0; i < numIncdSnps; i++) {
        SnpInfo *snp = incdSnpInfoVec[i];
        snpNameMap.insert(pair<string,int>(snp->rsID, i));
        
    }
    for (i = 0; i < numKeptLDBlocks; i++) {
        ldblock = keptLdBlockInfoVec[i]; 
        ldblock->snpNameVec.clear();
        ldblock->snpInfoVec.clear();
        // LOGGER << "ldblock id: " << ldblock->ID << endl;
        iter1 = snpNameMap.find(block2snp_1[keptLdBlock2AllLdBlcokMap.at(ldblock->ID ) ]);
        iter2 = snpNameMap.find(block2snp_2[keptLdBlock2AllLdBlcokMap.at(ldblock->ID ) ]);
        bool skip = false;
        if (iter1 == snpNameMap.end() || iter2 == snpNameMap.end() || iter1->second >= iter2->second) ldblock->kept = false;
        if(!ldblock->kept) continue;
        vector<int> snp_indx;
        for (j = iter1->second; j <= iter2->second; j++) {
            snp = incdSnpInfoVec[j];
            if(snp->rsID == lastSNPIDInPreLDblocks){continue;} // if snp is matched to two ldblock
            ldblock->snpNameVec.push_back(snp->rsID);
            ldblock->snpInfoVec.push_back(snp);
            snp->isInBlock = true;
            snp->block = ldblock->ID;
            if (j == iter2->second) { lastSNPIDInPreLDblocks = snp->rsID;}
        }

        ldblock->numSnpInBlock = ldblock->snpInfoVec.size();
        ldblock->startSnpIdx = ldblock->snpInfoVec[0]->index;
        ldblock->endSnpIdx = ldblock->snpInfoVec[ldblock->numSnpInBlock-1]->index;
    }
 
}

void Data::makeBlockLDmatrix(const string &bedFile, const string &LDmatType, const unsigned block, const string &dirname, const bool writeLdmTxt, int ldBlockRegionWind){
    LOGGER << "Making block LD matricies ..." << endl;
    
    struct stat sb;
    if (stat(dirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
        // Folder doesn't exist, create it
        string create_cmd = "mkdir " + dirname;
        system(create_cmd.c_str());
        LOGGER << "Created folder [" << dirname << "] to store LD matrices." << endl;
    }
    
    mapSnpsToBlocks();
    
    keptLdBlockInfoVec = makeKeptLDBlockInfoVec(ldBlockInfoVec);
    numKeptLDBlocks = (unsigned) keptLdBlockInfoVec.size();
    
#pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < numKeptLDBlocks; i++) {
        LDBlockInfo *ldblock = keptLdBlockInfoVec[i];

        string outBinfile = dirname + "/block" + ldblock->ID + ".ldm.bin";
        FILE *outbin = fopen(outBinfile.c_str(), "wb");
        ofstream outtxt;
        string outTxtfile;
        if (writeLdmTxt) {
            outTxtfile = dirname + "/block" + ldblock->ID + ".ldm.txt";
            outtxt.open(outTxtfile.c_str());
        }
        string outSnpfile = dirname + "/block" + ldblock->ID + ".snp.info";
        string outldmfile = dirname + "/block" + ldblock->ID + ".ldm.info";

        if(!i) LOGGER << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
            
        MatrixXf rval = generateLDmatrixPerBlock(bedFile, ldblock->snpNameVec).cast<float>();
        
        unsigned numSnpInBlock = ldblock->numSnpInBlock;
        uint64_t nElements = (uint64_t) numSnpInBlock * (uint64_t) numSnpInBlock;
        fwrite(rval.data(), sizeof(double), nElements, outbin);
        
        if (writeLdmTxt) {
            for (unsigned ii=0; ii<numSnpInBlock; ++ii){
                for (unsigned jj=0; jj<numSnpInBlock; ++jj) {
                    outtxt << ldblock->ID << "\t" << ldblock->snpNameVec[ii] << "\t" << ldblock->snpNameVec[jj] << "\t" << rval(ii,jj) << endl;
                }
            }
        }
        
        fclose(outbin);
        if (writeLdmTxt) outtxt.close();
        
        outputBlockLDmatrixInfo(*ldblock, outSnpfile, outldmfile);
        
        if(!(i%1)) LOGGER << " computed block " << ldblock->ID << "\r" << flush;

        if (block) {
            LOGGER << "Written the LD matrix into file [" << outBinfile << "]." << endl;
            if (writeLdmTxt) LOGGER << "Written the LD matrix into file [" << outTxtfile << "]." << endl;
            LOGGER << "Written the LD matrix SNP info into file [" << outSnpfile << "]." << endl;
            LOGGER << "Written the LD matrix ldm info into file [" << outldmfile << "]." << endl;
        }
    }
    
    if (!block) {
        LOGGER << "Written the LD matrix into folder [" << dirname << "/block*.ldm.bin]." << endl;
        if (writeLdmTxt) LOGGER << "Written the LD matrix into text file [" << dirname << "/block*.ldm.txt]." << endl;
        
        if (chromInfoVec.size() >= 22) {  // genome-wide build of LD matrices
            mergeLdmInfo(LDmatType, dirname);
        } else {
            LOGGER << "Written the LD matrix into folder [" << dirname << "/block*.snp.info]." << endl;
            LOGGER << "Written the LD matrix into folder [" << dirname << "/block*.ldm.info]." << endl;
        }
    }
}

void Data::outputBlockLDmatrixInfo(const LDBlockInfo &block, const string &outSnpfile, const string &outldmfile) const {
    // write snp info
    ofstream out1(outSnpfile.c_str());
    out1 << boost::format("%6s %15s %10s %10s %15s %6s %6s %12s %10s %10s\n")
    % "Chrom"
    % "ID"
    % "Index"
    % "GenPos"
    % "PhysPos"
    % "A1"
    % "A2"
    % "A1Freq"
    % "N"
    % "Block";
    SnpInfo *snp;
    for (unsigned i=0; i < block.numSnpInBlock; ++i) {
        snp = block.snpInfoVec[i];
        out1 << boost::format("%6s %15s %10s %10s %15s %6s %6s %12f %10s %10s\n")
        % snp->chrom
        % snp->rsID
        % snp->index
        % snp->genPos
        % snp->physPos
        % snp->a1
        % snp->a2
        % snp->af
        % numKeptInds
        % snp->block;
    }
    out1.close();

    // svd matrix for ld blocks here.
    ofstream out2(outldmfile.c_str());
    out2 << boost::format("%10s %6s %15s %15s %15s %15s %12s\n")
    % "Block"
    % "Chrom"
    % "StartSnpIdx"
    % "StartSnpID"
    % "EndSnpIdx"
    % "EndSnpID"
    % "NumSnps";
    out2 << boost::format("%10s %6s %15s %15s %15s %15s %12s\n")
    % block.ID
    % block.chrom
    % block.startSnpIdx
    % incdSnpInfoVec[block.startSnpIdx]->rsID
    % block.endSnpIdx
    % incdSnpInfoVec[block.endSnpIdx]->rsID
    % block.numSnpInBlock;
    out2.close();
}

void Data::impG(const unsigned block, double diag_mod){
    VectorXi numImpSnp;
    VectorXi numTypSnp;
    numImpSnp.setZero(numLDBlocks);
    numTypSnp.setZero(numLDBlocks);
    for (unsigned i = 0; i < numLDBlocks; i++ ){
        LDBlockInfo *ldblock = ldBlockInfoVec[i];
        if (!ldblock->kept) continue;
        for (unsigned j=0; j<ldblock->numSnpInBlock; ++j) {
            SnpInfo *snp = ldblock->snpInfoVec[j];
            if (snp->included) {
                ++numTypSnp[i];
            } else {
                ++numImpSnp[i];
            }
        }
    }
    unsigned totalNumImpSnp = numImpSnp.sum();
        
    LOGGER << "Imputing summary statistics for " << to_string(totalNumImpSnp) << " SNPs in the LD reference but not in the GWAS data file..." << endl;

    Gadget::Timer timer;
    timer.setTime();
    
#pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < numLDBlocks; i++ ){
        LDBlockInfo *ldblock = ldBlockInfoVec[i];
        if (!ldblock->kept) continue;
        
        if (numImpSnp[i]) {
            
            Stat::Normal normal;
            
            /// Step 1. construct LD
            MatrixXd LDPerBlock = eigenVecLdBlock[i] * eigenValLdBlock[i].asDiagonal() * eigenVecLdBlock[i].transpose();
            
            LDPerBlock.diagonal().array() += (double)diag_mod;
            /// Step 2. Construct the LD correlation matrix among the typed SNPs(LDtt) and the LD correlation matrix among the missing SNPs and typed SNPs (LDit).
            // Step 2.1 divide SNPs into typed and untyped SNPs
            VectorXi typedSnpIdx(numTypSnp[i]);
            VectorXi untypedSnpIdx(numImpSnp[i]);
            VectorXd zTypSnp(numTypSnp[i]);
            VectorXd nTypSnp(numTypSnp[i]);
            VectorXd varyTypSnp(numTypSnp[i]);
            for(unsigned j=0, idxTyp=0, idxImp=0; j < ldblock->numSnpInBlock; j++){
                SnpInfo *snp = ldblock->snpInfoVec[j];
                if(snp->included){
                    // typed snp
                    typedSnpIdx[idxTyp] = j;
                    zTypSnp[idxTyp] = snp->gwas_b / snp->gwas_se;
                    nTypSnp[idxTyp] = snp->gwas_n;
                    double hetj = 2.0 * snp->gwas_af * (1.0 - snp->gwas_af);
                    varyTypSnp[idxTyp] = hetj * (snp->gwas_n * snp->gwas_se * snp->gwas_se + snp->gwas_b * snp->gwas_b);
                    ++idxTyp;
                } else {
                    untypedSnpIdx[idxImp] = j;
                    ++idxImp;
                }
            }
            // Step 2.2 construct LDtt and LDit and Ztt.
            MatrixXd LDtt = LDPerBlock(typedSnpIdx,typedSnpIdx);
            MatrixXd LDit = LDPerBlock(untypedSnpIdx,typedSnpIdx);
            // Step 2.3 //  The Z score for the missing SNPs;
            VectorXd LDi_Z = LDtt.ldlt().solve(zTypSnp);
            VectorXd zImpSnp = LDit * LDi_Z;
            // Step 3. re-calcualte beta and se
            // if snp is missing use median to replace N
            std::sort(nTypSnp.data(), nTypSnp.data() + nTypSnp.size());
            double nMedian = nTypSnp[nTypSnp.size()/2];  // median
            // calcuate median of phenotypic variance
            std::sort(varyTypSnp.data(), varyTypSnp.data() + varyTypSnp.size());
            double varyMedian = varyTypSnp[varyTypSnp.size()/2];  // median
            // begin impute
            for(unsigned j = 0; j < numImpSnp[i]; j++){
                SnpInfo *snp = ldblock->snpInfoVec[untypedSnpIdx[j]];
                double base = sqrt(2.0 * snp->af *(1.0 - snp->af) * (nMedian + zImpSnp[j] * zImpSnp[j]));
                snp->gwas_b = zImpSnp[j] * sqrt(varyMedian)/base;
                snp->gwas_se = sqrt(varyMedian) / base;
                snp->gwas_n = nMedian;
                snp->gwas_af = snp->af;
                snp->gwas_pvalue = 2*(1.0-normal.cdf_01(abs(snp->gwas_b/snp->gwas_se)));
                snp->included = true;
                //LOGGER << "b " << snp->gwas_b << " se " << snp->gwas_se << " z " << snp->gwas_b/snp->gwas_se << " p " << snp->gwas_pvalue << endl;
            }
        }
        
        if(!(i%10)) LOGGER << " imputed block " << i << "\r" << flush;
        
        if (block) {
            string outfile = title + ".block" + ldblock->ID + ".imputed.ma";
            ofstream out(outfile.c_str());
            out << boost::format("%15s %10s %10s %15s %15s %15s %15s %15s\n") % "SNP" % "A1" % "A2" % "freq" % "b" % "se" % "p" % "N";
            for (unsigned i=0; i<ldblock->numSnpInBlock; ++i) {
                SnpInfo *snp = ldblock->snpInfoVec[i];
                out << boost::format("%15s %10s %10s %15s %15s %15s %15s %15s\n")
                % snp->rsID
                % snp->a1
                % snp->a2
                % snp->gwas_af
                % snp->gwas_b
                % snp->gwas_se
                % snp->gwas_pvalue
                % snp->gwas_n;
            }
            out.close();

            timer.getTime();
            LOGGER << "Imputation of summary statistics is completed (time used: " << timer.format(timer.getElapse()) << ")." << endl;
            LOGGER << "Summary statistics of all SNPs are save into file [" + outfile + "]." << endl;

        }

    }

    if (block) return;
    
    string outfile = title + ".imputed.ma";
    ofstream out(outfile.c_str());
    out << boost::format("%15s %10s %10s %15s %15s %15s %15s %15s\n") % "SNP" % "A1" % "A2" % "freq" % "b" % "se" % "p" % "N";
    for (unsigned i=0; i<numSnps; ++i) {
        SnpInfo *snp = snpInfoVec[i];
        out << boost::format("%15s %10s %10s %15s %15s %15s %15s %15s\n")
        % snp->rsID
        % snp->a1
        % snp->a2
        % snp->gwas_af
        % snp->gwas_b
        % snp->gwas_se
        % snp->gwas_pvalue
        % snp->gwas_n;
    }
    out.close();

    timer.getTime();
    LOGGER << "Imputation of summary statistics is completed (time used: " << timer.format(timer.getElapse()) << ")." << endl;
    LOGGER << "Summary statistics of all SNPs are save into file [" + outfile + "]." << endl;
}

void Data::getEigenDataForLDBlock(const string &bedFile, const string &ldBlockInfoFile, int ldBlockRegionWind, const string &filename, const float eigenCutoff){
    int i,j;
    vector<locus_bp> snpVec;
    SnpInfo *snp;

    map<int, string>  chrEndSnp;
    for (i = 1; i < numIncdSnps; i++) {
        snp = incdSnpInfoVec[i];
        if(incdSnpInfoVec[i]->chrom != incdSnpInfoVec[i-1]->chrom){
            chrEndSnp.insert(pair<int, string>(incdSnpInfoVec[i - 1]->chrom,incdSnpInfoVec[i - 1]->rsID ));
        }
    }
    chrEndSnp.insert(pair<int, string>(incdSnpInfoVec[numIncdSnps - 1]->chrom,incdSnpInfoVec[numIncdSnps - 1]->rsID ));
    /////////////////////////////////////////
    // Step 2. Map snps to blocks
    /////////////////////////////////////////
    vector<string> block2snp_1(numLDBlocks), block2snp_2(numLDBlocks);
    map<string,int> keptLdBlock2AllLdBlcokMap;
    vector<locus_bp>::iterator iter;
    map<int, string>::iterator chrIter;
    for (i = 0; i < numIncdSnps ; i++) {
        snp = incdSnpInfoVec[i];
        snpVec.push_back(locus_bp(snp->rsID, snp->chrom, snp->physPos ));
    }
#pragma omp parallel for private(iter, chrIter)
    for (i = 0; i < numLDBlocks; i++) {
        // find lowest snp_name in the block
        LDBlockInfo *ldblock = ldBlockInfoVec[i];

        iter = find_if(snpVec.begin(), snpVec.end(), locus_bp( ldblock->ID ,ldblock->chrom, ldblock->startPos - ldBlockRegionWind));
        if (iter != snpVec.end()) block2snp_1[i] = iter->locusName;
        else block2snp_1[i] = "NA";
    }
#pragma omp parallel for private(iter, chrIter)
    for (i = 0; i < numLDBlocks; i++) {
        LDBlockInfo *ldblock = ldBlockInfoVec[i];
        if (block2snp_1[i] == "NA") {
            block2snp_2[i] = "NA";
            continue;
        }
        iter = find_if(snpVec.begin(), snpVec.end(), locus_bp(ldblock->ID, ldblock->chrom, ldblock->endPos + ldBlockRegionWind));
        if (iter != snpVec.end()){
            if (iter->bp ==  ldblock->endPos + ldBlockRegionWind){
                block2snp_2[i] = iter->locusName;
            }else {
                if(iter!=snpVec.begin()){
                    iter--;
                    block2snp_2[i] = iter->locusName;
                }
                else block2snp_2[i] = "NA";
            }
        }
        else {
            chrIter = chrEndSnp.find(ldblock->chrom);
            if (chrIter == chrEndSnp.end()) block2snp_2[i] = "NA";
            else block2snp_2[i] = chrIter->second;
        }
    }
    int mapped = 0;
    for (i = 0; i < numLDBlocks; i++) {
        LDBlockInfo *ldblock = ldBlockInfoVec[i];
        if (block2snp_1[i] != "NA" && block2snp_2[i] != "NA")
        {
            mapped++;
            // ldblock->kept = true;
            keptLdBlock2AllLdBlcokMap.insert(pair<string, int>(ldblock->ID, i));
        } else {
            ldblock->kept = false;
        }
    }
    if (mapped < 1) LOGGER.e(0, "No SNP can be mapped to the provided LD block list. Please check the input data regarding chromosome and bp.");
    else cout << mapped << " LD blocks have at least one SNP." << endl;

    keptLdBlockInfoVec = makeKeptLDBlockInfoVec(ldBlockInfoVec);
    numKeptLDBlocks = (unsigned) keptLdBlockInfoVec.size();

    map<string, int>::iterator iter1, iter2;
    map<string, int> snpNameMap;
    VectorXf snpNumInldblock(numKeptLDBlocks);
    VectorXf eigenNumInldblock(numKeptLDBlocks);
    vector<VectorXf> cumsumNonNeg(numKeptLDBlocks);
    
    for (i = 0; i < numIncdSnps; i++) {
        SnpInfo *snp = incdSnpInfoVec[i];
        snpNameMap.insert(pair<string,int>(snp->rsID, i));
        
    }
      // eigenvector and eigenvalue
    //eigenValLdBlock.resize(numKeptLDBlocks);
    //eigenVecLdBlock.resize(numKeptLDBlocks);
    string lastSNPIDInPreLDblocks;
    string LDmatType = "block";
    string outfilename = filename + "." + LDmatType + ".eigen" + ".bin";
    FILE *out3 = fopen(outfilename.c_str(), "wb");

    bool readBedBool = true;
    for (i = 0; i < numKeptLDBlocks; i++) {
        LDBlockInfo *ldblock = keptLdBlockInfoVec[i];
        // cout << "ldblock id: " << ldblock->ID << endl;
        iter1 = snpNameMap.find(block2snp_1[keptLdBlock2AllLdBlcokMap.at(ldblock->ID ) ]);
        iter2 = snpNameMap.find(block2snp_2[keptLdBlock2AllLdBlcokMap.at(ldblock->ID ) ]);
        bool skip = false;
        if (iter1 == snpNameMap.end() || iter2 == snpNameMap.end() || iter1->second >= iter2->second) ldblock->kept = false;
        snpNumInldblock[i] = iter2->second - iter1->second + 1;
        // cout << "ldblock->kept: " << ldblock->kept << endl;
        if(!ldblock->kept) continue;
        vector<int> snp_indx;
        for (j = iter1->second; j <= iter2->second; j++) {
            if(snp->rsID == lastSNPIDInPreLDblocks){continue;} // if snp is matched to two ldblock
            snp_indx.push_back(j);
            ldblock->snpNameVec.push_back(incdSnpInfoVec[j]->rsID);
            if (j == iter2->second) { lastSNPIDInPreLDblocks = snp->rsID;}
        }
        if(readBedBool) {
            cout << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
            readBedBool = false;
        }
        MatrixXf eigenVec;
        VectorXf eigenVal;
        float sumPosEigVal = 0;
        MatrixXf rval = generateLDmatrixPerBlock(bedFile, ldblock->snpNameVec).cast<float>();
        // cout << "rval: " << rval << endl;
        // cout << "rval cols: " << rval.cols() << " rval rows: " << rval.rows() << endl;
        eigenDecomposition(rval, eigenCutoff,eigenVal, eigenVec, sumPosEigVal);
        // cout << "rval: " << rval.row(0) << endl;
        // cout << " Generate and save SVD of LD matrix from LD block " << i << "\r" << flush;
        // save svd matrix
        int32_t numEigenValue = eigenVal.size();
        int32_t numSnpInBlock = ldblock->snpNameVec.size();
        eigenNumInldblock[i] = numEigenValue;
        // save summary
        // 1, nrow of eigenVecGene[i]
        fwrite(&numSnpInBlock, sizeof(int32_t), 1, out3);
        //        cout << "rval: " << rval << endl;
        // cout << "eigenVec: "  << eigenVec << endl;
        // cout << "";
        // 2, ncol of eigenVecGene[i]
        fwrite(&numEigenValue, sizeof(int32_t), 1, out3);
        // 3. sum of all positive eigenvalues
        fwrite(&sumPosEigVal, sizeof(float), 1, out3);
        // 4. eigenCutoff
        fwrite(&eigenCutoff, sizeof(float), 1, out3);
        // 5. eigenvalues
        fwrite(eigenVal.data(), sizeof(float), numEigenValue, out3);
        // 6. eigenvectors;
        uint64_t nElements = (uint64_t) numSnpInBlock * (uint64_t) numEigenValue;
        fwrite(eigenVec.data(), sizeof(float), nElements, out3);
        cout << " Generate and save Eigen decomposition result for LD block " << i << ", number of SNPs " << numSnpInBlock << ", number of selected eigenvalues " << numEigenValue << "\r" << flush;
    }
    fclose(out3);
    cout << "To explain " << eigenCutoff*100 << "% variance in LD, on average " << int(eigenNumInldblock.mean()) << " eigenvalues are selected across LD blocks (mean number of SNPs is " << int(snpNumInldblock.mean()) << ")." << endl;

}

void Data::readBlockLdmBinaryAndDoEigenDecomposition(const string &dirname, const unsigned block, const float eigenCutoff, const bool writeLdmTxt){
    struct stat sb;
    if (stat(dirname.c_str(), &sb) != 0 || !S_ISDIR(sb.st_mode)) {
        // Folder doesn't exist, create it
        LOGGER.e(0,"Error: cannot find the folder [" + dirname + "]");
    }
    
    vector<int> numSnpInRegion;
    numSnpInRegion.resize(numLDBlocks);
    for(int i = 0; i < numLDBlocks;i++){
        LDBlockInfo *block = ldBlockInfoVec[i];
        numSnpInRegion[i] = block->numSnpInBlock;
    }
        
    keptLdBlockInfoVec = ldBlockInfoVec;
    
#pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < numLDBlocks; i++){
        if (block && i != block - 1) continue;
        
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[i];

        string infile = dirname + "/block" + blockInfo->ID + ".ldm.bin";
        FILE *fp = fopen(infile.c_str(), "rb");
        if(!fp){LOGGER.e(0,"Error: can not open the file [" + infile + "] to read.");}

        string outBinfile = dirname + "/block" + blockInfo->ID + ".eigen.bin";
        FILE *outbin = fopen(outBinfile.c_str(), "wb");

        string outTxtfile;
        ofstream outtxt;
        if (writeLdmTxt) {
            outTxtfile = dirname + "/block" + blockInfo->ID + ".eigen.txt";
            outtxt.open(outTxtfile.c_str());
        }
        
        int32_t blockSize = blockInfo->numSnpInBlock;
        
        MatrixXf ldm(blockSize, blockSize);
        uint64_t nElements = (uint64_t)blockSize * (uint64_t)blockSize;
                
        if(fread(ldm.data(), sizeof(float), nElements, fp) != nElements){
            cout << "fread(U.data(), sizeof(float), nElements, fp): " << fread(ldm.data(), sizeof(float), nElements, fp) << endl;
            cout << "nEle: " << nElements << " ldm.size: " << ldm.size() <<  " ldm.col: " << ldm.cols() << " row: " << ldm.rows() << endl;
            LOGGER.e(0,"In LD block " + blockInfo->ID + ",size error in " + outBinfile);
            // cout << "Read " << svdLDfile << " error (U)" << endl;
            // LOGGER.e(0,"read file error");
        }
        
        MatrixXf eigenVec;
        VectorXf eigenVal;
        float sumPosEigVal;  // sum of all positive eigenvalues
        // cout << "rval: " << rval << endl;
        // cout << "rval cols: " << rval.cols() << " rval rows: " << rval.rows() << endl;
        eigenDecomposition(ldm, eigenCutoff, eigenVal, eigenVec, sumPosEigVal);
        // cout << "rval: " << rval.row(0) << endl;
        // cout << " Generate and save SVD of LD matrix from LD block " << i << "\r" << flush;
        // save svd matrix

        blockInfo->sumPosEigVal = sumPosEigVal;
        
        int32_t numEigenValue = eigenVal.size();
        int32_t numSnpInBlock = blockSize;

        //cout << " Generate Eigen decomposition result for LD block " << i << ", number of SNPs " << numSnpInBlock << ", number of selected eigenvalues " << numEigenValue << endl;
        
        // save summary
        // 1, the number of SNPs in the block
        fwrite(&numSnpInBlock, sizeof(int32_t), 1, outbin);
        // 2, the number of eigenvalues at with the given cutoff
        fwrite(&numEigenValue, sizeof(int32_t), 1, outbin);
        // 3. sum of all the positive eigenvalues
        fwrite(&sumPosEigVal, sizeof(float), 1, outbin);
        // 4. eigenvalue cutoff based on the proportion of variance explained in LD
        fwrite(&eigenCutoff, sizeof(float), 1, outbin);
        // 5. the selected eigenvalues
        fwrite(eigenVal.data(), sizeof(float), numEigenValue, outbin);
        // 6. the selected eigenvector;
        nElements = (uint64_t) numSnpInBlock * (uint64_t) numEigenValue;
        fwrite(eigenVec.data(), sizeof(float), nElements, outbin);
        
        if (writeLdmTxt) {
            outtxt << "Block " << blockInfo->ID << endl;
            outtxt << "numSnps " << numSnpInBlock << endl;
            outtxt << "numEigenvalues " << numEigenValue << endl;
            outtxt << "SumPositiveEigenvalues " << sumPosEigVal << endl;
            outtxt << "EigenCutoff " << eigenCutoff << endl;
            outtxt << "Eigenvalues\n" << eigenVal.transpose() << endl;
            outtxt << "Eigenvectors\n" << eigenVec << endl;
            outtxt << endl;
        }
        
        fclose(outbin);
        if (writeLdmTxt) outtxt.close();

        if(!(i%1)) cout << " computed block " << blockInfo->ID << "\r" << flush;
        
        if (block) {
            cout << "Written the eigen data for block LD matrix into file [" << outBinfile << "]." << endl;
            if (writeLdmTxt) cout << "Written the eigen data for block LD matrix into file [" << outTxtfile << "]." << endl;
        }
    }

    if (!block) {
        cout << "Written the eigen data for block LD matrix into file [" << dirname << "/block*.eigen.bin]." << endl;
        if (writeLdmTxt) cout << "Written the eigen data for block LD matrix into file [" << dirname << "/block*.eigen.txt]." << endl;
    }

}


void Data::readEigenMatrixBinaryFile(const string &dirname, const double eigenCutoff){
    if (!Gadget::directoryExist(dirname)) {
        LOGGER.e(0," cannot find the folder [" + dirname + "]");
    }
    
    vector<int>numSnpInRegion;
    LDBlockInfo * block;
    numSnpInRegion.resize(numLDBlocks);
    for(int i = 0; i < numLDBlocks;i++){
        block = ldBlockInfoVec[i];
        numSnpInRegion[i] = block->numSnpInBlock;
    }
    eigenValLdBlock.resize(numLDBlocks);
    eigenVecLdBlock.resize(numLDBlocks);
        
#pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < numLDBlocks; i++){
        LDBlockInfo * block;
        block = ldBlockInfoVec[i];
        
        if (!block->kept) continue;
        
        int32_t cur_m = 0;
        int32_t cur_k = 0;
        float sumPosEigVal = 0;
        float oldEigenCutoff =0;
        
        string infile = dirname + "/block" + block->ID + ".eigen.bin";
        FILE *fp = fopen(infile.c_str(), "rb");
        if(!fp){LOGGER.e(0,"can not open the file [" + infile + "] to read.");}

        // 1. marker number
        if(fread(&cur_m, sizeof(int32_t), 1, fp) != 1){
            LOGGER.e(0,"Read " + infile + " error (m)");
        }    
        if(cur_m != numSnpInRegion[i]){
            LOGGER.e(0,"In LD block " + block->ID + ", inconsistent marker number to marker information in " + infile);
        }
        // 2. ncol of eigenVec (number of eigenvalues)
        if(fread(&cur_k, sizeof(int32_t), 1, fp) != 1){
            LOGGER.e(0,"In LD block " + block->ID + ", error about number of eigenvalues in  " + infile);
        }
        // 3. sum of all positive eigenvalues
        if(fread(&sumPosEigVal, sizeof(float), 1, fp) != 1){
            LOGGER.e(0,"In LD block " + block->ID + ", error about the sum of positive eigenvalues in " + infile);
        }
        // 4. eigenCutoff
        if(fread(&oldEigenCutoff, sizeof(float), 1, fp) != 1){
            LOGGER.e(0,"In LD block " + block->ID + ", error about eigen cutoff used in " + infile);
        }
        // 5. eigenvalues
        VectorXf lambda(cur_k);
        if(fread(lambda.data(), sizeof(float), cur_k, fp) != cur_k){
            LOGGER.e(0,"In LD block " + block->ID + ",size error about eigenvalues in " + infile);
        }
        // 6. eigenvector
        MatrixXf U(cur_m, cur_k);
        uint64_t nElements = (uint64_t)cur_m * (uint64_t)cur_k;
        if(fread(U.data(), sizeof(float), nElements, fp) != nElements){
            LOGGER << "fread(U.data(), sizeof(double), nElements, fp): " << fread(U.data(), sizeof(float), nElements, fp) << endl;
            LOGGER << "nEle: " << nElements << " U.size: " << U.size() <<  " U.col: " << U.cols() << " row: " << U.rows() << endl;
            LOGGER.e(0,"In LD block " + block->ID + ",size error about eigenvectors in " + infile);
        }
        fclose(fp);
        bool haveValue = false;
        int revIdx = 0;
        if(oldEigenCutoff < eigenCutoff & i == 0){
            LOGGER << "Warning: current proportion of variance in LD block is set as " + to_string(eigenCutoff)+ ". But the proportion of variance is set as "<< to_string(oldEigenCutoff) + " in "  + infile + ".\n";
        }
        if (eigenCutoff < oldEigenCutoff) {
            truncateEigenMatrix(sumPosEigVal, eigenCutoff, lambda.cast<double>(), U.cast<double>(), eigenValLdBlock[i], eigenVecLdBlock[i]);
        } else {
            eigenValLdBlock[i] = lambda.cast<double>();
            eigenVecLdBlock[i] = U.cast<double>();
        }
        block->sumPosEigVal = sumPosEigVal;
        block->eigenvalues = lambda.cast<double>();
    }
    LOGGER << "GWAS LDM of" << numLDBlocks << " GWAS LD blocks to be include from [" << dirname << "]." << endl;
}

void Data::readEigenMatrixBinaryFileAndMakeWandQ(const string &dirname, const double eigenCutoff, const vector<VectorXd> &GWASeffects, const double nGWAS, const bool noscale, const bool makePseudoSummary){
    if (!Gadget::directoryExist(dirname)) {
        LOGGER.e(0, " cannot find the folder [" + dirname + "]");
    }
    
    vector<int>numSnpInRegion(numKeptLDBlocks);
    
    for(int i = 0; i < numKeptLDBlocks;i++){
        LDBlockInfo *block = keptLdBlockInfoVec[i];
        numSnpInRegion[i] = block->numSnpInBlock;
    }
    eigenValLdBlock.resize(numLDBlocks);
    eigenVecLdBlock.resize(numLDBlocks);
    wcorrBlocks.resize(numKeptLDBlocks);
    numSnpsBlock.resize(numKeptLDBlocks);
    numEigenvalBlock.resize(numKeptLDBlocks);
    Qblocks.resize(numKeptLDBlocks);

    
    //Constructing pseudo summary statistics for training and validation data sets, with 90% sample size for training and 10% for validation
    
    float n_trn, n_val;
    if (makePseudoSummary) {
        n_trn = 0.9*float(numKeptInds);
        n_val = numKeptInds - n_trn;
        pseudoGwasNtrn = n_trn;

        pseudoGwasEffectTrn.resize(numKeptLDBlocks);
        pseudoGwasEffectVal.resize(numKeptLDBlocks);
        b_val.setZero(numIncdSnps);
    }

// #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < numKeptLDBlocks; i++){
        LDBlockInfo *block = keptLdBlockInfoVec[i];
        int32_t cur_m = 0;
        int32_t cur_k = 0;
        float sumPosEigVal = 0;
        float oldEigenCutoff =0;
        
        string infile = dirname + "/block" + block->ID + ".eigen.bin";
        FILE *fp = fopen(infile.c_str(), "rb");
        if(!fp){LOGGER.e(0, " can not open the file [" + infile + "] to read.");}

        // 1. marker number
        if(fread(&cur_m, sizeof(int32_t), 1, fp) != 1){
            LOGGER.e(0, "Read " + infile + " error (m)");
        }
                
        if(cur_m != numSnpInRegion[i]){
            LOGGER.e(0, "In LD block " + block->ID + ", inconsistent marker number to marker information in " + infile);
        }
        // 2. ncol of eigenVec (number of eigenvalues)
        if(fread(&cur_k, sizeof(int32_t), 1, fp) != 1){
            LOGGER.e(0, "In LD block " + block->ID + ", error about number of eigenvalues in  " + infile);
        }
        // 3. sum of all positive eigenvalues
        if(fread(&sumPosEigVal, sizeof(float), 1, fp) != 1){
            LOGGER.e(0, "In LD block " + block->ID + ", error about the sum of positive eigenvalues in " + infile);
        }
        // 4. eigenCutoff
        if(fread(&oldEigenCutoff, sizeof(float), 1, fp) != 1){
            LOGGER.e(0, "In LD block " + block->ID + ", error about eigen cutoff used in " + infile);
        }
        // 5. eigenvalues
        VectorXf lambda(cur_k);
        if(fread(lambda.data(), sizeof(float), cur_k, fp) != cur_k){
            LOGGER.e(0, "In LD block " + block->ID + ",size error about eigenvalues in " + infile);
        }
        // 6. eigenvector
        MatrixXf U(cur_m, cur_k);
        uint64_t nElements = (uint64_t)cur_m * (uint64_t)cur_k;
        if(fread(U.data(), sizeof(float), nElements, fp) != nElements){
            cout << "fread(U.data(), sizeof(float), nElements, fp): " << fread(U.data(), sizeof(float), nElements, fp) << endl;
            cout << "nEle: " << nElements << " U.size: " << U.size() <<  " U.col: " << U.cols() << " row: " << U.rows() << endl;
            LOGGER.e(0, "In LD block " + block->ID + ",size error about eigenvectors in " + infile);
        }
        bool haveValue = false;
        int revIdx = 0;
        if(oldEigenCutoff < eigenCutoff & i == 0){
            cout << "Warning: current proportion of variance in LD block is set as " + to_string(eigenCutoff)+ ". But the proportion of variance is set as "<< to_string(oldEigenCutoff) + " in "  + infile + ".\n";
            // LOGGER.e(0, "");
        }
        if (eigenCutoff < oldEigenCutoff) {
            truncateEigenMatrix(sumPosEigVal, eigenCutoff, lambda.cast<double>(), U.cast<double>(), eigenValLdBlock[i], eigenVecLdBlock[i]);
        } else {
            eigenValLdBlock[i] = lambda.cast<double>();
            eigenVecLdBlock[i] = U.cast<double>();
        }
        block->sumPosEigVal = sumPosEigVal;
        block->eigenvalues = lambda.cast<double>();
        
        // make w and Q
        VectorXd sqrtLambda = eigenValLdBlock[i].array().sqrt();
        wcorrBlocks[i] = (1.0/sqrtLambda.array()).matrix().asDiagonal() * (eigenVecLdBlock[i].transpose() * GWASeffects[i] );
        //MatrixXf tmpQblocks = sqrtLambda.asDiagonal() * eigenVecLdBlock[i].transpose();
        //MatrixDat matrixDat = MatrixDat(block->snpNameVec, tmpQblocks);
        Qblocks[i] = sqrtLambda.asDiagonal() * eigenVecLdBlock[i].transpose();
        
        if (noscale) {
            VectorXd Dsqrt(block->numSnpInBlock);
            for (unsigned j=0; j<block->numSnpInBlock; ++j) {
                SnpInfo *snp = block->snpInfoVec[j];
                Dsqrt[j] = sqrt(snp->twopq);
            }
            Qblocks[i] = Qblocks[i] * Dsqrt.asDiagonal();
        }

        numSnpsBlock[i] = Qblocks[i].cols();
        numEigenvalBlock[i] = Qblocks[i].rows();
        
        // make pseudo summary data
        if (makePseudoSummary) {
            long size = eigenValLdBlock[i].size();
            VectorXd rnd(size);
            for (unsigned j=0; j<size; ++j) {
                rnd[j] = Stat::snorm();
            }
            
            pseudoGwasEffectTrn[i] = gwasEffectInBlock[i] + sqrt(1.0/n_trn - 1.0/nGWASblock[i]) * eigenVecLdBlock[i] * (eigenValLdBlock[i].array().sqrt().matrix().asDiagonal() * rnd);

            pseudoGwasEffectVal[i] = nGWASblock[i]/n_val * gwasEffectInBlock[i] - n_trn/n_val * pseudoGwasEffectTrn[i];
            b_val.segment(block->startSnpIdx, block->numSnpInBlock) = pseudoGwasEffectVal[i];
        }
        
        // eigenVecLdBlock[i].resize(0,0);
    }

    nGWASblock.resize(numKeptLDBlocks);
    for (unsigned i = 0; i < numKeptLDBlocks; i++){
        LDBlockInfo *ldblock = keptLdBlockInfoVec[i];
        nGWASblock[i] = nGWAS;
    }
}

void Data::UseLDBlockEigenMakeWAndQblocks(const vector<VectorXd> &GWASeffects, const double nGWAS, const bool noscale,const bool makePseudoSummary){
    ///// 
    float n_trn, n_val;
    if (makePseudoSummary) {
        n_trn = 0.9*float(numKeptInds);
        n_val = numKeptInds - n_trn;
        pseudoGwasNtrn = n_trn;
        pseudoGwasEffectTrn.resize(numKeptLDBlocks);
        pseudoGwasEffectVal.resize(numKeptLDBlocks);
        b_val.setZero(numIncdSnps);
    }
    wcorrBlocks.resize(numKeptLDBlocks);
    numSnpsBlock.resize(numKeptLDBlocks);
    numEigenvalBlock.resize(numKeptLDBlocks);
    Qblocks.resize(numKeptLDBlocks);
    VectorXd sqrtLambda;

    for(int i = 0,k = 0; i < numLDBlocks; i++){
        LDBlockInfo *block = keptLdBlockInfoVec[i];
        if(!block->kept) continue;
        sqrtLambda = eigenValLdBlock[i].array().sqrt();
        wcorrBlocks[k] = (1.0/sqrtLambda.array()).matrix().asDiagonal() * (eigenVecLdBlock[i].transpose() * GWASeffects[i] );
    
        Qblocks[k] = sqrtLambda.asDiagonal() * eigenVecLdBlock[i].transpose();
        if (noscale) {
            VectorXd Dsqrt(block->numSnpInBlock);
            for (unsigned j=0; j<block->numSnpInBlock; ++j) {
                SnpInfo *snp = block->snpInfoVec[j];
                Dsqrt[j] = sqrt(snp->twopq);
            }
            Qblocks[k] = Qblocks[i] * Dsqrt.asDiagonal();
        }

        numSnpsBlock[k] = Qblocks[i].cols();
        numEigenvalBlock[k] = Qblocks[i].rows();

        // make pseudo summary data
        if (makePseudoSummary) {
            long size = eigenValLdBlock[i].size();
            VectorXd rnd(size);
            for (unsigned j=0; j<size; ++j) {
                rnd[j] = Stat::snorm();
            }
            pseudoGwasEffectTrn[k] = gwasEffectInBlock[i] + sqrt(1.0/n_trn - 1.0/nGWASblock[i]) * eigenVecLdBlock[i] * (sqrtLambda.matrix().asDiagonal() * rnd);
            pseudoGwasEffectVal[k] = nGWASblock[i]/n_val * gwasEffectInBlock[i] - n_trn/n_val * pseudoGwasEffectTrn[k];
            b_val.segment(block->startSnpIdx, block->numSnpInBlock) = pseudoGwasEffectVal[i];
        }
        k++;
        if(k == numKeptLDBlocks) break;
    }

    nGWASblock.resize(numKeptLDBlocks);
    for (unsigned i = 0; i < numKeptLDBlocks; i++){
        LDBlockInfo *ldblock = keptLdBlockInfoVec[i];
        nGWASblock[i] = nGWAS;
    }
}

void Data::truncateEigenMatrix(const float sumPosEigVal, const double eigenCutoff, const VectorXd &oriEigenVal, const MatrixXd &oriEigenVec, VectorXd &newEigenVal, MatrixXd &newEigenVec){
    int revIdx = oriEigenVal.size();
    VectorXd cumsumNonNeg(revIdx);
    cumsumNonNeg.setZero();
    revIdx = revIdx -1;
    if(oriEigenVal(revIdx) < 0) LOGGER << "Error, all eigenvector are negative" << endl;
    cumsumNonNeg(revIdx) = oriEigenVal(revIdx);
    revIdx = revIdx -1;
    
    while(oriEigenVal(revIdx) > 1e-10 ){
        cumsumNonNeg(revIdx) = oriEigenVal(revIdx) + cumsumNonNeg(revIdx + 1);
        revIdx =revIdx - 1;
        if(revIdx <= 0) break;
    }
    // LOGGER << "revIdx: " << revIdx << endl;
    // LOGGER << "size: " << oriEigenVal.size()  << " eigenVal: " << oriEigenVal.head(10) << endl;
    // LOGGER << "cumsumNoNeg: " << cumsumNonNeg << endl;
    cumsumNonNeg = cumsumNonNeg/sumPosEigVal;  // calcualte the cumulative proportion of variance explained
    bool haveValue = false;
    // LOGGER << "cumsumNonNeg: " << cumsumNonNeg << endl;
    // LOGGER << "revIdx : " << revIdx << endl;
    // LOGGER << "revIdx: "  << revIdx  << endl;
    for (revIdx = revIdx + 1; revIdx < oriEigenVal.size(); revIdx ++ ){
        // LOGGER << "revIdx: " << revIdx << " cumsumNonNeg: " << cumsumNonNeg(revIdx) << endl;
        if(eigenCutoff >= cumsumNonNeg(revIdx) ){
            revIdx = revIdx -1;
            haveValue = true;
            break;
        }
    }
    // LOGGER << "revIdx : " << revIdx << endl;
    if(!haveValue) revIdx = oriEigenVal.size() - 1;
    // LOGGER << "cumsumNonNeg: " << cumsumNonNeg.size() << endl;
    // LOGGER << "cumsumNoNeg: " << cumsumNonNeg << endl;
    // LOGGER << "revIdx: "  << revIdx  << endl;

    newEigenVec = oriEigenVec.rightCols(oriEigenVal.size() - revIdx);
    newEigenVal = oriEigenVal.tail(oriEigenVal.size() - revIdx);
    // LOGGER << "eigenValAdjusted size: " << eigenValAdjusted.size() << " eigenValue eventually: " << eigenValAdjusted << endl;
    //eigenvalueNum = eigenVal.size() - revIdx;
    // LOGGER << endl;

}

void Data::readEigenMatrix(const string &dirname, const double eigenCutoff,const bool qcBool){
    LOGGER  << "............................" << endl;
    LOGGER << "Reading GWAS low-rank LD matrices..." << endl;
    LOGGER  << "............................" << endl;
    if(boost::filesystem::is_directory(dirname)) {
        LOGGER << "Reading GWAS low-rank LD matrices (SBayesRC format) ..." << endl;
        readSBayesRCBlockLdmInfoFile(dirname + "/ldm.info");
        readSBayesRCBlockLdmSnpInfoFile(dirname + "/snp.info");
        if(qcBool) readEigenMatrixBinaryFile(dirname,eigenCutoff);
    } else {
        string filename = dirname + ".eigen.ldblock";
        string LDMatType = "";
        readEigenMatLDBlockSnpInfoFile( filename + ".snp.info" ,LDMatType);
        readEigenMatLDBlockInfoFile(filename + ".info");
        readEigenMatBinFile(dirname,eigenCutoff,"ldblock");
    }
}

vector<LDBlockInfo*> Data::makeKeptLDBlockInfoVec(const vector<LDBlockInfo*> &ldBlockInfoVec){
    vector<LDBlockInfo*> keptLDBlock;
    ldblockNames.clear();
    LDBlockInfo * ldblock = NULL;
    for (unsigned i=0, j=0; i< numLDBlocks; ++i) {
        ldblock = ldBlockInfoVec[i];
        if(ldblock->kept) {
            ldblock->index = j++;  // reindex inds
            keptLDBlock.push_back(ldblock);
            ldblockNames.push_back(ldblock->ID);
        }
    }
    return keptLDBlock;
}

void Data::buildMMEigen(const string &dirname, const bool sampleOverlap, const double eigenCutoff, const bool noscale,const bool gwasInfoOnly){
    // check all available LD Blocks.
    includeMatchedBlocks();
    for (unsigned i=0; i<numIncdSnps; ++i) {
        SnpInfo *snp = incdSnpInfoVec[i];
        if (snp->gwas_b == -999) {
            LOGGER.e(0," SNP " + snp->rsID + " in the LD reference has no summary data. Run --impute-summary first.");
        }
    }
    LOGGER  << "............................" << endl;
    LOGGER << "Building GWAS various input parameters based on " << numKeptLDBlocks << " LD Blocks with " << numIncdSnps << " SNPs ... " << endl;
    LOGGER  << "............................" << endl;
    if(gwasInfoOnly) {
        LOGGER << "Construct various maps to indicate the relationship between the complex trait and GWAS." << endl;
        ConstructGwasEqtlGeneMaps();
    }
    LOGGER << "Scale GWAS effect size assuming unit GWAS trait phenotypic variance (Jian Yang et.al. (2012))..." << endl;
    scaleGwasEffects();

    LOGGER << "Construct corrected GWAS trait phenotype vector based on GWAS low-rank LD matrices (Lloyd-Jones, Zeng et.al. (2019))..." << endl;

    // read eigen LD matrix here 
    if(boost::filesystem::is_directory(dirname)) {
        readEigenMatrixBinaryFileAndMakeWandQ(dirname, eigenCutoff, gwasEffectInBlock, numKeptInds, noscale, true);
     } else {
        UseLDBlockEigenMakeWAndQblocks(gwasEffectInBlock, numKeptInds, noscale);
     }

    QblocksDat.clear();
    SnpInfo *snp;
    for (unsigned i=0; i<numKeptLDBlocks; ++i) {
        LDBlockInfo* ldblock = keptLdBlockInfoVec[i];
        MatrixDat matrixDat = MatrixDat(ldblock->snpNameVec, Qblocks[i]);
        QblocksDat.push_back(matrixDat);
        
    }

    if(gwasInfoOnly) {
        LOGGER  << "............................" << endl;
        LOGGER << "Summary of model input parameters:" << endl;
        LOGGER  << "............................" << endl;
        LOGGER << boost::format("%60s %8s %8s\n") %"" %"Mean" %"SD";
        LOGGER << boost::format("%60s %8.3f %8.3f\n") %"GWAS SNP Phenotypic variance" %Gadget::calcMean(varySnp) %sqrt(Gadget::calcVariance(varySnp));
        LOGGER << boost::format("%60s %8.3f %8.3f\n") %"GWAS SNP heterozygosity" %Gadget::calcMean(snp2pq) %sqrt(Gadget::calcVariance(snp2pq));
        LOGGER << boost::format("%60s %8.0f %8.0f\n") %"GWAS SNP sample size" %Gadget::calcMean(n) %sqrt(Gadget::calcVariance(n));
        LOGGER << boost::format("%60s %8.0f %8.0f\n") %"GWAS SNP allele flip number (ref: GWAS LD)" %numGWASFlip %0;
        LOGGER << boost::format("%60s %8.3f %8.3f\n") %"GWAS SNP effect (in genotype SD unit)" %Gadget::calcMean(b) %sqrt(Gadget::calcVariance(b));
        LOGGER << boost::format("%60s %8.3f %8.3f\n") %"GWAS SNP SE" %Gadget::calcMean(se) %sqrt(Gadget::calcVariance(se));
        LOGGER << boost::format("%60s %8.3f %8.3f\n") %"GWAS LD block size" %Gadget::calcMean(numSnpsBlock) %sqrt(Gadget::calcVariance(numSnpsBlock));
        LOGGER << boost::format("%60s %8.3f %8.3f\n") %"GWAS block rank" %Gadget::calcMean(numEigenvalBlock) %sqrt(Gadget::calcVariance(numEigenvalBlock));
        LOGGER << "............................" << endl;
        LOGGER << "Begin to model inference process..." << endl;
        LOGGER  << "............................" << endl;
    }
    lowRankModel = true;
}


void Data::includeMatchedBlocks(const bool haveEqtlInfo){
    // this step is to construct gwasSnp2geneVec
//    LOGGER << "Matching blocks..." << endl;
    SnpInfo * snp;
    LDBlockInfo * ldblock;
    
    for (unsigned i=0; i<numLDBlocks; ++i){
        ldblock = ldBlockInfoVec[i];
        ldblock->block2GwasSnpVec.clear();
        ldblock->snpInfoVec.clear();
    }
    for (unsigned j=0; j<numIncdSnps; ++j){
        snp = incdSnpInfoVec[j];
        ldblock = ldBlockInfoMap[snp->block];
        ldblock->block2GwasSnpVec.push_back(j);
        ldblock->snpInfoVec.push_back(snp);
        if(haveEqtlInfo){
            if(snp->iseQTL){
                ldblock->block2EqtlBoolVec.push_back(true);
            } else {
                ldblock->block2EqtlBoolVec.push_back(false);
            }
        }
    }
    for (unsigned i=0; i<numLDBlocks; ++i){
        ldblock = ldBlockInfoVec[i];
        if(ldblock->block2GwasSnpVec.size() == 0){
            ldblock->kept = false;
        } else {
            ldblock->startSnpIdx = ldblock->block2GwasSnpVec[0];
            ldblock->endSnpIdx = ldblock->block2GwasSnpVec[ldblock->numSnpInBlock-1];
            ldblock->kept = true;
        }
    }
        
    keptLdBlockInfoVec = makeKeptLDBlockInfoVec(ldBlockInfoVec);
    numKeptLDBlocks = (unsigned) keptLdBlockInfoVec.size();

    ldblock2gwasSnpMap.clear();
//    LOGGER << "Construct map from ld to snp" << endl;
    for(unsigned i = 0; i < numKeptLDBlocks; i++){
        ldblock = keptLdBlockInfoVec[i];
        ldblock2gwasSnpMap.insert(pair<int, vector<int> > (i,ldblock->block2GwasSnpVec));
    }
    
    // LOGGER << numKeptLDBlocks << " GWAS LD blocks are included." << endl;
}


//void Data::constructWandQ(const double eigenCutoff, const bool noscale){
//    VectorXd nMinusOne;
//    snp2pq.resize(numIncdSnps);
//    D.resize(numIncdSnps);
//   // ZPZdiag.resize(numIncdSnps);
//    ZPy.resize(numIncdSnps);
//    b.resize(numIncdSnps);
//    n.resize(numIncdSnps);
//    nMinusOne.resize(numIncdSnps);
//    se.resize(numIncdSnps);
//    tss.resize(numIncdSnps);
//    SnpInfo *snp;
//    for (unsigned i=0; i<numIncdSnps; ++i) {
//        snp = incdSnpInfoVec[i];
//        snp->af = snp->gwas_af;
//        snp2pq[i] = snp->twopq = 2.0f*snp->gwas_af*(1.0-snp->gwas_af);
//        if(snp2pq[i]==0) LOGGER << "Error: SNP " << snp->rsID << " af " << snp->af << " has 2pq = 0." << endl;
//        D[i] = snp2pq[i]*snp->gwas_n;
//        b[i] = snp->gwas_b * sqrt(snp2pq[i]); // scale the marginal effect so that it's in per genotype SD unit
//        n[i] = snp->gwas_n;
//        nMinusOne[i] = snp->gwas_n - 1;
//        se[i]= snp->gwas_se * sqrt(snp2pq[i]);
//        tss[i] = D[i]*(n[i]*se[i]*se[i] + b[i]*b[i]);
//        ZPy[i] = n[i]*b[i];
//        //  D[i] = 1.0/(se[i]*se[i]+b[i]*b[i]/snp->gwas_n);  // NEW!
//        //  snp2pq[i] = snp->twopq = D[i]/snp->gwas_n;       // NEW!
//    }
//    LOGGER << endl;
//    LDBlockInfo * ldblock;
//    wcorrBlocks.resize(numKeptLDBlocks);
//    // Qblocks.resize(numKeptLDBlocks);
//    numSnpsBlock.resize(numKeptLDBlocks);
//    numEigenvalBlock.resize(numKeptLDBlocks);
//    Qblocks.clear();
//    VectorXd sqrtLambda;
//    // save gwas marginal effect into block
//    gwasEffectInBlock.resize(numKeptLDBlocks);
//    for (unsigned i = 0; i < numKeptLDBlocks; i++){
//        ldblock = keptLdBlockInfoVec[i];
//        gwasEffectInBlock[i] = b(ldblock->block2GwasSnpVec);
//        // calculate wbcorr and Qblocks
//        sqrtLambda = eigenValLdBlock[i].array().sqrt();
//        //LOGGER << eigenVecLdBlock[i].transpose().rows() << " " << eigenVecLdBlock[i].transpose().cols() << " " << gwasEffectInBlock[i].size() << endl;
//        wcorrBlocks[i] = (1.0/sqrtLambda.array()).matrix().asDiagonal() * (eigenVecLdBlock[i].transpose() * gwasEffectInBlock[i] );
//        // LOGGER << "eigenVecLdBlock[i]: " << eigenVecLdBlock[i] << endl;
//        // LOGGER << "wcorrBlocks[i]: " << wcorrBlocks[i] << endl;
//        // LOGGER << "sqrtLambda: " << sqrtLambda << endl;
//        // LOGGER << gwasEffectInBlock[i] << endl;
//        MatrixXd tmpQblocks = sqrtLambda.asDiagonal() * eigenVecLdBlock[i].transpose();
//        MatrixDat matrixDat = MatrixDat(ldblock->snpNameVec,tmpQblocks );
//        // LOGGER << "Qblock: " << endl;
//        // LOGGER << matrixDat.values << endl;
//        Qblocks.push_back(matrixDat);
//        numSnpsBlock[i] = Qblocks[i].ncol;
//        numEigenvalBlock[i] = Qblocks[i].nrow;
//    }
//    
//    //b.array() -= b.mean();  // DO NOT CENTER b
//    // estimate phenotypic variance based on the input allele frequencies in GWAS
//
//    //  Vp_buf = h_buf * N_buf * se_buf * se_buf + h_buf * b_buf * b_buf * N_buf / (N_buf - 1.0);
//    //VectorXd ypySrt = D.array()*n.array()*se.array().square() + D.array() * b.array().square() * n.array() / nMinusOne.array();
//    VectorXd ypySrt = D.array()*(n.array()*se.array().square()+b.array().square());
//    VectorXd varpSrt = ypySrt.array()/n.array();
//    std::sort(ypySrt.data(), ypySrt.data() + ypySrt.size());
//    std::sort(varpSrt.data(), varpSrt.data() + varpSrt.size());
//    ypy = ypySrt[ypySrt.size()/2];  // median
//    varPhenotypic = varpSrt[varpSrt.size()/2];
//    //LOGGER << "varPhenotypic: " << varPhenotypic << endl;
//    VectorXd nSrt = n;
//    std::sort(nSrt.data(), nSrt.data() + nSrt.size());
//    numKeptInds = nSrt[nSrt.size()/2]; // median
//    
//    nGWASblock.resize(numKeptLDBlocks);
//    for (unsigned i=0; i<numKeptLDBlocks; ++i) {
//        nGWASblock[i] = numKeptInds;
//    }
//
//    for (unsigned i=0; i<numIncdSnps; ++i) {
//        snp = incdSnpInfoVec[i];
//        D[i] = varPhenotypic/(se[i]*se[i]+b[i]*b[i]/snp->gwas_n);  // NEW!
//        snp2pq[i] = snp->twopq = D[i]/snp->gwas_n;       // NEW!
//        tss[i] = D[i]*(n[i]*se[i]*se[i] + b[i]*b[i]);
//        // Need to adjust R and C models X'X matrix depending scale of genotypes or not
//        if (noscale == true) {
//            D[i] = snp2pq[i]*snp->gwas_n;
//        } else {
//            D[i] = snp->gwas_n;
//        }
//    }
//    // LOGGER << endl << "snp2pq: " << endl << snp2pq << endl;
//    //ypy = numKeptInds;
//    // NEW END
//    // data summary
//    LOGGER << "\nData summary:" << endl;
//    LOGGER << boost::format("%40s %8s %8s\n") %"" %"mean" %"sd";
//    LOGGER << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP Phenotypic variance" %Gadget::calcMean(varpSrt) %sqrt(Gadget::calcVariance(varpSrt));
//    LOGGER << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP heterozygosity" %Gadget::calcMean(snp2pq) %sqrt(Gadget::calcVariance(snp2pq));
//    LOGGER << boost::format("%40s %8.0f %8.0f\n") %"GWAS SNP sample size" %Gadget::calcMean(n) %sqrt(Gadget::calcVariance(n));
//    LOGGER << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP effect (in genotype SD unit)" %Gadget::calcMean(b) %sqrt(Gadget::calcVariance(b));
//    LOGGER << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP SE" %Gadget::calcMean(se) %sqrt(Gadget::calcVariance(se));
//    LOGGER << boost::format("%40s %8.3f %8.3f\n") %"LD block size" %Gadget::calcMean(numSnpsBlock) %sqrt(Gadget::calcVariance(numSnpsBlock));
//    LOGGER << boost::format("%40s %8.3f %8.3f\n") %"LD block rank" %Gadget::calcMean(numEigenvalBlock) %sqrt(Gadget::calcVariance(numEigenvalBlock));
//}


void Data::mergeLdmInfo(const string &outLDmatType, const string &dirname) {
    
    if (outLDmatType != "block") LOGGER.e(0, " --merge-ldm-info only works for block LD matrices at the moment!");
    
    string dir_path = dirname; // Replace with your folder path
    string search_str = outLDmatType;
    DIR* dirp = opendir(dir_path.c_str());
    
    if (dirp == NULL) {
        LOGGER.e(0, " opening directory [" + dirname + "]");
    }
    
    // find out all ldm in the folder
    vector<string> file_list;
    
    dirent* dp;
    while ((dp = readdir(dirp)) != NULL) {
        string file_name = dp->d_name;
        if (file_name.find(search_str) != string::npos) {
            file_list.push_back(file_name);
        }
    }
    
    closedir(dirp);
    
    set<unsigned> blockIdxSet;
    blockIdxSet.clear();
    
    for (vector<string>::iterator it = file_list.begin(); it != file_list.end(); ++it) {
        size_t block_pos = it->find(search_str);
        size_t dot_pos = it->find_first_of(".", block_pos);
        if (block_pos != string::npos && dot_pos != string::npos) {
            string block_num_str = it->substr(block_pos + search_str.size(), dot_pos - block_pos - search_str.size());
            int block_num = atoi(block_num_str.c_str());
            blockIdxSet.insert(block_num);
        }
    }
    
    if (blockIdxSet.size() == 0) {
        LOGGER.e(0, " there is no info file to merge in folder [" + dirname + "].");
    }
    
    unsigned nldm = blockIdxSet.size();
    
    set<unsigned>::iterator it = blockIdxSet.begin();
    
    string outSnpInfoFile = dirname + "/snp.info";
    string outldmInfoFile = dirname + "/ldm.info";
    
    ofstream out1(outSnpInfoFile.c_str());
    out1 << boost::format("%6s %15s %10s %10s %15s %6s %6s %12s %10s %10s\n")
    % "Chrom"
    % "ID"
    % "Index"
    % "GenPos"
    % "PhysPos"
    % "A1"
    % "A2"
    % "A1Freq"
    % "N"
    % "Block";
    
    ofstream out2(outldmInfoFile.c_str());
    out2 << boost::format("%10s %6s %15s %15s %15s %15s %12s\n")
    % "Block"
    % "Chrom"
    % "StartSnpIdx"
    % "StartSnpID"
    % "EndSnpIdx"
    % "EndSnpID"
    % "NumSnps";
    
    
    
    unsigned snpIdx = 0;
    unsigned ldmIdx = 0;
    
    for (unsigned i=0; i<nldm; ++i) {
        
        string snpInfoFile = dirname + "/block" + to_string(*it) + ".snp.info";
        string ldmInfoFile = dirname + "/block" + to_string(*it) + ".ldm.info";
        
        // read snp info file
        ifstream in1(snpInfoFile.c_str());
        if (!in1) LOGGER.e(0, " can not open the file [" + snpInfoFile + "] to read.");
        cout << "Reading SNP info from file [" + snpInfoFile + "]." << endl;
        string header;
        string id, allele1, allele2;
        int chr, physPos,ld_n;
        float genPos;
        float allele1Freq;
        int idx = 0;
        int index;
        string blockID;
        map<string, int> snpID2index;
        getline(in1, header);
        while (in1 >> chr >> id >> index >> genPos >> physPos >> allele1 >> allele2 >> allele1Freq >> ld_n >> blockID) {
            out1 << boost::format("%6s %15s %10s %10s %15s %6s %6s %12f %10s %10s\n")
            % chr
            % id
            % snpIdx
            % genPos
            % physPos
            % allele1
            % allele2
            % allele1Freq
            % ld_n
            % blockID;
            snpID2index[id] = snpIdx;
            ++snpIdx;
        }
        in1.close();
        
        
        // read ldm info file
        ifstream in2(ldmInfoFile.c_str());
        if (!in2) LOGGER.e(0, " can not open the file [" + ldmInfoFile + "] to read.");
        cout << "Reading LDM info from file [" + ldmInfoFile + "]." << endl;
        
        int blockStart, blockEnd, snpNum;
        string startSnpID, endSnpID;
        getline(in2, header);
        while (in2 >> id >> chr >> blockStart >> startSnpID >> blockEnd >> endSnpID >> snpNum) {
            out2 << boost::format("%10s %6s %15s %15s %15s %15s %12s\n")
            % id
            % chr
            % snpID2index[startSnpID]
            % startSnpID
            % snpID2index[endSnpID]
            % endSnpID
            % snpNum;
            
            ++ldmIdx;
        }
        in2.close();
        
        ++it;
        
        remove(snpInfoFile.c_str());
        remove(ldmInfoFile.c_str());
    }
    
    out1.close();
    out2.close();
    
    cout << "Written " << snpIdx << " SNPs info into file [" + outSnpInfoFile + "]." << endl;
    cout << "Written " << ldmIdx << " LDMs info into file [" + outldmInfoFile + "]." << endl;
    
}


void Data::mergeBlockGwasSummary(const string &gwasSummaryFile, const string &title) {
    string dir_path = "."; // Replace with your folder path
    string search_str1 = gwasSummaryFile + ".block";
    string search_str2 = ".ma";
    DIR* dirp = opendir(dir_path.c_str());
    
    if (dirp == NULL) {
        LOGGER.e(0,"Error opening directory [" + dir_path + "]");
    }
    
    // find out all .ma files in the folder
    vector<string> file_list;
    
    dirent* dp;
    while ((dp = readdir(dirp)) != NULL) {
        string file_name = dp->d_name;
        if (file_name.find(search_str1) != string::npos && file_name.find(search_str2) != string::npos) {
            file_list.push_back(file_name);
        }
    }
    
    closedir(dirp);
    
    map<unsigned, string> blockIdxMap;
    blockIdxMap.clear();
    
    for (vector<string>::iterator it = file_list.begin(); it != file_list.end(); ++it) {
        size_t block_pos = it->find(search_str1);
        size_t dot_pos = it->find_first_of(".", block_pos);
        if (block_pos != string::npos && dot_pos != string::npos) {
            string block_num_str = it->substr(block_pos + search_str1.size(), dot_pos - block_pos - search_str1.size());
            int block_num = atoi(block_num_str.c_str());
            blockIdxMap[block_num] = *it;
        }
    }
    
    if (blockIdxMap.size() == 0) {
        LOGGER.e(0,"there is no info file to merge in folder [" + dir_path + "].");
    }
    
    unsigned nBlk = blockIdxMap.size();
    
    LOGGER << "Merging GWAS summary statistics files across " + to_string(nBlk) + " blocks..." << endl;
    
    map<unsigned, string>::iterator it = blockIdxMap.begin();
    
    string outMaFile = dir_path + "/" + title + ".ma";
    
    ofstream out(outMaFile.c_str());
    out << boost::format("%15s %10s %10s %15s %15s %15s %15s %15s\n") % "SNP" % "A1" % "A2" % "freq" % "b" % "se" % "p" % "N";
    
    unsigned snpIdx = 0;
    
    for (unsigned i=0; i<nBlk; ++i) {
        
        string mafile = it->second;
        
        // read snp info file
        ifstream in(mafile.c_str());
        if (!in) LOGGER.e(0,"can not open the file [" + mafile + "] to read.");
        LOGGER << "Reading summary statistics from file [" + mafile + "]." << endl;

        string header;
        getline(in, header);

        string id, allele1, allele2, freq, b, se, pval, n;
        while (in >> id >> allele1 >> allele2 >> freq >> b >> se >> pval >> n) {
            out << boost::format("%15s %10s %10s %15s %15s %15s %15s %15s\n")
            % id
            % allele1
            % allele2
            % freq
            % b
            % se
            % pval
            % n;
            ++snpIdx;
        }
        in.close();
        
        ++it;
    }
    
    out.close();
    
    LOGGER << "Written " << snpIdx << " SNPs info into file [" + outMaFile + "]." << endl;
}

void Data::constructPseudoSummaryData(){
    LOGGER << "Constructing pseudo summary statistics for training and validation data sets, with 90% sample size for training and 10% for validation." << endl;
    
    pseudoGwasEffectTrn.resize(numKeptLDBlocks);
    pseudoGwasEffectVal.resize(numKeptLDBlocks);
    
    double n_trn = 0.9*double(numKeptInds);
    double n_val = numKeptInds - n_trn;
    pseudoGwasNtrn = n_trn;
    b_val.resize(numIncdSnps);

    for (unsigned i=0; i<numKeptLDBlocks; ++i) {
        LDBlockInfo* block = keptLdBlockInfoVec[i];
        
        long size = eigenValLdBlock[i].size();
        VectorXd rnd(size);
        for (unsigned j=0; j<size; ++j) {
            rnd[j] = Stat::snorm();
        }
        
        pseudoGwasEffectTrn[i] = gwasEffectInBlock[i] + sqrt(1.0/n_trn - 1.0/nGWASblock[i]) * eigenVecLdBlock[i] * (eigenValLdBlock[i].array().sqrt().matrix().asDiagonal() * rnd);

        pseudoGwasEffectVal[i] = nGWASblock[i]/n_val * gwasEffectInBlock[i] - n_trn/n_val * pseudoGwasEffectTrn[i];
        b_val.segment(block->startSnpIdx, block->numSnpInBlock) = pseudoGwasEffectVal[i];
    }
    
    
}

void Data::constructWandQ(const vector<VectorXd> &GWASeffects, const double nGWAS, const bool noscale) {
    wcorrBlocks.resize(numKeptLDBlocks);
    numSnpsBlock.resize(numKeptLDBlocks);
    numEigenvalBlock.resize(numKeptLDBlocks);
    Qblocks.clear();

    for (unsigned i = 0; i < numKeptLDBlocks; i++){
        LDBlockInfo *ldblock = keptLdBlockInfoVec[i];
        // calculate wbcorr and Qblocks
        VectorXd sqrtLambda = eigenValLdBlock[i].array().sqrt();
        //LOGGER << eigenVecLdBlock[i].transpose().rows() << " " << eigenVecLdBlock[i].transpose().cols() << " " << gwasEffectInBlock[i].size() << endl;
        wcorrBlocks[i] = (1.0/sqrtLambda.array()).matrix().asDiagonal() * (eigenVecLdBlock[i].transpose() * GWASeffects[i] );
        // LOGGER << "eigenVecLdBlock[i]: " << eigenVecLdBlock[i] << endl;
        // LOGGER << "wcorrBlocks[i]: " << wcorrBlocks[i] << endl;
        // LOGGER << "sqrtLambda: " << sqrtLambda << endl;
        // LOGGER << gwasEffectInBlock[i] << endl;
//        MatrixXd tmpQblocks = sqrtLambda.asDiagonal() * eigenVecLdBlock[i].trßanspose();
//        MatrixDat matrixDat = MatrixDat(ldblock->snpNameVec, tmpQblocks);
//        // LOGGER << "Qblock: " << endl;
//        // LOGGER << matrixDat.values << endl;
//        Qblocks.push_back(matrixDat);
        Qblocks[i] = sqrtLambda.asDiagonal() * eigenVecLdBlock[i].transpose();
        
        if (noscale) {
            VectorXd Dsqrt(ldblock->numSnpInBlock);
            for (unsigned j=0; j<ldblock->numSnpInBlock; ++j) {
                SnpInfo *snp = ldblock->snpInfoVec[j];
                Dsqrt[j] = sqrt(snp->twopq);
            }
            Qblocks[i] = Qblocks[i] * Dsqrt.asDiagonal();
        }
        
        numSnpsBlock[i] = Qblocks[i].cols();
        numEigenvalBlock[i] = Qblocks[i].rows();
        
        eigenVecLdBlock[i].resize(0,0);
    }
    
    nGWASblock.resize(numKeptLDBlocks);
    for (unsigned i = 0; i < numKeptLDBlocks; i++){
        LDBlockInfo *ldblock = keptLdBlockInfoVec[i];
        nGWASblock[i] = nGWAS;
    }
}

void Data::scaleGwasEffects(){
    snp2pq.resize(numIncdSnps);
    scalingGWASFactorVec.resize(numIncdSnps);
    b.resize(numIncdSnps);
    n.resize(numIncdSnps);
    D.resize(numIncdSnps);
    se.resize(numIncdSnps);
    tss.resize(numIncdSnps); // only used in SBayesC
    ZPy.resize(numIncdSnps);
    varySnp.resize(numIncdSnps);

    SnpInfo *snp;
    int lineIdx = 0;
    auto startTime = std::chrono::steady_clock::now();
    for (unsigned i=0; i<numIncdSnps; ++i) {
        Gadget::showProgressBar(i, numIncdSnps, startTime,"Scale GWAS effects");
        double varpsPG=0,varpsDivid2pqPG = 0;
        snp = incdSnpInfoVec[i];
        snp->af = snp->gwas_af;
        snp2pq[i] = snp->twopq = 2.0f*snp->gwas_af*(1.0-snp->gwas_af);
        if(snp2pq[i]==0) LOGGER << "Error: SNP " << snp->rsID << " af " << snp->af << " has 2pq = 0." << endl;
        n[i] = snp->gwas_n;
        D[i] = snp2pq[i]*snp->gwas_n;
        b[i] = snp->gwas_b; // * sqrt(snp2pq[i]); // scale the marginal effect so that it's in per genotype SD unit
        se[i]= snp->gwas_se ; //* sqrt(snp2pq[i]);
        tss[i] = n[i]*(n[i]*se[i]*se[i] + b[i]*b[i]);
        ZPy[i] = n[i]*b[i];

        // start to scale
        if(lineIdx == 0) LOGGER.w(0,"Assuming non-scaled GWAS effect and se is used. Please double-check it.");
        // calculate scaling factor
        scalingGWASFactorVec[i] = sqrt(1/(n[i]* se[i] * se[i] + b[i]*b[i]));
        snp->scaleFactor = scalingGWASFactorVec[i];
        // check beta and se scaled or not
        // When generating beta effect under standardized genotype based on beta effect from GWAS using genotypes at 0/1/2 scale
        // one way is to time sqrt(snp2pq), in this way we also need to set phenotypic variance as the median of each SNP-based
        // phenotypic variance.
        // Another way is assuming phenotypic variance as 1, and then timeing a per-SNP constant scalingGWASFactorVec[i] 
        // if(false){
            // varpsPG =  D[i]*(n[i]*se[i]*se[i] + b[i]*b[i])/n[i];
            // varpsDivid2pqPG =  D[i]*(n[i]*se[i]*se[i] + b[i]*b[i])/(n[i] * snp2pq[i]);
            // b[i] = snp->gwas_b  * sqrt(snp2pq[i]); // scale the marginal effect so that it's in per genotype SD unit
            // se[i]= snp->gwas_se * sqrt(snp2pq[i]);
            // varySnp[i] = varpsPG;
        // }
        b[i] = snp->gwas_b  * (snp->scaleFactor); // scale the marginal effect so that it's in per genotype SD unit
        se[i]= snp->gwas_se * (snp->scaleFactor);
        varySnp[i] = 1.0;
        /////////////////////////////////////
        lineIdx ++;
    }
    // Estimate sample size
    VectorXd nSrt = n;
    std::sort(nSrt.data(), nSrt.data() + nSrt.size());
    numKeptInds = Gadget::findMedian(nSrt);
    varPhenotypic = Gadget::findMedian(varySnp);
    // LOGGER << "varPhenotypic: " << varPhenotypic << endl;

    // divide gwas marginal effect into blocks
    if (numKeptLDBlocks) {
        gwasEffectInBlock.resize(numKeptLDBlocks);
        nGWASblock.resize(numKeptLDBlocks);
        for (unsigned i = 0; i < numKeptLDBlocks; i++){
            LDBlockInfo *ldblock = keptLdBlockInfoVec[i];
            gwasEffectInBlock[i] = b(ldblock->block2GwasSnpVec);
            nGWASblock[i] = numKeptInds;
        }
    }
}
