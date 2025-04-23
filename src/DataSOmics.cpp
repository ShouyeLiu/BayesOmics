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
///////////////////////////////////////////////////
/// Read SMR BESD format
double Data::adjeQTLSE(double beta,double p)
{
    if(p==1){
        return 1e10;
    }else{
        // double z2=qchisq(p, 1);
        double z2=Stat::ChiSq::qchisq(p, 1);
        return fabs(beta/sqrt(z2));
    }
}

void Data::readGeneSampleSizeFile(const string &geneSamSizeFile){
    ifstream in(geneSamSizeFile.c_str());
    if (!in) LOGGER.e(0, " can not open the file [" + geneSamSizeFile + "] to read.");
    LOGGER << "Reading gene sample size from [" + geneSamSizeFile + "]." << endl;
    unsigned sampleSize;
    string snpID, geneID, geneOri;
    map<string, GeneInfo *>::iterator iterGene;
    map<string, EqtlInfo *>::iterator iterEqtl;
    unsigned idx = 0;
    Gadget::Tokenizer colData;
    string inputStr;
    string sep(", \t");
    string id;
    // check header 
    getline(in,inputStr);
    colData.getTokens(inputStr, sep);
    if(colData.size() == 2){
        LOGGER << "Gene sample size file has two columns and assume all SNPs within the same gene have the same sample size." << endl;
    } else if(colData.size() ==3){
        LOGGER << "Gene sample size file has three columns and assume all SNPs within the same gene may have a least one different sample size." << endl;
    } else {
        LOGGER.e(0,"There must have two or three columns in the [" + geneSamSizeFile + "].");
    }
    // now 
    while (getline(in,inputStr)) {
        colData.getTokens(inputStr, sep);
        if(colData.size() == 2){
            // have two columns: GeneID\tN
            geneID = colData[0];
            sampleSize = std::stod(colData[1]);
            iterGene = geneInfoMap.find(geneID);
            if(iterGene != geneInfoMap.end()){
                iterGene->second->sampleSize = sampleSize;
            }
        } else if(colData.size() ==3){
            // have three columns: GeneID\tSNPID\tN
            geneID = colData[0];
            snpID = colData[1];
            sampleSize = std::stod(colData[2]);
            iterGene = geneInfoMap.find(geneID);
            iterEqtl = eqtlInfoMap.find(snpID);
            if(iterGene != geneInfoMap.end() && iterEqtl != eqtlInfoMap.end()){
                iterGene->second->cisSnpSampleSizeMap.insert(pair<string, int>(iterEqtl->second->rsID,sampleSize));
            } 
        } else {
            // something wrong
            LOGGER.e(0,"the column number is not 2 or 3 in the line " + to_string(idx) + ".");
        }
        idx++;
    }
    in.close();
    LOGGER << idx << " genes with updated sample size to be included from [" + geneSamSizeFile + "]." << endl;
}

void Data::includeGene(const string &includeGeneFile){
    ifstream in(includeGeneFile.c_str());
    if (!in) LOGGER.e(0,"Can not open the file [" + includeGeneFile + "] to read.");
    for (unsigned i=0; i<numGenes ; ++i) {
        geneInfoVec[i]->kept = false;
    }
    map<string, GeneInfo *>::iterator it, end = geneInfoMap.end();
    string id;
    int numKeptGenesLocal = 0;
    while (in >> id) {
        it = geneInfoMap.find(id);
        if (it != end) {
            it->second->kept = true;
            numKeptGenesLocal ++;
        }
    }
    in.close();
    LOGGER << numKeptGenesLocal << " genes with updated sample size to be included from [" + includeGeneFile + "]." << endl;
}

void Data::excludeGene(const string &excludeGeneFile){
    ifstream in(excludeGeneFile.c_str());
    if (!in) LOGGER.e(0,"Can not open the file [" + excludeGeneFile + "] to read.");
    map<string, GeneInfo *>::iterator it, end = geneInfoMap.end();
    string id;
    int idx = 0;
    while (in >> id) {
        it = geneInfoMap.find(id);
        if (it != end) {
            it->second->kept = false;
            idx ++;
        }
    }
    in.close();
    LOGGER << to_string(idx) << " gene(s) were removed." << endl;
}

void Data::includeSpecificGene(const string &includeSpecificGeneID){
    for (unsigned i=0; i<numGenes ; ++i) {
        geneInfoVec[i]->kept = false;
    }
    map<string, GeneInfo *>::iterator it, end = geneInfoMap.end();
    string id =includeSpecificGeneID;
    it = geneInfoMap.find(id);
    if (it != end) {
            it->second->kept = true;
    } else {
        LOGGER.e(0,"Cannot find geneID: " + includeSpecificGeneID + " in the dataset, please check gene names.");
    }
}

////////////////////////////////////////////////////////////////////
/////////  Generate eigen LD matrix for gene region  ///////////////
////////////////////////////////////////////////////////////////////
///////////// Step 2.3 read LD matrix of gene regions info (3 functions)
int Data::readEigenMatGeneInfoFile(vector<GeneInfo *> &geneInfoVecLD , map<string, GeneInfo *> &geneInfoMapLD, const string &infoFile){
    // Read bim file: recombination rate is defined between SNP i and SNP i-1
    ifstream in(infoFile.c_str());
    if (!in) LOGGER.e(0, " can not open the file [" + infoFile + "] to read.");
    LOGGER << "Reading xQTL LDM gene info from file [" + infoFile + "]." << endl;
    // vector<GeneInfo*> geneInfoVecTmp;
    // map<string, GeneInfo*> geneInfoMapTmp;
    map<string,int> gene2snpMap;
    string header;
    string id;
    int  chr, start, end, windows, snpNum;
    int idx = 0;
    int snpCount =  1;
    string snp;
    GeneInfo *gene;
    map<string, GeneInfo*>::iterator iterGene;
    map<string, EqtlInfo*>::iterator iterEqtl;
    // Step 1. read gene info from ld gene matrix. In this step, if gene in LD file does not belong to 
    // gene set in BESD file. gene->kept will be set as false;
    getline(in, header);
    while (in >>chr>>id >>start>>end>> windows >> snp >> snpNum ) {
        if (gene2snpMap.insert(pair<string, int>(id + "_" + snp , snpCount)).second == false) {
            LOGGER.e(0, " Duplicate LDBlock-SNP pair found: \"" + id + "_" + snp + "\".");
        } else{
            // Step 1.1 build a new geneInfo firstly
            if(snpCount == 1){
                gene = new GeneInfo(idx++, id, chr);
                if(geneInfoMap.find(gene->ensemblID) == geneInfoMap.end()){
                    // set BESD genelist as benchmark, this LD gene is not included in BESD genelist;
                    gene->kept = false;
                }
                gene->start = start;
                gene->end   = end;
                gene->typSnpIdxBesd.clear(); // used to extract eqtl beta and se values;
                gene->typSnpIdxGeneLD.clear(); // used to extract eqtl with beta idx in gene ld 
                gene->impSnpIdxGeneLD.clear(); // used to extract eqtl idx need to impute
                gene->cisSnpNameVec.clear();
            }
            // Step 1.2 find wheter given gene eigen snp exists in GWAS ld reference or not
            if(snpInfoMap.find(snp) == snpInfoMap.end()){ 
                // if not, we need to remove this gene to make sure the consistence between GWAS LD and gene LD.
                // LOGGER.w(0,"gene " + gene->ensemblID + " was removed since eQTL " + snp + " from gene can not be found from GWAS LD reference.");
                gene->kept = false;
            } 
            // Step 1.3 continue to add eQTLs into gene
            if(snpCount < snpNum){
                gene->cisSnpNameVec.push_back(snp);
            }
            // Step 1.4 compare besd eqtl and gene eigen eqtl;
            if(gene->kept){
                iterGene = geneInfoMap.find(id);
                iterEqtl = eqtlInfoMap.find(snp);
                if(iterEqtl == eqtlInfoMap.end()){
                // besd eqtl cannot be found in gene eigen LD. we need to impute this.
                    gene->impSnpIdxGeneLD.push_back(snpCount - 1);
                } else {
                    // besd eqtl can be found in gene eigen LD, but check if besd eqtl is included or not
                    if(!iterEqtl->second->included){
                        gene->impSnpIdxGeneLD.push_back(snpCount - 1);
                    } else {
                        gene->typSnpIdxBesd.push_back(iterGene->second->cisSnpID2IdxMapInGene[snp]);
                        gene->typSnpIdxGeneLD.push_back(snpCount - 1);
                    }
                }
            }
            // Step 1.5 add the last eQTL into gene 
            if(snpCount == snpNum){
                gene->cisSnpNameVec.push_back(snp);
                gene->numSnpInGene = snpNum;
                // if the number of missing SNPs is to larger, we need to remove this gene
                if(gene->kept && double (gene->impSnpIdxGeneLD.size())/ (double) snpNum > SNPMISSRATE ){gene->kept = false;}
                if(gene->impSnpIdxGeneLD.size() != 0 ){gene->impCisSnpEffBool = true;}
                geneInfoVecLD.push_back(gene);
                if (geneInfoMapLD.insert(pair<string, GeneInfo*>(id, gene)).second == false) {
                    LOGGER.e(0, " Duplicate LD block ID found: \"" + id + "\".");}
                snpCount = 1;
                continue;  // snpCount ++ not works
            }
            snpCount++;
        } // end of non-dup
    } //  end of while loop
    in.close();
    LOGGER << geneInfoVecLD.size() << " genes to be included from [" + infoFile + "]." << endl;
    // Step 2. if genes from eQTL summary statistics are not found in genelist from ld matrix,
    // gene->kept will be set as false.
    for(int i = 0;i < numGenes;i++){
        gene = geneInfoVec[i];
        iterGene = geneInfoMapLD.find(gene->ensemblID);
        if(iterGene == geneInfoMapLD.end()){
            // set LD genelist as benchmark, this BESD gene is not included in LD genelist
            gene->kept = false;
            continue;
        }      
    }
    return geneInfoVecLD.size();
}
void Data::readEigenMatGeneSnpInfoFile(vector<EqtlInfo *> &eqtlInfoVecLD,map<string, EqtlInfo *> &eqtlInfoMapLD, const string &snpInfoFile){
    ifstream in(snpInfoFile.c_str());
    if (!in) LOGGER.e(0, " can not open the file [" + snpInfoFile + "] to read.");
    LOGGER << "Reading xQTL LDM SNP info from file [" + snpInfoFile + "]." << endl;
    // vector<EqtlInfo*> eqtlInfoVecTmp;
    // map<string, EqtlInfo*> eqtlInfoMapTmp;
    EqtlInfo *eqtl;
    map<string, EqtlInfo*>::iterator iterEqtl;
    map<string, SnpInfo*>::iterator iterSnp;
    string header;
    string id, allele1, allele2;
    int chr, physPos,ld_n;
    double genPos;
    double allele1Freq;
    int idx = 0;
    getline(in, header);
    int imputeAf = 0;
    while (in >> chr >> id >> genPos >> physPos >> allele1 >> allele2>>allele1Freq>>ld_n) {
        eqtl = new EqtlInfo(idx++, id, allele1, allele2, chr, genPos, physPos,allele1Freq);
        // Step 2.1 match LD eqtl to BESD eqtl
        iterEqtl = eqtlInfoMap.find(eqtl->rsID);
        if(iterEqtl == eqtlInfoMap.end()){
            // gene ld eqtl cannot be found in besd eqtl. this need to be imputed.
            eqtl->included = false;
        }else{
            eqtl->eqtl_n = iterEqtl->second->eqtl_n;
            if(eqtl->physPos != iterEqtl->second->physPos){
                LOGGER.e(0, "Physical position of eqtl from gene LD [ phypos: " + to_string(eqtl->physPos) +  "] and BESD [ phypos: " + to_string(iterEqtl->second->physPos) +  "] is not consistent. Please check the version of their human genome reference.");
            }
            // allele flip between gene LD eqtl and BESD eqtl;

            // check allele frequency
            if(iterEqtl->second->af == -9 || abs(iterEqtl->second->af ) >= 1){
                if(imputeAf == 0){
                    LOGGER << "Warning: allele frequency in eQTL esi file is missing, use allele frequency from gene eigen matrix info." << endl;
                    imputeAf ++;
                }
            } else {
                if(eqtl->a1 == iterEqtl->second->a1 && eqtl->a2 == iterEqtl->second->a2){
                    eqtl->af = iterEqtl->second->af;
                } else if (eqtl->a2 == iterEqtl->second->a1 && eqtl->a1 == iterEqtl->second->a2) {
                    eqtl->af = 1- iterEqtl->second->af;
                    eqtl->flipped = true;
                }   
            } // end of allele frequency
        }
        eqtl->ld_n  = ld_n;
        eqtlInfoVecLD.push_back(eqtl);
        if (eqtlInfoMapLD.insert(pair<string, EqtlInfo*>(id, eqtl)).second == false) {
            LOGGER.e(0, " Duplicate SNP ID found: \"" + id + "\".");
        }
    }
    in.close();
    LOGGER << eqtlInfoVecLD.size() << " SNPs to be included from [" + snpInfoFile + "]." << endl;
    // remove unmatch eqtl from eQTL summary
    // int imputeAf = 0;
    // for(int i = 0; i < numeQTLs; i++){
    //     eqtl = eqtlInfoVec[i];

    //     /// step 1. remove eqtl without matching to gwas snplist. 
    //     // Already done this in the readEqtlSummaryFile function
    //     // if(snpInfoMap.find(eqtl->rsID) == snpInfoMap.end()){eqtl->included = false; continue;}
    //     // step 2.2 match BESD eqtl to LD eqtl
    //     iterEqtl = eqtlInfoMapLD.find(eqtl->rsID);
    //     if(iterEqtl == eqtlInfoMapLD.end() ){
    //         eqtl->included = false;
    //         continue;
    //     }
    //     if(eqtl->af == -9 || abs(eqtl->af) >= 1){
    //         if(imputeAf == 0){
    //             LOGGER << "Warning: allele frequency in eQTL esi file is missing, use allele frequency from gene eigen matrix info." << endl;
    //             imputeAf ++;
    //         }
    //         if(eqtl->a1 == iterEqtl->second->a1 && eqtl->a2 == iterEqtl->second->a2){
    //                 eqtl->af = iterEqtl->second->af;
    //         } else if (eqtl->a2 == iterEqtl->second->a1 && eqtl->a1 == iterEqtl->second->a2) {eqtl->af = 1- iterEqtl->second->af;}   
    //         }
    //     }
    // return eqtlInfoVecTmp;
}
void Data::readEigenMatGeneBinFile(vector<GeneInfo*> &geneInfoVecLD,map<string, GeneInfo *> &geneInfoMapLD,vector<EqtlInfo*> &eqtlInfoVecLD,map<string, EqtlInfo *> eqtlInfoMapLD,const string &geneEigenMatrixFile, const double eigenCutoff,const double diag_mod){
    FILE *fp = fopen(geneEigenMatrixFile.c_str(), "rb");
    if(!fp){LOGGER.e(0, " can not open the file [" + geneEigenMatrixFile + "] to read.");}
    GeneInfo * gene;  // gene from eqtl summary statistics
    GeneInfo * geneLD;  // gene from per gene ld matrix
    map<string, GeneInfo*>::iterator iterGene;
    map<string, SnpInfo*>::iterator iterSnp;
    map<string, EqtlInfo*>::iterator iterEqtl;
    int numSnpInRegion;
    vector<int> commonSnpIdx;
    vector<string> commonSnpID;
    set<string> eQTLUniqSet; // used to calculate num of eqtls in genic region.
    std::vector<string>::iterator iterVec;
    eQTLUniqSet.clear();
    int totalGeneNum=0;
    //////////////////////////////////////////////////////////////////////////////////
    //////// Step 1. Match gene and eQTLs between BESD and LD gene reference /////////
    //////////////////////////////////////////////////////////////////////////////////
    // #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < geneInfoVecLD.size(); i++){
        //////////////////////////////////////////////////////////////////////////////////
        //////// Step 1.1 read gene LD matrix one by one /////////
        //////////////////////////////////////////////////////////////////////////////////
        geneLD = geneInfoVecLD[i];
        if (!geneLD->kept) continue;

        numSnpInRegion = geneLD->numSnpInGene;

        int32_t cur_mknum = 0;
        int32_t cur_k = 0;
        float sumPosEigVal = 0;
        float oldEigenCutoff =0;

        // 1. marker number
        if(fread(&cur_mknum, sizeof(int32_t), 1, fp) != 1){
            LOGGER.e(0,"Read " + geneEigenMatrixFile + " error (m)");
        }
        if(cur_mknum != numSnpInRegion){
            LOGGER << "Current number: " << cur_mknum << " numSnpInRegion: " << numSnpInRegion << endl;
            LOGGER.e(0,"In region  " + to_string(i) + ", inconsistent marker number to marker information in " + geneEigenMatrixFile);
        }
        // 2. ncol of eigenVec (number of eigenvalues)
        if(fread(&cur_k, sizeof(int32_t), 1, fp) != 1){
            LOGGER.e(0,"In region  " + to_string(i) + ", error about number of eigenvalues in  " + geneEigenMatrixFile);
        }
        // 3. sum of all positive eigenvalues
        if(fread(&sumPosEigVal, sizeof(float), 1, fp) != 1){
            LOGGER.e(0,"In region  " + to_string(i) + ", error about sumPosEigVal in " + geneEigenMatrixFile);
        }
        // 4. eigenCutoff
        if(fread(&oldEigenCutoff, sizeof(float), 1, fp) != 1){
            LOGGER.e(0,"In region  " + to_string(i) + ", error about oldEigenCutoff used in " + geneEigenMatrixFile);
        }
        // 5. eigenvalues
        VectorXf lambda(cur_k);
        if(fread(lambda.data(), sizeof(float), cur_k, fp) != cur_k){
            LOGGER.e(0,"In region  " + to_string(i) + ",size error about eigenvalues in " + geneEigenMatrixFile);
        }
        //6. read eigenvector
        MatrixXf U(cur_mknum, cur_k);
        uint64_t nElements = (uint64_t)cur_mknum * (uint64_t)cur_k;
        if(fread(U.data(), sizeof(float), nElements, fp) != nElements){
            LOGGER << "fread(U.data(), sizeof(double), nElements, fp): " << fread(U.data(), sizeof(float), nElements, fp) << endl;
            LOGGER << "nEle: " << nElements << " U.size: " << U.size() <<  " U.col: " << U.cols() << " row: " << U.rows() << endl;
            LOGGER.e(0,"In region  " + to_string(i) + ",size error about eigenvectors in " + geneEigenMatrixFile);
        }
        // 8. select eigenvectors based on proportion
        MatrixXd eigenVector;
        VectorXd eigenValue;
        bool haveValue = false;
        int revIdx = 0;
        if(oldEigenCutoff != eigenCutoff & i == 0){
            LOGGER << "Warning: current proportion of variance in region  is set as " + to_string(eigenCutoff)+ ". But the proportion of variance is set as "<< to_string(oldEigenCutoff) + " in "  + geneEigenMatrixFile + ".\n";
        }
        if(eigenCutoff < oldEigenCutoff){
            truncateEigenMatrix(sumPosEigVal, eigenCutoff, lambda.cast<double>(), U.cast<double>(), eigenValue, eigenVector);
        } else {
            eigenValue = lambda.cast<double>();
            eigenVector = U.cast<double>();
        }

        //////////////////////////////////////////////////////////////////////////////////
        //////// Step 1.2 gene LD matrix QC  /////////
        //////////////////////////////////////////////////////////////////////////////////
        // remove genes if genes from ld matrix are not found in genelist from eQTL summary statistics
        iterGene = geneInfoMap.find(geneLD->ensemblID);
        if(iterGene != geneInfoMap.end()){
            gene = iterGene->second; // gene from besd
            if(!gene->kept) geneLD->kept = false;
            if(!geneLD->kept) gene->kept = false;
            if (!gene->kept) continue;
        }
        //////////////////////////////////////////////////////////////////////////////////////////
        //////// Step 1.4 Here we use gene eigen matrix as benchmark, exclude redundant Snps and 
        ////////      impute cis-eQTLs if gaps between eQTLs occur.
        //////////////////////////////////////////////////////////////////////////////////////////
        if(geneLD->impCisSnpEffBool){
            Stat::Normal normal;
            /// Step 1. construct LD
            MatrixXd LDPerGene = eigenVector * eigenValue.asDiagonal() * eigenVector.transpose();
            
            LDPerGene.diagonal().array() += (double)diag_mod;
            /// Step 2. Construct the LD correlation matrix among the typed SNPs(LDtt) and the LD correlation matrix among the missing SNPs and typed SNPs (LDit).
            // Step 2.1 divide SNPs into typed and untyped SNPs
            // VectorXi typedSnpIdx(gene->typSnpIdxGeneLD.size());
            // VectorXi untypedSnpIdx(gene->impSnpIdxGeneLD.size());
            unsigned numTypSnpIdxGeneLD = geneLD->typSnpIdxGeneLD.size();
            double localBeta, localSE;
            VectorXd zTypSnp(numTypSnpIdxGeneLD);
            VectorXd nTypSnp(numTypSnpIdxGeneLD);
            VectorXd varyTypSnp(numTypSnpIdxGeneLD);

            for(unsigned j=0;j < numTypSnpIdxGeneLD; j++){
                iterEqtl = eqtlInfoMapLD.find(geneLD->cisSnpNameVec[geneLD->typSnpIdxGeneLD[j]]);
                if(iterEqtl->second->included){
                    int snpBesdIdx = gene->cisSnpID2IdxMapInGene[geneLD->cisSnpNameVec[j]]; // Use this to extract eqtl beta and se from besd format.
                    // typed snp
                    localBeta = gene->eQTLMarginEffect[snpBesdIdx];
                    localSE = gene->eQTLMarginEffectSE[snpBesdIdx];
                    zTypSnp[j] = localBeta/ localSE;
                    nTypSnp[j] = iterEqtl->second->eqtl_n;
                    double hetj = 2.0 * iterEqtl->second->af * (1.0 - iterEqtl->second->af);
                    varyTypSnp[j] = hetj * (iterEqtl->second->eqtl_n * localSE * localSE + localBeta * localBeta);
                }
            }

            // Step 2.2 construct LDtt and LDit and Ztt.
            MatrixXd LDtt = LDPerGene(geneLD->typSnpIdxGeneLD,geneLD->typSnpIdxGeneLD);
            MatrixXd LDit = LDPerGene(geneLD->impSnpIdxGeneLD,geneLD->typSnpIdxGeneLD);
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
            VectorXd eQTLMarginEffImp(geneLD->impSnpIdxGeneLD.size());
            VectorXd eQTLMarginEffSEImp(geneLD->impSnpIdxGeneLD.size());

            for(unsigned j = 0; j < geneLD->impSnpIdxGeneLD.size(); j++){
                iterEqtl = eqtlInfoMapLD.find(geneLD->cisSnpNameVec[geneLD->impSnpIdxGeneLD[j]]);
                double base = sqrt(2.0 * iterEqtl->second->af *(1.0 - iterEqtl->second->af) * (nMedian + zImpSnp[j] * zImpSnp[j]));
                eQTLMarginEffImp(j) = zImpSnp[j] * sqrt(varyMedian)/base;
                eQTLMarginEffSEImp(j) = sqrt(varyMedian) / base;

                iterEqtl->second->eqtl_n = nMedian;
                // eqtl->af = snp->af;
                // snp->gwas_pvalue = 2*(1.0-normal.cdf_01(abs(snp->gwas_b/snp->gwas_se)));
                iterEqtl->second->included = true;
                //LOGGER << "b " << snp->gwas_b << " se " << snp->gwas_se << " z " << snp->gwas_b/snp->gwas_se << " p " << snp->gwas_pvalue << endl;
            }
            geneLD->eQTLMarginEffect.resize(geneLD->numSnpInGene);
            geneLD->eQTLMarginEffectSE.resize(geneLD->numSnpInGene);
            // typed
            geneLD->eQTLMarginEffect(geneLD->typSnpIdxGeneLD) = gene->eQTLMarginEffect(geneLD->typSnpIdxBesd);
            geneLD->eQTLMarginEffectSE(geneLD->typSnpIdxGeneLD) = gene->eQTLMarginEffectSE(geneLD->typSnpIdxBesd);
            // imp
            geneLD->eQTLMarginEffect(geneLD->impSnpIdxGeneLD) =  eQTLMarginEffImp;
            geneLD->eQTLMarginEffect(geneLD->impSnpIdxGeneLD) =  eQTLMarginEffSEImp;
        } else {
            geneLD->eQTLMarginEffect = gene->eQTLMarginEffect(geneLD->typSnpIdxBesd);
            geneLD->eQTLMarginEffectSE = gene->eQTLMarginEffectSE(geneLD->typSnpIdxBesd);
        }
    }
    // find common genes between per LD reference and eQTL summary statistics
    // keptGeneInfoVec  = makeKeptGeneInfoVec(geneInfoVec);
    // numKeptGenes = (unsigned) keptGeneInfoVec.size();
    geneInfoVec = geneInfoVecLD;
    geneInfoMap = geneInfoMapLD;
    eqtlInfoVec = eqtlInfoVecLD;
    eqtlInfoMap = eqtlInfoMapLD;
    keptGeneInfoVec  = makeKeptGeneInfoVec(geneInfoVec);
    numKeptGenes = (unsigned) keptGeneInfoVec.size();
    LOGGER << numKeptGenes  << " genes are matched between gene eigen  LD reference and eQTL summary data." << endl;
    if(numKeptGenes == 0) LOGGER.e(0," there is no gene matched between gene and eigen  LD reference and eQTL summary data. Please check your flag and files.");
    includeMatchedEqtl();
}

/// Generate eigen LD matrix based on gene information file 
MatrixXd Data::generateLDmatrixPerBlock(const vector<int> &snplistIdx){
    int numSnpInRange = snplistIdx.size();
    
    if (numSnpInRange == 0) LOGGER.e(0,"No SNP is retained for analysis.");
    if (numKeptInds == 0) LOGGER.e(0,"No individual is retained for analysis.");
    if (numKeptInds < 2) LOGGER.e(0, " Cannot calculate LD matrix with number of individuals < 2.");
    
    // Here we use the following equation to calculate LD
    // LD = X'*X/N = X'_p * [x_1,..,x_k,..,x_K ]/N
    MatrixXd denseZPZ;
    denseZPZ.setZero(numSnpInRange, numSnpInRange);
    if(numSnpInRange == 1) {
        denseZPZ(0,0) = 1.0;
        return denseZPZ;
    }
    MatrixXd ZP(numSnpInRange, numKeptInds);  // SNP x Ind
    VectorXd Dtmp;
    Dtmp.setZero(numSnpInRange);
    MatrixXd ZZ = Z(all,snplistIdx);
    ZP = ZZ.transpose();
    // Center and scale ZP
    for(unsigned j = 0; j < numSnpInRange; j++){
        Dtmp[j] = Gadget::calcVariance(ZP.row(j))*numKeptInds;
        ZP.row(j) = (ZP.row(j).array() - ZP.row(j).mean())/sqrt(Dtmp[j]);
    }
    // Center and sccale Zk and then calculate denseZPZ
    VectorXd Zk(numKeptInds); // Ind x 1
    Dtmp.setZero(numSnpInRange);
    for(unsigned k = 0; k < numSnpInRange; k++){
        Zk = ZZ.col(k);
        Dtmp[k] = Gadget::calcVariance(Zk)*numKeptInds;
        Zk = (Zk.array() - Zk.mean())/sqrt(Dtmp[k]);
        denseZPZ.col(k) = ZP * Zk;
    }    
    return denseZPZ;
}

void Data::calcGeneLDEigenDecomBasedBesd(const string &besdFile,const string &eqtlSummaryQueryFile, const string &geneListFile, const string specificGeneID, const string &filename, const double cisRegionWind, const double eigenCutoff, bool debugBool ){
    /////////////////////////////////////////
    // Step 1. Map eQTL to gene based on besd fromat 
    /////////////////////////////////////////
    bool smrBesdBool = false;
    bool hasSnpInfo = true;
    if(!besdFile.empty()){
        readEsiFile(besdFile + ".esi",smrBesdBool,hasSnpInfo);  // snp info
        includeMatchedEqtl();
        readEpiFile(besdFile + ".epi"); // probe (gene) info 
        bool makeGeneLDEigenBool = true;
        readBesdFile(besdFile + ".besd",hasSnpInfo, makeGeneLDEigenBool,smrBesdBool);
    } else if(!geneListFile.empty()){
        geneInfoVec.clear();
        geneInfoMap.clear();
        bool geneBool = readMultiEigenMatInfoFile(geneListFile, geneInfoVec,geneInfoMap);
        numGenes = geneInfoVec.size();
    } else if(!eqtlSummaryQueryFile.empty()){
        bool hasGeneLdmInfo = false;
        geneInfoVec.clear(); geneInfoMap.clear();
        eqtlInfoVec.clear(); eqtlInfoMap.clear();
        readQeuryGZFormat(eqtlSummaryQueryFile + ".query.gz",geneInfoVec,geneInfoMap,eqtlInfoVec,eqtlInfoMap,smrBesdBool,hasSnpInfo,hasGeneLdmInfo);
        // summary data first becasue above function don't have summary data
        numeQTLs = eqtlInfoVec.size();
        numGenes = geneInfoVec.size();
    } else {
        LOGGER.e(0," wrong type for gene-snp info when constructing gene LD matrix.");
    }

        // here we only use one gene.
    if(specificGeneID.size()) includeSpecificGene(specificGeneID);

    //////// summarise
    incdEqtlInfoVec = makeIncdEqtlInfoVec(eqtlInfoVec);
    numIncdEqtls = (unsigned) incdEqtlInfoVec.size();
    // resize kept genes
    keptGeneInfoVec  = makeKeptGeneInfoVec(geneInfoVec);
    numKeptGenes = (unsigned) keptGeneInfoVec.size();

    LOGGER  << numKeptGenes  << " Genes with " << numIncdEqtls << " eQTLs are matched between gene eigen  LD reference and eQTL summary data." << endl;
    if(numKeptGenes == 0) LOGGER.e(0, " There is no gene matched between gene and eigen  LD reference and eQTL summary data. Please check your flag and files.");

    // GeneInfo *gene;
    vector<VectorXd> cumsumNonNeg(numKeptGenes);
    vector<string> genelists(numKeptGenes);
    bool readBedBool = true;
    string LDMatType = "gene";
    unsigned currentGeneNum = 0;
    vector<string> cisSnplistInGene;
    SnpInfo * snp;
    string filenameFull = filename + ".eigen." + LDMatType + ".bin";
    FILE *out3 = fopen(filenameFull.c_str(), "wb");
    // debug mode; use to generate ldm
    FILE * rvalOut;
    if(debugBool) {
        string rvalFile = filename + ".ldm." + LDMatType + ".bin";
        rvalOut = fopen(rvalFile.c_str(), "wb");
        LOGGER << "Debug mode is used and LD matrix will be also generated." << endl;
    }
    ///////////////
    float totalSteps = numKeptGenes;
    LOGGER << "The proportion of variance is set as " << eigenCutoff << endl;
    float eigenCutoffFloat = static_cast<float>(eigenCutoff);
    /// read bed file.
    if(debugBool){
        string outFile; std::ofstream file1;
        /// ind snplist
        outFile = ".ldind";
        file1.open((filename + outFile).c_str());
        file1 << "indID" << endl;
        IndInfo *indInfo = NULL;
        for (unsigned i = 0,indIdx =0; i < numInds; i++) {
            indInfo = indInfoVec[i];
            file1 << indInfo->famID << endl; 
        }
        file1.close();
        /// gwas snplist
        outFile = ".ldSnp";
        file1.open((filename + outFile).c_str());
        file1 << "snpID" << endl;
        for(unsigned i =0; i < numIncdSnps;i++){
            SnpInfo * snp = incdSnpInfoVec[i];
            file1 << snp->rsID << endl;
        }
        file1.close();
        outFile = ".ldgeno";
        file1.open((filename + outFile).c_str());
        file1 << Z.format(Eigen::IOFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n"));
        file1.close();
    }
    vector<int> cisSnpIdxPerGene;
    auto startTime = std::chrono::steady_clock::now();
    for (unsigned i = 0; i < numKeptGenes; i++) {
        // LOGGER  << " Generate and save eigen of LD matrix from gene region " << i + 1 << "/" << numKeptGenes << "\r" << flush;  
        string taksName = " Generate and save eigen of LD matrix from gene region";
        Gadget::showProgressBar(i, numKeptGenes, startTime,taksName);

        GeneInfo *gene = keptGeneInfoVec[i]; 
        // here we need to re-order cis-snps
        gene->cisSnpOrderedByGWASLD.clear();
        cisSnpIdxPerGene.clear();
        for(unsigned j = 0; j < numIncdSnps; j++ ){
            snp = incdSnpInfoVec[j];
            if(gene->cisSnpID2IdxMapInGene.find(snp->rsID) != gene->cisSnpID2IdxMapInGene.end()){
                snp->iseQTL = true;
                gene->cisSnpOrderedByGWASLD.push_back(snp->rsID);
                cisSnpIdxPerGene.push_back(j);
            }  
        }
        // if(gene->cisSnpNameVec.size() < 2 || cisSnpIdxPerGene.size() < 2){
        //     LOGGER.w(0,"Gene " + gene->ensemblID + " will be removed due to its xQTL number is less than 2.");
        //     gene->kept = false;
        //     continue;
        // }
        MatrixXf eigenVec;
        VectorXf eigenVal;
        float sumPosEigVal = 0;
        MatrixXf rval = generateLDmatrixPerBlock(cisSnpIdxPerGene).cast<float>();
        int32_t numSnpInBlock = cisSnpIdxPerGene.size();
        if(debugBool) {
            ////////
            string outFile; std::ofstream file1;
            /// gwas snplist
            outFile = ".ldsnp";
            file1.open((filename + gene->ensemblID + outFile).c_str());
            file1 << "snpID" << endl;
            for(unsigned i =0; i < gene->cisSnpOrderedByGWASLD.size();i++){
                file1 << gene->cisSnpOrderedByGWASLD[i] << endl;
            }
            file1.close();
            outFile = ".ldgeno";
            file1.open((filename + gene->ensemblID + outFile).c_str());
            MatrixXd ZZ = Z(all,cisSnpIdxPerGene);
            file1 << ZZ.format(Eigen::IOFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n"));
            file1.close();

            outFile = ".ldgeno.ldm";
            file1.open((filename + gene->ensemblID + outFile).c_str());
            MatrixXd Zcorr = ZZ.transpose() * ZZ;
            file1 << Zcorr.format(Eigen::IOFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n"));
            file1.close();

            outFile = ".ldgeno.ldm.corr";
            file1.open((filename + gene->ensemblID + outFile).c_str());
            file1 << rval.format(Eigen::IOFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n"));
            file1.close();
            //////////////////////////////
            // save number 
            fwrite(&numSnpInBlock, sizeof(int32_t), 1, rvalOut);
            // r correlation matrix
            uint64_t nRval = (uint64_t) numSnpInBlock * (uint64_t) numSnpInBlock;
            fwrite(rval.data(), sizeof(float), nRval,rvalOut );
        }
        eigenDecomposition(rval.cast<float>(), eigenCutoffFloat,eigenVal, eigenVec,sumPosEigVal);
        // save eigen  matrix
        int32_t numEigenValue = eigenVal.size();
        fwrite(&numSnpInBlock, sizeof(int32_t), 1, out3);
        // 2, ncol of eigenVecGene[i]
        fwrite(&numEigenValue, sizeof(int32_t), 1, out3);
        // 3. sum of eigen values
        fwrite(&sumPosEigVal, sizeof(float), 1, out3);
        //4. eigenCutoff
        fwrite(&eigenCutoffFloat, sizeof(float), 1, out3);
        // 5. save eigenvalue
        fwrite(eigenVal.data(), sizeof(float), numEigenValue, out3);
        // 6. save eigenvector;
        uint64_t nElements = (uint64_t) numSnpInBlock * (uint64_t) numEigenValue;
        fwrite(eigenVec.data(), sizeof(float), nElements, out3);
    }
    LOGGER  << "Total "<< numIncdEqtls <<  " xQTLs info from " << numKeptGenes << " moleculars are written into file [" << filenameFull << "]." << endl;
    fclose(out3); 
    if(debugBool){fclose(rvalOut);}
    // here remove genes with xQTL less than 2
    keptGeneInfoVec  = makeKeptGeneInfoVec(geneInfoVec);
    numKeptGenes = (unsigned) keptGeneInfoVec.size();
    genelists.resize(0); // here we don't need to use snplist or genelist,
    
    outputEigenMatIndInfoFromLDMat(LDMatType, filename,false);
}

void Data::calcBlockLDEigenDecom(const string &bedFile, const string &ldBlockInfoFile, int ldBlockRegionWind,const string &filename,const double eigenCutoff, bool outInfoOnly){
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
    if (mapped < 1) LOGGER.e(0, "No GWAS LD block can be mapped to the SNP data. Please check the input data regarding chromosome and bp.");
    else LOGGER << mapped << " GWAS LD block(s) have been mapped to SNP data." << endl;

    keptLdBlockInfoVec = makeKeptLDBlockInfoVec(ldBlockInfoVec);
    numKeptLDBlocks = (unsigned) keptLdBlockInfoVec.size();

    map<string, int>::iterator iter1, iter2;
    map<string, int> snpNameMap;
    vector<int> snpNumInldblock(numKeptLDBlocks);
    vector<VectorXd> cumsumNonNeg(numKeptLDBlocks);
    
    for (i = 0; i < numIncdSnps; i++) {
        snp = incdSnpInfoVec[i];
        snpNameMap.insert(pair<string,int>(snp->rsID, i));
        
    }
    string lastSNPIDInPreLDblocks;
    string LDMatType = "ldblock";
    string filenameFull = filename + ".eigen." + LDMatType + ".bin";
    FILE *out3 = fopen(filenameFull.c_str(), "wb");
    LOGGER << "The proportion of variance is set as " << eigenCutoff << endl;
    float eigenCutoffFloat = static_cast<float>(eigenCutoff);
    vector<string> ldblockNames, snplists;
    bool readBedBool = true;
    for (i = 0; i < numKeptLDBlocks; i++) {
        ldblock = keptLdBlockInfoVec[i]; 
        iter1 = snpNameMap.find(block2snp_1[keptLdBlock2AllLdBlcokMap.at(ldblock->ID ) ]);
        iter2 = snpNameMap.find(block2snp_2[keptLdBlock2AllLdBlcokMap.at(ldblock->ID ) ]);
        bool skip = false;
        if (iter1 == snpNameMap.end() || iter2 == snpNameMap.end() || iter1->second >= iter2->second) ldblock->kept = false;
        snpNumInldblock[i] = iter2->second - iter1->second + 1;
        if(!ldblock->kept) continue;
        vector<int> snp_indx;
        for (j = iter1->second; j <= iter2->second; j++) {
            snp = incdSnpInfoVec[j];
            if(snp->rsID == lastSNPIDInPreLDblocks){continue;} // if snp is matched to two ldblock
            snp_indx.push_back(j);
            ldblock->snpNameVec.push_back(snp->rsID);
            snp->isInBlock = true;
            if (j == iter2->second) { lastSNPIDInPreLDblocks = snp->rsID;}
        }
        /////////////////////////////////
        /////// save LD block info
        LOGGER  << " Generate and save eigen of LD matrix from GWAS LD block region " << i + 1 << "/" << numKeptLDBlocks << "\r" << flush;  
        MatrixXf eigenVec;
        VectorXf eigenVal;
        float sumPosEigVal = 0;
        MatrixXf rval = generateLDmatrixPerBlock(snp_indx).cast<float>();

        eigenDecomposition(rval, eigenCutoffFloat,eigenVal, eigenVec,sumPosEigVal);
        // save eigen  matrix
        int32_t numEigenValue = eigenVal.size();
        // int32_t numSnpInBlock = ldblock->snpNameVec.size();
        int32_t numSnpInBlock = ldblock->snpNameVec.size();
        // 1, nrow of eigenVecGene[i] 
        // cout << "numSnpInBlock: " << numSnpInBlock << endl;
        fwrite(&numSnpInBlock, sizeof(int32_t), 1, out3);
        // 2, ncol of eigenVecGene[i]
        fwrite(&numEigenValue, sizeof(int32_t), 1, out3);
        // 3. sum of eigen values
        fwrite(&sumPosEigVal, sizeof(float), 1, out3);
        //4. eigenCutoff
        fwrite(&eigenCutoffFloat, sizeof(float), 1, out3);
        // 5. save eigenvalue
        fwrite(eigenVal.data(), sizeof(float), numEigenValue, out3);
        // 6. save eigenvector;
        uint64_t nElements = (uint64_t) numSnpInBlock * (uint64_t) numEigenValue;
        fwrite(eigenVec.data(), sizeof(float), nElements, out3);

    }
    LOGGER  << "Total "<< numIncdSnps <<  " SNPs info from " << numKeptLDBlocks << " LD Blocks are written into file [" << out3 << "]." << endl;
    fclose(out3); 
    vector<string> genelists;genelists.resize(0); // here we don't need to use snplist or genelist,
    outputEigenMatIndInfoFromLDMat(LDMatType, filename,false);
}
void Data::calcAndOutputEigenMatBinFromLDMat(const string bedFile, const vector<string> snplists,const string filename,const string LDMatType,const double eigenCutoff){
        string filenameFull = filename + ".eigen." + LDMatType + ".bin";
        FILE *out3 = fopen(filenameFull.c_str(), "wb");
        MatrixXf eigenVec;
        VectorXf eigenVal;
        float sumPosEigVal = 0;
        MatrixXd rval = generateLDmatrixPerBlock(bedFile, snplists);    
        eigenDecomposition(rval.cast<float>(), eigenCutoff,eigenVal, eigenVec,sumPosEigVal);
        // save eigen  matrix
        int32_t numEigenValue = eigenVal.size();
        // int32_t numSnpInBlock = ldblock->snpNameVec.size();
        int32_t numSnpInBlock = snplists.size();
        // 1, nrow of eigenVecGene[i] 
        // cout << "numSnpInBlock: " << numSnpInBlock << endl;
        fwrite(&numSnpInBlock, sizeof(int32_t), 1, out3);
        // 2, ncol of eigenVecGene[i]
        fwrite(&numEigenValue, sizeof(int32_t), 1, out3);
        // 3. sum of eigen values
        fwrite(&sumPosEigVal, sizeof(float), 1, out3);
        //4. eigenCutoff
        fwrite(&eigenCutoff, sizeof(float), 1, out3);
        // 5. save eigenvalue
        fwrite(eigenVal.data(), sizeof(float), numEigenValue, out3);
        // 6. save eigenvector;
        uint64_t nElements = (uint64_t) numSnpInBlock * (uint64_t) numEigenValue;
        fwrite(eigenVec.data(), sizeof(float), nElements, out3);
        fclose(out3); 
}
void Data::outputEigenMatIndInfoFromLDMat(const string LDMatType, const string filename,const bool &mergeBool)  {
    SnpInfo *snp;
    EqtlInfo *eqtl;
    GeneInfo * gene;
    map<string, SnpInfo*>::iterator iterSnp;
    map<string, EqtlInfo*>::iterator iterEqtl;
    map<string, GeneInfo*>::iterator iterGene;
    string outfilename;
    int numSnpsInLDmatType = 0;
    int numLDmatTypes = 0;
    
    if(LDMatType == "gene"){
        outfilename = filename + ".eigen." + LDMatType;
    } else if(LDMatType == "ldblock"){
        outfilename = filename + ".eigen." + LDMatType;
    }
    string outfile1 = outfilename + ".snp.info";
    string outfile2 = outfilename + ".info";
	LOGGER  << incdEqtlInfoVec.size() << endl;
    // write snp info 
    ofstream out1(outfile1.c_str());
    out1 << boost::format("%6s %15s %10s %15s %6s %6s %12s  %10s\n")
    % "Chrom"
    % "ID"
    % "GenPos"
    % "PhysPos"
    % "A1"
    % "A2"
    % "A1Freq"
    % "N";
    if(LDMatType == "ldblock"){
        for (unsigned i=0; i < numIncdSnps ; ++i) {
            snp = incdSnpInfoVec[i];
            if(!snp->isInBlock) continue;
            out1 << boost::format("%6s %15s %10s %15s %6s %6s %12f %10s \n")
                % snp->chrom
                % snp->rsID
                % snp->genPos
                % snp->physPos
                % snp->a1
                % snp->a2
                % snp->af
                % numKeptInds;
            numSnpsInLDmatType ++;
        }
    }else if (LDMatType == "gene"){
        if(!mergeBool){
            // this loop is used to print snpinfo frm LD reference genotypes.
            // when calculating low-rank LD, only blockinfo and eqtllist is used to find
            // corresponding snplists in LD reference
            for (unsigned i=0; i < numIncdSnps ; ++i) {
                snp = incdSnpInfoVec[i];
                if(!snp->iseQTL) continue;
                out1 << boost::format("%6s %15s %10s %15s %6s %6s %12f %10s \n")
                    % snp->chrom
                    % snp->rsID
                    % snp->genPos
                    % snp->physPos
                    % snp->a1
                    % snp->a2
                    % snp->af
                    % numKeptInds;
            numSnpsInLDmatType ++;
            }
        } else {
            /// this loop will be used when mergeing multiple low-rank gene ld.
            /// When merging multiple gene LDs, no snpInfo was constructed
            for (unsigned i=0; i < numIncdEqtls ; ++i) {
                eqtl = incdEqtlInfoVec[i];
                out1 << boost::format("%6s %15s %10s %15s %6s %6s %12f %10s \n")
                    % eqtl->chrom
                    % eqtl->rsID
                    % eqtl->genPos
                    % eqtl->physPos
                    % eqtl->a1
                    % eqtl->a2
                    % eqtl->af
                    % eqtl->ld_n;
                numSnpsInLDmatType ++;
            }     
        } // end of if statement

    }else {

    }
    out1.close();
    LOGGER  << "Total "<< numSnpsInLDmatType <<  " SNPs info from " + LDMatType + " region are written into file [" << outfile1 << "]." << endl;
    ofstream out2(outfile2.c_str());
    if(LDMatType == "gene"){
        // write snp info 
        out2 << boost::format("%6s %15s %10s %15s %15s %15s %15s\n")
        % "Chrom"
        % "probeID"
        % "start"
        % "end"
        % "windows"
        % "snpInGene"
        % "NumSnpInGene";

        if(!mergeBool){
            // when generate gene LD matrix, we need to use variable cisSnpOrderedByGWASLD that reordered based on gwas LD order
            for (unsigned i=0; i < numKeptGenes ; ++i) {
                gene = keptGeneInfoVec[i];
                for(unsigned j = 0; j < gene->cisSnpOrderedByGWASLD.size() ; ++j){
                    out2 << boost::format("%6s %15s %10s %15s %15s %15s %15s \n")
                        % gene->chrom
                        % gene->ensemblID
                        % gene->start
                        % gene->end
                        % gene->windows
                        % gene->cisSnpOrderedByGWASLD[j]
                        % gene->cisSnpOrderedByGWASLD.size();
                }
            numLDmatTypes++;
            }
        } else {
            /// when merging multiple LD, we only need to use cisSnpNameVec for each genes
            for (unsigned i=0; i < numKeptGenes ; ++i) {
                gene = keptGeneInfoVec[i];
                for(unsigned j = 0; j < gene->cisSnpNameVec.size() ; ++j){
                    out2 << boost::format("%6s %15s %10s %15s %15s %15s %15s \n")
                        % gene->chrom
                        % gene->ensemblID
                        % gene->start
                        % gene->end
                        % gene->windows
                        % gene->cisSnpNameVec[j]
                        % gene->cisSnpNameVec.size();
                }
            numLDmatTypes++;
            }
        } // end of if statement
    } else if(LDMatType == "ldblock")  {
        // eigen  matrix for ld blocks here.
        out2 << boost::format("%6s %15s %10s %15s %15s %15s\n")
        % "Chrom"
        % "LDBLOCK"
        % "start"
        % "end"
        % "snpInLdBlock"
        % "NumSnpInLdBlock";
        LDBlockInfo * ldblock;
        for (unsigned i=0; i < numKeptLDBlocks ; ++i) {
            ldblock = keptLdBlockInfoVec[i];
            // LOGGER  << "ldblock id: " << ldblock->ID << endl;
            for(unsigned j = 0; j < ldblock->snpNameVec.size() ; ++j){
                out2 << boost::format("%6s %15s %10s %15s %15s %15s \n")
                    % ldblock->chrom
                    % ldblock->ID
                    % ldblock->startPos
                    % ldblock->endPos
                    % ldblock->snpNameVec[j]
                    % ldblock->snpNameVec.size();
                numLDmatTypes ++;
            }
        }
    }else {
    }
    out2.close();
    LOGGER  << "Total " << numLDmatTypes << " " << LDMatType << " region info are written into file [" << outfile2 << "]." << endl;
}

///////////////////////////////////////////////////////////////////////////////////////
///////////// Step 2.1 read LD matrix of ld blocks info (5 functions)     /////////////
///////////////////////////////////////////////////////////////////////////////////////
void Data::readEigenMatrixGene(const string &infoFile, const double eigenCutoff,const string LDMatType, const unsigned includeChr, const double diag_mod){
    //// Step 1. read xQTL LD snp info firstly and report an error when allele difference between GWAS and xQTL LD occur.
    LOGGER  << "............................" << endl;
    LOGGER << "Reading molecular(e.g., gene/protein) low-rank LD matrices..." << endl;
    LOGGER  << "............................" << endl;
    readEigenMatGeneSnpInfoFile(infoFile + ".eigen.gene.snp.info",includeChr);  // here we need to allele check
    includeMatchedEqtl();
    //// Step 2. read geneinfo
    readEigenMatGeneInfoFile(infoFile + ".eigen.gene.info");  //  this contains eqtl sample size
    keptGeneInfoVec  = makeKeptGeneInfoVec(geneInfoVec);
    numKeptGenes = (unsigned) keptGeneInfoVec.size();
    // Step 3. read xQTL low-rank eigen matrix 
    readEigenMatBinFile(infoFile,eigenCutoff,"gene");
}

void Data::readEqtlSummaryFileFromLD2BESD(const string &besdFile, const string &eqtlSummaryQueryFile, const unsigned includeChr, const double afDiff, const double mafmin, const double mafmax, const double pValueThreshold, const bool imputeN,const double eigenCutoff){           
    //// summary of gene and xQTLs 
    /// now remove gene ld blocks without gene summary stats
    keptGeneInfoVec  = makeKeptGeneInfoVec(geneInfoVec);
    numKeptGenes = keptGeneInfoVec.size();
    incdEqtlInfoVec = makeIncdEqtlInfoVec(eqtlInfoVec);
    numIncdEqtls = incdEqtlInfoVec.size();
    /////////////////////////////////////////////
    vector<GeneInfo *> geneInfoVecBESD;
    map<string, GeneInfo *> geneInfoMapBESD;
    vector<EqtlInfo *> eqtlInfoVecBESD;
    map<string, EqtlInfo *> eqtlInfoMapBESD;
    bool epiBool = false, esiBool = false, besdBool = false;
    bool hasSnpInfo = true;
    bool hasGeneLdmInfo = true;
    bool smrBesdBool = false;
    // Step 4. read BESD format xQTL info to extract b and se
    LOGGER  << "............................" << endl;
    LOGGER << "Reading xQTL summary data..." << endl;
    LOGGER  << "............................" << endl;
    if(!besdFile.empty()){
        epiBool = readMultiEpiFile(besdFile + ".epi", geneInfoVecBESD,geneInfoMapBESD,hasSnpInfo,hasGeneLdmInfo);
        esiBool = readMultiEsiFile(besdFile + ".esi",eqtlInfoVecBESD, eqtlInfoMapBESD,smrBesdBool,hasSnpInfo,hasGeneLdmInfo);
        if(epiBool && esiBool){
            besdBool = readMultiBesdFile(besdFile + ".besd",geneInfoVecBESD,geneInfoMapBESD,eqtlInfoVecBESD, eqtlInfoMapBESD,smrBesdBool,hasSnpInfo,hasGeneLdmInfo);
        } else{
            LOGGER.e(0,"Something is wrong about BESD format.");
        }
    } else if(!eqtlSummaryQueryFile.empty()) {
        readQeuryGZFormat(eqtlSummaryQueryFile + ".query.gz",geneInfoVecBESD,geneInfoMapBESD,eqtlInfoVecBESD,eqtlInfoMapBESD ,smrBesdBool,hasSnpInfo,hasGeneLdmInfo);
    } else {
        LOGGER.e(0,"Wrong format for xqtl information");
    }
    // Step 5. Now we have gene info from xQTL LD and gene info from BESD. we need to align them and then impute missing xQTL effect if needed.
    alignxQTLGeneAndBESDGeneInfo(geneInfoVecBESD,geneInfoMapBESD,eqtlInfoVecBESD, eqtlInfoMapBESD);
}

void Data::alignxQTLGeneAndBESDGeneInfo(vector<GeneInfo *> &geneInfoVecBESD,map<string, GeneInfo *> &geneInfoMapBESD,
        vector<EqtlInfo *> &eqtlInfoVecBESD,map<string, EqtlInfo *> &eqtlInfoMapBESD, double diag_mod){
        // align xQTL effect and se
        GeneInfo * gene, *geneBESD; // gene is info from xQTL LD; geneBESD from besd format file. here we use gene as benchmark
        EqtlInfo * eqtl, *eqtlBESD;
        map<string, GeneInfo*>::iterator iterGene;
        map<string, EqtlInfo*>::iterator iterEqtl;
        string eqtlID;
        vector<int> gene2CisSnpVecBESDLocal;
        int numGeneImp = 0;
        int numGeneWithLargeMissingRate = 0;
        int numMissSNPPerGene = 0;
        int missingXqtlRatePerGene = 0;
        auto startTime = std::chrono::steady_clock::now();
        for(unsigned j = 0; j < numKeptGenes; j++){
            Gadget::showProgressBar(j, numKeptGenes, startTime,"Match xQTL info to GWAS info");
            gene = keptGeneInfoVec[j];
            if(gene->ensemblID == "ENSG00000124532"){
                cout << endl;
            }
            iterGene = geneInfoMapBESD.find(gene->ensemblID);
            if(iterGene == geneInfoMapBESD.end()){
                gene->kept = false; // cannot find LD gene in BESD gene set;
                continue;
            }
            geneBESD = iterGene->second; // found in LD gene set
            // now we need to check xQTL missing situation and reorder xQTL b and se based on gene LD 
            gene->typSnpIdxBesd.clear(); // used to extract eqtl beta and se values;
            gene->typSnpIdxGeneLD.clear(); // used to extract eqtl with beta idx in gene ld 
            gene->impSnpIdxGeneLD.clear(); // used to extract eqtl idx need to impute
            gene2CisSnpVecBESDLocal.clear();
            numMissSNPPerGene = 0;
            for(unsigned i = 0; i < gene->cisSnpNameVec.size();i ++)
            {
                eqtlID = gene->cisSnpNameVec[i];
                if(geneBESD->cisSnpID2IdxMapInGene.find(eqtlID) == geneBESD->cisSnpID2IdxMapInGene.end()){
                    // missing;
                    gene->impSnpIdxGeneLD.push_back(i);
                    gene->impCisSnpEffBool = true;
                    numMissSNPPerGene ++;
                } else {
                    // found
                    // extract se and b index from besd data 
                    int eqtlIdxBESD = geneBESD->cisSnpID2IdxMapInGene.find(eqtlID)->second;
                    gene2CisSnpVecBESDLocal.push_back(eqtlIdxBESD);
                    gene->typSnpIdxGeneLD.push_back(i);
                    // extract sample size and allele frequency to besd;
                    eqtl = eqtlInfoMap.find(eqtlID)->second; // find eqtl in xqtl ld;
                    eqtlBESD = eqtlInfoMapBESD.find(eqtlID)->second; // find eqtl in besd;
                    if(eqtlBESD->eqtl_n + 999 < 1e-6 ) eqtl->eqtl_n = eqtlBESD->eqtl_n; // use sample size from besd format
                    if(eqtlBESD->af > 0 ) eqtl->af = eqtlBESD->af; // use af from besd file 
                    eqtl->flipped = eqtlBESD->flipped; // use flipped status from besd file
                }
            } // end of loop cisSnpNameVec
            /// calculate missing rate per gene
            if(gene->impCisSnpEffBool) missingXqtlRatePerGene = (double) numMissSNPPerGene/ (double) gene->cisSnpNameVec.size();
            if( missingXqtlRatePerGene > MISSING_XQTL_RATE_PER_GENE){
                gene->kept = false; // too many xQTLs are missing, remove this gene;
                numGeneWithLargeMissingRate ++;
                continue;
            }
            if(gene->impCisSnpEffBool){
                imputexQTLEffect(gene,geneBESD,eqtlInfoMapBESD,gene2CisSnpVecBESDLocal,diag_mod);
                numGeneImp++;
            } else {
                gene->eQTLMarginEffect = geneBESD->eQTLMarginEffect(gene2CisSnpVecBESDLocal);
                gene->eQTLMarginEffectSE = geneBESD->eQTLMarginEffectSE(gene2CisSnpVecBESDLocal);
            } 
            /// add additional sample size
            if(geneBESD->cisSnpSampleSizeMap.size()!=0) {
                gene->cisSnpSampleSizeMap = geneBESD->cisSnpSampleSizeMap;
            }
        } // end of gene->impCisSnpEffBool
        // summary the 
        keptGeneInfoVec  = makeKeptGeneInfoVec(geneInfoVec);
        numKeptGenes = keptGeneInfoVec.size();
        // incdEqtlInfoVec = makeIncdEqtlInfoVec(eqtlInfoVec);
        // numIncdEqtls = incdEqtlInfoVec.size();
        makeIncdEqtlInfoVecBasedOnGeneInfoVec();



        if(numGeneWithLargeMissingRate) LOGGER << to_string(numGeneWithLargeMissingRate) + " gene(s) are removed due to larger missing xQTL rate ( >=" + to_string(MISSING_XQTL_RATE_PER_GENE) + ")" << endl;
        if(numGeneImp) LOGGER << "Imputate " + to_string(numGeneImp) + " gene(s) with rate  ( <" + to_string(MISSING_XQTL_RATE_PER_GENE) + ")" << endl;
        LOGGER << to_string(numKeptGenes) + " genes and " + to_string(numIncdEqtls) + " xQTLs are left after matching xQTL sumstats to molecualr LD reference and checking xQTL missing rate for each gene." << endl;
    }

void Data::imputexQTLEffect(GeneInfo * &gene,GeneInfo * &geneBESD,map<string, EqtlInfo *> &eqtlInfoMapBESD,vector<int> gene2CisSnpVecBESDLocal, double diag_mod) {
    // GeneInfo * gene, *geneBESD; // gene is info from xQTL LD
    EqtlInfo * eqtl, *eqtlBESD;
    string eqtlID;
    map<string, EqtlInfo*>::iterator iterEqtl;
    Stat::Normal normal;
    /// Step 1. construct LD
    MatrixXd LDPerGene = gene->eigenVector * gene->eigenValue.asDiagonal() * gene->eigenVector.transpose();
    LDPerGene.diagonal().array() += (double)diag_mod;
    /// Step 2. Construct the LD correlation matrix among the typed SNPs(LDtt) and the LD correlation matrix among the missing SNPs and typed SNPs (LDit).
    // Step 2.1 divide SNPs into typed and untyped SNPs
    // VectorXi typedSnpIdx(gene->typSnpIdxGeneLD.size());
    // VectorXi untypedSnpIdx(gene->impSnpIdxGeneLD.size());
    unsigned numTypSnpIdxGeneLD = gene->typSnpIdxGeneLD.size();
    double localBeta, localSE;
    VectorXd zTypSnp(numTypSnpIdxGeneLD);
    VectorXd nTypSnp(numTypSnpIdxGeneLD);
    VectorXd varyTypSnp(numTypSnpIdxGeneLD);
    
    for(unsigned j=0;j < numTypSnpIdxGeneLD; j++){
        eqtlID = gene->cisSnpNameVec[gene->typSnpIdxGeneLD[j]];
        iterEqtl = eqtlInfoMapBESD.find(eqtlID); // eqtl info from besd
        eqtlBESD = iterEqtl->second;
        if(eqtlBESD->included){
            int snpBesdIdx = geneBESD->cisSnpID2IdxMapInGene[eqtlBESD->rsID]; // Use this to extract eqtl beta and se from besd format.
            // typed snp
            localBeta = gene->eQTLMarginEffect[snpBesdIdx];
            localSE = gene->eQTLMarginEffectSE[snpBesdIdx];
            zTypSnp[j] = localBeta/ localSE;
            nTypSnp[j] = eqtlBESD->eqtl_n;
            double hetj = 2.0 * eqtlBESD->af * (1.0 - eqtlBESD->af);
            varyTypSnp[j] = hetj * (eqtlBESD->eqtl_n * localSE * localSE + localBeta * localBeta);
        }
    }

    // Step 2.2 construct LDtt and LDit and Ztt.
    MatrixXd LDtt = LDPerGene(gene->typSnpIdxGeneLD,gene->typSnpIdxGeneLD);
    MatrixXd LDit = LDPerGene(gene->impSnpIdxGeneLD,gene->typSnpIdxGeneLD);
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
    VectorXd eQTLMarginEffImp(gene->impSnpIdxGeneLD.size());
    VectorXd eQTLMarginEffSEImp(gene->impSnpIdxGeneLD.size());

    for(unsigned j = 0; j < gene->impSnpIdxGeneLD.size(); j++){
        eqtlID = gene->cisSnpNameVec[gene->impSnpIdxGeneLD[j]];
        eqtl = eqtlInfoMap.find(eqtlID)->second;
        double base = sqrt(2.0 * eqtl->af *(1.0 - eqtl->af) * (nMedian + zImpSnp[j] * zImpSnp[j]));
        eQTLMarginEffImp(j) = zImpSnp[j] * sqrt(varyMedian)/base;
        eQTLMarginEffSEImp(j) = sqrt(varyMedian) / base;
        eqtl->eqtl_n = nMedian;
        eqtl->included = true;
        //LOGGER << "b " << snp->gwas_b << " se " << snp->gwas_se << " z " << snp->gwas_b/snp->gwas_se << " p " << snp->gwas_pvalue << endl;
    }
    gene->eQTLMarginEffect.resize(gene->numSnpInGene);
    gene->eQTLMarginEffectSE.resize(gene->numSnpInGene);
    // typed
    gene->eQTLMarginEffect(gene->typSnpIdxGeneLD) = gene->eQTLMarginEffect(gene2CisSnpVecBESDLocal);
    gene->eQTLMarginEffectSE(gene->typSnpIdxGeneLD) = gene->eQTLMarginEffectSE(gene2CisSnpVecBESDLocal);
    // imp
    gene->eQTLMarginEffect(gene->impSnpIdxGeneLD) =  eQTLMarginEffImp;
    gene->eQTLMarginEffect(gene->impSnpIdxGeneLD) =  eQTLMarginEffSEImp;
}


void Data::cleanUpUselessParameters(){
    eQTLEffMeanAcrossGenes.resize(0);
    sampleSizeAcrossGene.clear();
    eigenValGene.clear();
    eigenVecGene.clear();
    eigenValLdBlock.clear();
    eigenVecLdBlock.clear();
}

long Data::estimateSamSize(VectorXd &beta,VectorXd &se)
{
    long sampleSize;
    VectorXd Y(se.size());
    MatrixXd X(beta.size(),2);
    for(int i=0;i<se.size();i++){
        Y(i)=se[i]*se[i];
        X(i,0)=1.0;
        X(i,1)=beta[i]*beta[i]-Y(i);
    }  
    MatrixXd XtX_i=(X.transpose()*X).inverse();
    VectorXd w_hat=XtX_i*X.transpose()*Y;
    // cout << "what : " << w_hat << endl;
    sampleSize = -1/w_hat[1];
    return sampleSize;
}

