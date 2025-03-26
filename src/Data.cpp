//
//  data.cpp
//  gctb
//
//  Created by Jian Zeng on 14/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#include "Data.hpp"

// most read file methods are adopted from GCTA with modification

bool SnpInfo::isProximal(const SnpInfo &snp2, const double genWindow) const {
    return chrom == snp2.chrom && fabs(genPos - snp2.genPos) < genWindow;
}

bool SnpInfo::isProximal(const SnpInfo &snp2, const unsigned physWindow) const {
    return chrom == snp2.chrom && abs(physPos - snp2.physPos) < physWindow;
}

void AnnoInfo::getSnpInfo() {
    vector<SnpInfo*> incdSnpVec;
    unsigned numIncdSnps = 0;
    SnpInfo *snp;
    for (unsigned i=0; i<size; ++i) {
        snp = memberSnpVec[i];
        if (snp->included) {
            incdSnpVec.push_back(snp);
            snp->annoIdx.resize(snp->numAnnos);   // the index of SNP in the annotation
            for (unsigned j=0; j<snp->numAnnos; ++j) {
                if (snp->annoVec[j] == this) {
                    snp->annoIdx[j] = numIncdSnps;
                }
            }
            ++numIncdSnps;
        }
    }
    memberSnpVec = incdSnpVec;
    size = numIncdSnps;
    snp2pq.resize(size);
    for (unsigned i=0; i<size; ++i) {
        snp2pq[i] = memberSnpVec[i]->twopq;
    }
}

void AnnoInfo::print() {
    LOGGER << boost::format("%6s %30s %12s %8.3f\n")
    % (idx+1)
    % label
    % size
    % fraction;
}

void Data::readFamFile(const string &famFile){
    // ignore phenotype column
    ifstream in(famFile.c_str());
    if (!in)  LOGGER.e(0,"Can not open the file [" + famFile + "] to read.");
    LOGGER.i(0, "Reading PLINK FAM file from [" + famFile + "].");
    indInfoVec.clear();
    indInfoMap.clear();
    string fid, pid, dad, mom, sex, phen;
    unsigned idx = 0;
    while (in >> fid >> pid >> dad >> mom >> sex >> phen) {
        IndInfo *ind = new IndInfo(idx++, fid, pid, dad, mom, atoi(sex.c_str()));
        indInfoVec.push_back(ind);
        if (indInfoMap.insert(pair<string, IndInfo*>(ind->catID, ind)).second == false) {
            LOGGER.e(0,"Duplicate individual ID found: \"" + fid + "\t" + pid + "\".");
        }
    }
    in.close();
    numInds = (unsigned) indInfoVec.size();
    LOGGER.i(0, to_string(numInds) + " individuals to be included from [" + famFile + "].");
}

void Data::readBimFile(const string &bimFile) {
    // Read bim file: recombination rate is defined between SNP i and SNP i-1
    ifstream in(bimFile.c_str());
    if (!in) LOGGER.e(0,"Can not open the file [" + bimFile + "] to read.");
    LOGGER << "Reading PLINK BIM file from [" + bimFile + "]." << endl;
    snpInfoVec.clear();
    snpInfoMap.clear();
    string id, allele1, allele2;
    unsigned chr, physPos;
    double genPos;
    unsigned idx = 0;
    while (in >> chr >> id >> genPos >> physPos >> allele1 >> allele2) {
        SnpInfo *snp = new SnpInfo(idx++, id, allele1, allele2, chr, genPos, physPos);
        snpInfoVec.push_back(snp);
        chromosomes.insert(snp->chrom);
        if (snpInfoMap.insert(pair<string, SnpInfo*>(id, snp)).second == false) {
            LOGGER.e(0,"Duplicate SNP ID found: \"" + id + "\".");
        }
    }
    in.close();
    numSnps = (unsigned) snpInfoVec.size();
    LOGGER << numSnps << " SNPs to be included from [" + bimFile + "]." << endl;
}


void Data::readBedFile(const bool noscale, const string &bedFile){
    unsigned i = 0, j = 0;
    if (numIncdSnps == 0) LOGGER.e(0,"Error: No SNP is retained for analysis.");
    if (numKeptInds == 0) LOGGER.e(0,"Error: No individual is retained for analysis.");

    Z.resize(numKeptInds, numIncdSnps);
    ZPZdiag.resize(numIncdSnps);
    snp2pq.resize(numIncdSnps);

    // /////////// debug param ////////////
    // Zori.resize(numKeptInds, numIncdSnps);
    // Zcen.resize(numKeptInds, numIncdSnps);
    // ////////////////////////////////////

    // Read bed file
    FILE *in = fopen(bedFile.c_str(), "rb");
    if (!in) LOGGER.e(0,"Error: can not open the file [" + bedFile + "] to read.");
    cout << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
    char header[3];
    fread(header, sizeof(header), 1, in);
    if (!in || header[0] != 0x6c || header[1] != 0x1b || header[2] != 0x01) {
        cerr << "Error: Incorrect first three bytes of bed file: " << bedFile << endl;
        exit(1);
    }

    // Read genotypes
    SnpInfo *snpInfo = NULL;
    IndInfo *indInfo = NULL;
    unsigned snp = 0;
    unsigned nmiss=0;
    double sum=0.0, mean=0.0;

    const int bedToGeno[4] = {2, -9, 1, 0};
    unsigned size = (numInds+3)>>2;
    int genoValue;
    unsigned long long skip = 0;

    for (j = 0, snp = 0; j < numSnps; j++) {  // code adopted from BOLT-LMM with modification
        snpInfo = snpInfoVec[j];
        sum = 0.0;
        nmiss = 0;

        if (!snpInfo->included) {
//            in.ignore(size);
            skip += size;
            continue;
        }
        if (skip) fseek(in, skip, SEEK_CUR);
        skip = 0;
        char *bedLineIn = new char[size];
        fread(bedLineIn, 1, size, in);
        for (i = 0; i < numInds; i++) {
            indInfo = indInfoVec[i];
            if (!indInfo->kept) continue;
            genoValue = bedToGeno[(bedLineIn[i>>2]>>((i&3)<<1))&3];
            Z(indInfo->index, snp) = genoValue;
            if (genoValue == -9) ++nmiss;   // missing genotype
            else sum += genoValue;
        }
        delete[] bedLineIn;    
        // fill missing values with the mean
        mean = sum/double(numKeptInds - nmiss);
        if (nmiss) {
            for (i=0; i<numKeptInds; ++i) {
                if (Z(i,snp) == -9) Z(i,snp) = mean;
            }
        }
        // compute allele frequency
        snpInfo->af = 0.5f*mean;
        snp2pq[snp] = snpInfo->twopq = 2.0f*snpInfo->af*(1.0-snpInfo->af);
        // Zori.col(snp) = Z.col(snp);
        Z.col(snp).array() -= mean; // center column by 2p rather than the real mean
        // Zcen.col(snp) = Z.col(snp);
        if (!noscale) {
            Z.col(snp).array() /= sqrt(snp2pq[snp]);  // standardise to have variance one
            snp2pq[snp] = Gadget::calcVariance(Z.col(snp));
        }
        
        Z.col(snp).array() *= RinverseSqrt.array();

        if (++snp == numIncdSnps) break;
    }
    fclose(in);
    ZPZdiag = Z.colwise().squaredNorm();
    cout << "snp2pq sum: " << snp2pq.mean() << endl;
    
    cout << "Genotype data for " << numKeptInds << " individuals and " << numIncdSnps << " SNPs are included from [" + bedFile + "]." << endl;
}

void Data::readPhenotypeFile(const string &phenFile, const unsigned mphen) {
    if (phenFile.empty()) return;
    // NA: missing phenotype
    ifstream in(phenFile.c_str());
    if (!in) LOGGER.e(0,"Can not open the phenotype file [" + phenFile + "] to read.");
    LOGGER << "Reading phenotypes from [" + phenFile + "]." << endl;
    map<string, IndInfo*>::iterator it, end=indInfoMap.end();
    IndInfo *ind = NULL;
    Gadget::Tokenizer colData;
    string inputStr;
    string sep(" \t");
    string id;
    unsigned line=0;
    while (getline(in,inputStr)) {
        if (inputStr.empty()) {continue;}
        colData.getTokens(inputStr, sep);
        id = colData[0] + ":" + colData[1];
        it = indInfoMap.find(id);
        if (it != end && colData[mphen+1] != "NA" && colData[mphen+1] != "-9") {
            ind = it->second;
            ind->phenotype = std::stod(colData[mphen+1].c_str());
            ++line;
            // cout << "line: " << line << " id: " << id  <<  " " << ind->phenotype  << endl;
        }
    }
    in.close();
    LOGGER << "Non-missing phenotypes of trait " << mphen << " of " << line << " individuals are included from [" + phenFile + "]." << endl;
}

void Data::readCovariateFile(const string &covarFile){
    if (covarFile.empty()) return;
    ifstream in(covarFile.c_str());
    if (!in) LOGGER.e(0,"Can not open the file [" + covarFile + "] to read.");
    map<string, IndInfo*>::iterator it, end=indInfoMap.end();
    IndInfo *ind = NULL;
    Gadget::Tokenizer colData;
    string inputStr;
    string sep(" \t");
    string id;
    unsigned line=0;
    unsigned numCovariates=0;
    while (getline(in,inputStr)) {
        if (inputStr.empty()) {continue;}
        colData.getTokens(inputStr, sep);
        if (line==0) {
            numCovariates = (unsigned)colData.size() - 2;
            numFixedEffects = numCovariates + 1;
            fixedEffectNames.resize(numFixedEffects);
            fixedEffectNames[0] = "Intercept";
            for (unsigned i=0; i<numCovariates; ++i) {
                fixedEffectNames[i+1] = colData[i+2];
            }
            ++line;
        }
        id = colData[0] + ":" + colData[1];
        it = indInfoMap.find(id);
        if (it != end) {
            ind = it->second;
            ind->covariates.resize(numCovariates + 1);  // plus intercept
            ind->covariates[0] = 1;
            for (unsigned i=2; i<colData.size(); ++i) {
                ind->covariates[i-1] = std::stod(colData[i].c_str());
            }
            ++line;
        }
    }
    in.close();
    
    LOGGER << "Read " << numCovariates << " covariates from [" + covarFile + "]." << endl;
}

void Data::readRandomCovariateFile(const string &covarFile){
    if (covarFile.empty()) return;
    ifstream in(covarFile.c_str());
    if (!in) LOGGER.e(0,"Can not open the file [" + covarFile + "] to read.");
    map<string, IndInfo*>::iterator it, end=indInfoMap.end();
    IndInfo *ind = NULL;
    Gadget::Tokenizer colData;
    string inputStr;
    string sep(" \t");
    string id;
    unsigned line=0;
    while (getline(in,inputStr)) {
        colData.getTokens(inputStr, sep);
        if (line==0) {
            numRandomEffects = (unsigned)colData.size() - 2;
            randomEffectNames.resize(numRandomEffects);
            for (unsigned i=0; i<numRandomEffects; ++i) {
                randomEffectNames[i] = colData[i+2];
            }
            ++line;
        }
        id = colData[0] + ":" + colData[1];
        it = indInfoMap.find(id);
        if (it != end) {
            ind = it->second;
            ind->randomCovariates.resize(numRandomEffects);
            for (unsigned i=0; i<numRandomEffects; ++i) {
                ind->randomCovariates[i] = std::stod(colData[i+2].c_str());
            }
            ++line;
        }
    }
    in.close();
    
    LOGGER << "Read " << numRandomEffects << " covariates as random effects from [" + covarFile + "]." << endl;
}

void Data::readResidualDiagFile(const string &resDiagFile){
    if (resDiagFile.empty()) return;
    ifstream in(resDiagFile.c_str());
    if (!in) LOGGER.e(0,"Can not open the file [" + resDiagFile + "] to read.");
    map<string, IndInfo*>::iterator it, end=indInfoMap.end();
    IndInfo *ind = NULL;
    string fid, pid, resDiag;
    string id;
    unsigned line=0;
    while (in >> fid >> pid >> resDiag) {
        id = fid + ":" + pid;
        it = indInfoMap.find(id);
        if (it != end) {
            ind = it->second;
            ind->rinverse = 1.0/std::stod(resDiag.c_str());
            ++line;
        }
    }
    in.close();
    weightedRes = true;
    
    LOGGER << "Read residual diagonal values for " << line << " individuals from [" + resDiagFile + "]." << endl;
}

void Data::keepMatchedInd(const string &keepIndFile, const unsigned keepIndMax,bool ldBool){  // keepIndFile is optional
    map<string, IndInfo*>::iterator it, end=indInfoMap.end();
    IndInfo *ind = NULL;
    vector<string> keep;
    keep.reserve(numInds);
    unsigned cnt=0;

    if(!ldBool){
        if (numFixedEffects) {
            for (unsigned i=0; i<numInds; ++i) {
                ind = indInfoVec[i];
                ind->kept = false;
                if (numRandomEffects) {
                    if (ind->phenotype!=-9 && ind->covariates.size() && ind->randomCovariates.size()) {
                        if (keepIndMax > cnt++)
                            keep.push_back(ind->catID);
                    }
                }
                else {
                    if (ind->phenotype!=-9 && ind->covariates.size()) {
                        if (keepIndMax > cnt++)
                            keep.push_back(ind->catID);
                    }
                }
            }
        }
        else {
            numFixedEffects = 1;
            fixedEffectNames = {"Intercept"};
            for (unsigned i=0; i<numInds; ++i) {
                ind = indInfoVec[i];
                ind->kept = false;
                if (numRandomEffects) {
                    if (ind->phenotype!=-9 && ind->randomCovariates.size()) {
                        ind->covariates.resize(1);
                        ind->covariates << 1;
                        if (keepIndMax > cnt++)
                            keep.push_back(ind->catID);
                    }
                }
                else {
                    if (ind->phenotype!=-9) {
                        ind->covariates.resize(1);
                        ind->covariates << 1;
                        if (keepIndMax > cnt++)
                            keep.push_back(ind->catID);
                    }
                }
            }
        }
    }

    if (!keepIndFile.empty()) {
        ifstream in(keepIndFile.c_str());
        if (!in) LOGGER.e(0,"Can not open the file [" + keepIndFile + "] to read.");
        LOGGER.i(0, "Reading kept individuals with GWAS phenotype from [" + keepIndFile + "].");
        string fid, pid;
        keep.clear();
        while (in >> fid >> pid) {
            keep.push_back(fid + ":" + pid);
        }
        in.close();
        // set false for all inds
        for (unsigned i=0; i<numInds; ++i) {
            ind = indInfoVec[i];
            ind->kept = false;
        }

    }
    
    numKeptInds = 0;
    for (unsigned i=0; i<keep.size(); ++i) {
        it = indInfoMap.find(keep[i]);
        if (it == end) {
            Gadget::Tokenizer token;
            token.getTokens(keep[i], ":");
            LOGGER.w(0," Individual " + token[0] + " " + token[1] + " from file [" + keepIndFile + "] does not exist!");
        } else {
            ind = it->second;
            ind->kept = true;
        }
    }
    

    keptIndInfoVec = makeKeptIndInfoVec(indInfoVec);
    numKeptInds =  (unsigned) keptIndInfoVec.size();

    if (numKeptInds == 0) LOGGER.e(0,"No individual is retained for analysis, please check individual difference between phenotype and keep file" + keepIndFile);
    

    RinverseSqrt.setOnes(numKeptInds); // used when reading bed format
    Rsqrt.setOnes(numKeptInds);
    for (unsigned i=0; i<numKeptInds; ++i) {
        RinverseSqrt[i] = sqrt(keptIndInfoVec[i]->rinverse);
        Rsqrt[i] = 1.0/RinverseSqrt[i];
    }

    if(ldBool){
        LOGGER << numKeptInds << " matched individuals are kept in the GWAS phenotype." << endl;
        return;
    }
    
    y.setZero(numKeptInds);
    for (unsigned i=0; i<numKeptInds; ++i) {
        y[i] = keptIndInfoVec[i]->phenotype;
    }
    y.array() *= RinverseSqrt.array();
    ypy = (y.array()-y.mean()).square().sum();
    // varPhenotypic = ypy/(numKeptInds - 1);
    varPhenotypic = Gadget::calcVariance(y);

    if(numFixedEffects){
        X.resize(numKeptInds, numFixedEffects);
        for (unsigned i=0; i<numKeptInds; ++i) {
            X.row(i) = keptIndInfoVec[i]->covariates.array() * RinverseSqrt[i];
        }
        XPXdiag = X.colwise().squaredNorm();
    }
    
    if (numRandomEffects) {
        W.resize(numKeptInds, numRandomEffects);
        for (unsigned i=0; i<numKeptInds; ++i) {
            W.row(i) = keptIndInfoVec[i]->randomCovariates.array() * RinverseSqrt[i];
        }
        WPWdiag = W.colwise().squaredNorm();
    }
    
    LOGGER << numKeptInds << " matched individuals are kept in the GWAS phenotype." << endl;
}

void Data::initVariances(const double heritability, const double propVarRandom){
//    double varPhenotypic = ypy/numKeptInds;;
    varGenotypic = varPhenotypic * heritability;
    varResidual  = varPhenotypic - varGenotypic;
    varRandom    = varPhenotypic * propVarRandom;
    // LOGGER <<ypy<<" "<<numKeptInds<<" "<< " vary: " << varPhenotypic<<" varg: " <<varGenotypic << " vare: " <<varResidual << " " << varRandom << endl;
}

void Data::includeSnp(const string &includeSnpFile){
    ifstream in(includeSnpFile.c_str());
    if (!in) LOGGER.e(0,"Can not open the file [" + includeSnpFile + "] to read.");
    for (unsigned i=0; i<numSnps; ++i) {
        snpInfoVec[i]->included = false;
    }
    map<string, SnpInfo*>::iterator it, end = snpInfoMap.end();
    string id;
    while (in >> id) {
        it = snpInfoMap.find(id);
        if (it != end) {
            it->second->included = true;
        }
    }
    in.close();
}

void Data::excludeSnp(const string &excludeSnpFile){
    ifstream in(excludeSnpFile.c_str());
    if (!in) LOGGER.e(0,"Can not open the file [" + excludeSnpFile + "] to read.");
    map<string, SnpInfo*>::iterator it, end = snpInfoMap.end();
    string id;
    while (in >> id) {
        it = snpInfoMap.find(id);
        if (it != end) {
            it->second->included = false;
        }
    }
    in.close();
}

void Data::includeChr(const unsigned chr){
    if (!chr) return;
    for (unsigned i=0; i<numSnps; ++i){
        SnpInfo *snp = snpInfoVec[i];
        if (snp->chrom != chr) snp->included = false;
    }
}
void Data::includeChreQTL(const unsigned chr){
    if (!chr) return;
    for (unsigned i=0; i<numeQTLs; ++i){
        EqtlInfo *eqtl = eqtlInfoVec[i];
        if (eqtl->chrom != chr) eqtl->included = false;
    }
    LOGGER << "xQTLs from chromosomes [" + to_string(chr) + " ] are kept." << endl;
}
void Data::includeChrGene(const unsigned chr){
    if (!chr) return;
    for (unsigned i=0; i< numGenes; ++i){
        GeneInfo *gene = geneInfoVec[i];
        if (gene->chrom != chr) gene->kept = false;
    }
}


void Data::includeBlockID(const string includeBlockID){
    if (includeBlockID.empty()) return;
    unsigned keptNumBlock = 0;
    for (unsigned i=0; i<numLDBlocks; ++i) {
        LDBlockInfo *blockInfo = ldBlockInfoVec[i];
        if(blockInfo->ID == includeBlockID){
            keptNumBlock++;
        } else {
            blockInfo->kept = false;
            for(unsigned j =0; j < blockInfo->snpInfoVec.size();j++){
                SnpInfo *snp = blockInfo->snpInfoVec[j];
                snp->included = false;
            }
        }
    }
    if(keptNumBlock){
        LOGGER << keptNumBlock << " LD Block(s) was kept. " << endl;
    }else {
        LOGGER.e(0,"No LD Block was left.");
    }
}

void Data::includeBlock(const unsigned block){
    if (!block) return;
    if (block > ldBlockInfoVec.size()) LOGGER.e(0," Reqest to include block " + to_string(block) + " but there are only " + to_string(ldBlockInfoVec.size()) + " in total!");
    for (unsigned i=0; i<numLDBlocks; ++i) {
        LDBlockInfo *blockInfo = ldBlockInfoVec[i];
        if (block != i+1) blockInfo->kept = false;
    }
    LDBlockInfo *blockInfo = ldBlockInfoVec[block-1];
    unsigned cnt = 0;
        
    if (blockInfo->numSnpInBlock) {
        for (unsigned i=0; i<numSnps; ++i){
            SnpInfo *snpInfo = snpInfoVec[i];
            snpInfo->included = false;
        }
        for (unsigned i=0; i<blockInfo->numSnpInBlock; ++i) {
            SnpInfo *snpInfo = blockInfo->snpInfoVec[i];
            snpInfo->included = true;
            ++cnt;
        }
    }
    else {
        for (unsigned i=0; i<numSnps; ++i){
            SnpInfo *snpInfo = snpInfoVec[i];
            if(!snpInfo->included) continue;
            if (snpInfo->chrom != blockInfo->chrom) {
                snpInfo->included = false;
            } else if (snpInfo->physPos < blockInfo->startPos) {
                snpInfo->included = false;
                // cout << "start here" << endl;
            }else if (snpInfo->physPos > blockInfo->endPos) {
                snpInfo->included = false;
                // cout << "end here" << endl;
            }else {
                snpInfo->included = true;
                ++cnt;
            }
        }
    }
    LOGGER << "Included " << cnt << " SNPs in block " << blockInfo->ID << endl;
}

void Data::includeSkeletonSnp(const string &skeletonSnpFile){
    ifstream in(skeletonSnpFile.c_str());
    if (!in) LOGGER.e(0,"can not open the file [" + skeletonSnpFile + "] to read.");
    map<string, SnpInfo*>::iterator it, end = snpInfoMap.end();
    string id;
    SnpInfo *snp;
    numSkeletonSnps = 0;
    while (in >> id) {
        it = snpInfoMap.find(id);
        if (it != end) {
            snp = it->second;
            snp->included = true;
            snp->skeleton = true;
            ++numSkeletonSnps;
        }
    }
    LOGGER << numSkeletonSnps << " skeleton SNPs are included." << endl;
    in.close();
}

void Data::excludeMHC(){
    long cnt = 0;
    for (unsigned i=0; i<numSnps; ++i) {
        SnpInfo *snp = snpInfoVec[i];
        if (snp->chrom == 6) {
            if (snp->physPos > 28e6 && snp->physPos < 34e6) {
                snp->included = false;
                ++cnt;
            }
        }
    }
    LOGGER << cnt << " SNPs in the MHC region (Chr6:28-34Mb) are excluded." << endl;
}

void Data::excludeAmbiguousSNP(){
    long cnt = 0;
    for (unsigned i=0; i<numSnps; ++i) {
        SnpInfo *snp = snpInfoVec[i];
        if ((snp->a1 == "T" && snp->a2 == "A") ||
            (snp->a1 == "A" && snp->a2 == "T") ||
            (snp->a1 == "G" && snp->a2 == "C") ||
            (snp->a1 == "C" && snp->a2 == "G")) {
            snp->included = false;
            ++cnt;
        }
    }
    LOGGER << cnt << " SNPs with ambiguous nucleotides, i.e. A/T or G/C, are excluded." << endl;
}

void Data::excludeSNPwithMaf(const double mafmin, const double mafmax){
    unsigned cntmin=0, cntmax=0;
    for (unsigned i=0; i<numSnps; ++i) {
        SnpInfo *snp = snpInfoVec[i];
        double maf = snp->af < 0.5 ? snp->af : 1.0 - snp->af;
        if (mafmin && maf < mafmin) {
            snp->included = false;
            ++cntmin;
        }
        if (mafmax && maf > mafmax) {
            snp->included = false;
            ++cntmax;
        }
    }
    if (mafmin) LOGGER << cntmin << " SNPs with MAF below " << mafmin << " are excluded." << endl;
    if (mafmax) LOGGER << cntmax << " SNPs with MAF above " << mafmax << " are excluded." << endl;
}

void Data::excludeRegion(const string &excludeRegionFile){
    ifstream in(excludeRegionFile.c_str());
    if (!in) LOGGER.e(0,"can not open the file [" + excludeRegionFile + "] to read.");

    map<int, vector<pair<long, long> > > regions;
    int chr;
    long regStart, regEnd;
    while (in >> chr >> regStart >> regEnd) {
        regions[chr].push_back(pair<long, long>(regStart, regEnd));
    }
    
    unsigned cnt = 0;
    map<int, vector<pair<long, long> > >::iterator it, end = regions.end();
    for (unsigned i=0; i<numSnps; ++i) {
        SnpInfo *snp = snpInfoVec[i];
        it = regions.find(snp->chrom);
        if (it != end) {
            long size = it->second.size();
            for (unsigned j=0; j<size; ++j) {
                regStart = it->second[j].first;
                regEnd   = it->second[j].second;
                if (snp->physPos > regStart && snp->physPos < regEnd) {
                    snp->included = false;
                    ++cnt;
                }
            }
        }
    }
    LOGGER << cnt << " SNPs are excluded due to --exclude-region." << endl;
}

void Data::reindexSnp(vector<SnpInfo*> snpInfoVec){
    SnpInfo *snp;
    for (unsigned i=0, idx=0; i<snpInfoVec.size(); ++i) {
        snp = snpInfoVec[i];
        if (snp->included) {
            snp->index = idx++;
        } else {
            snp->index = -9;
        }
    }
}

void Data::includeMatchedSnp(){
    reindexSnp(snpInfoVec);  // reindex for MPI purpose in terms of full snplist
    fullSnpFlag.resize(numSnps);
    for (int i=0; i<numSnps; ++i) fullSnpFlag[i] = snpInfoVec[i]->included; // for output purpose
    incdSnpInfoVec = makeIncdSnpInfoVec(snpInfoVec);
    numIncdSnps = (unsigned) incdSnpInfoVec.size();
    reindexSnp(incdSnpInfoVec);
    snp2pq.resize(numIncdSnps);
    
    map<int, vector<SnpInfo*> > chrmap;
    map<int, vector<SnpInfo*> >::iterator it;
    gwasSnpIdx2snpIDMap.clear();
    for (unsigned i=0; i<numIncdSnps; ++i) {
        SnpInfo *snp = incdSnpInfoVec[i];
        if (chrmap.find(snp->chrom) == chrmap.end()) {
            chrmap[snp->chrom] = *new vector<SnpInfo*>;
        }
        chrmap[snp->chrom].push_back(snp);

        gwasSnpIdx2snpIDMap.insert(pair<int,string>(i,snp->rsID));
    }
    numChroms = (unsigned) chrmap.size();
    chromInfoVec.clear();
    for (it=chrmap.begin(); it!=chrmap.end(); ++it) {
        int id = it->first;
        vector<SnpInfo*> &vec = it->second;
        ChromInfo *chr = new ChromInfo(id, (unsigned)vec.size(), vec[0]->index, vec.back()->index,vec[0]->physPos,vec.back()->physPos);
        chr->startSnpID = vec[0]->rsID;
        chr->endSnpID = vec.back()->rsID; 
        chromInfoVec.push_back(chr);
        chromInfoMap.insert(pair<int, ChromInfo*>(id, chr));
        //LOGGER << "size chrom " << id << ": " << vec.back()->physPos - vec[0]->physPos << endl;
    }

    LOGGER << numIncdSnps << " SNPs on " << numChroms << " chromosomes (ID: ";
    for(unsigned k = 0; k < chromInfoVec.size(); k ++){
        LOGGER << chromInfoVec[k]->id << " ";
    }
    LOGGER << ") are included." << endl;
    
//    if (numAnnos) setAnnoInfoVec();

//    if (numAnnos) {
//        if (myMPI::rank==0) LOGGER << "\nAnnotation info:" << endl;
//        numSnpAnnoVec.resize(numAnnos);
//        for (unsigned i=0; i<numAnnos; ++i) {
//            AnnoInfo *anno = annoInfoVec[i];
//            numSnpAnnoVec[i] = anno->size;
//            anno->getSnpInfo();
//            anno->fraction = double(anno->size)/double(numIncdSnps);
//            anno->print();
//        }
//        if (myMPI::rank==0) LOGGER << endl;
//    }
}

vector<SnpInfo*> Data::makeIncdSnpInfoVec(const vector<SnpInfo*> &snpInfoVec){
    vector<SnpInfo*> includedSnps;
    includedSnps.reserve(numSnps);
    snpEffectNames.reserve(numSnps);
    SnpInfo *snp = NULL;
    for (unsigned i=0; i<numSnps; ++i) {
        snp = snpInfoVec[i];
        if(snp->included) {
            //snp->index = j++;  // reindex snps
            includedSnps.push_back(snp);
            snpEffectNames.push_back(snp->rsID);
        }
    }
    if(includedSnps.size() == 0){
        LOGGER.e(0,"Zero gwas snps are kept.");
    }
    return includedSnps;
}

vector<IndInfo*> Data::makeKeptIndInfoVec(const vector<IndInfo*> &indInfoVec){
    vector<IndInfo*> keptInds;
    keptInds.reserve(numInds);
    gwasIndNames.clear(); // for debug
    IndInfo *ind = NULL;
    for (unsigned i=0, j=0; i<numInds; ++i) {
        ind = indInfoVec[i];
        if(ind->kept) {
            ind->index = j++;  // reindex inds
            keptInds.push_back(ind);
            gwasIndNames.push_back(ind->famID ); // for debug
        }
    }
    if(keptInds.size() == 0){
        LOGGER.e(0,"Zero individuals in the GWAS data are kept.");
    }
    return keptInds;
}

void Data::computeAlleleFreq(const MatrixXd &Z, vector<SnpInfo*> &incdSnpInfoVec, VectorXd &snp2pq){
    LOGGER << "Computing allele frequencies ..." << endl;
    snp2pq.resize(numIncdSnps);
    SnpInfo *snp = NULL;
    for (unsigned i=0; i<numIncdSnps; ++i) {
        snp = incdSnpInfoVec[i];
        snp->af = 0.5f*Z.col(i).mean();
        snp2pq[i] = snp->twopq = 2.0f*snp->af*(1.0-snp->af);
    }
}

void Data::getWindowInfo(const vector<SnpInfo*> &incdSnpInfoVec, const unsigned windowWidth, VectorXi &windStart, VectorXi &windSize){
    LOGGER << "Creating windows (window width: " + to_string(static_cast<long long>(windowWidth/1e6)) + "Mb) ..." << endl;
    int i=0, j=0;
    windStart.setZero(numIncdSnps);
    windSize.setZero(numIncdSnps);
    SnpInfo *snpi, *snpj;
    for (i=0; i<numIncdSnps; ++i) {
        snpi = incdSnpInfoVec[i];
        snpi->resetWindow();
        for (j=i; j>=0; --j) {
            snpj = incdSnpInfoVec[j];
            if (snpi->isProximal(*snpj, windowWidth/2)) {
                snpi->windStart = snpj->index;
                snpi->windSize++;
            } else break;
        }
        for (j=i+1; j<numIncdSnps; ++j) {
            snpj = incdSnpInfoVec[j];
            if (snpi->isProximal(*snpj, windowWidth/2)) {
                snpi->windSize++;
            } else break;
        }
        if(!(i%10000))
            LOGGER << "SNP " << i << " Window Size " << snpi->windSize << endl;
        if (!snpi->windSize) {
            LOGGER.e(0," SNP " + snpi->rsID + " has zero SNPs in its window!");
        }
        windStart[i] = snpi->windStart;
        windSize [i] = snpi->windSize;
        snpi->windEnd = snpi->windStart + snpi->windSize - 1;
    }
}

void Data::getNonoverlapWindowInfo(const unsigned windowWidth){
    if (!windowWidth) LOGGER.e(0," Did you forget to set window width by --wind [Mb]?");
    unsigned window = 0;
    unsigned currChr = incdSnpInfoVec[0]->chrom;
    unsigned long startPos = incdSnpInfoVec[0]->physPos;
    vector<int> windStartVec = {0};
    SnpInfo *snp;
    for (unsigned i=0; i<numIncdSnps; ++i) {
        snp = incdSnpInfoVec[i];
        if (snp->physPos - startPos > windowWidth || snp->chrom > currChr) {
            currChr = snp->chrom;
            startPos = snp->physPos;
            windStartVec.push_back(i);
            ++window;
        }
        snp->window = window;
    }
    numWindows = windStartVec.size();
    makeWindows = true;
    windStart = VectorXi::Map(&windStartVec[0], numWindows);
    windSize.setZero(numWindows);

    for (unsigned i=0; i<numWindows; ++i) {
        if (i != numWindows-1)
            windSize[i] = windStart[i+1] - windStart[i];
        else
            windSize[i] = numIncdSnps - windStart[i];
    }
    LOGGER << "Created " << numWindows << " non-overlapping " << windowWidth/1e3 << "kb windows with average size of " << windSize.sum()/double(numWindows) << " SNPs." << endl;
}

void Data::outputSnpResults(const VectorXd &posteriorMean, const VectorXd &posteriorSqrMean, const VectorXd &pip, const bool noscale, const string &filename) const {
    ofstream out(filename.c_str());
    out << boost::format("%6s %20s %6s %12s %6s %6s %12s %12s %12s %8s")
    % "Id"
    % "Name"
    % "Chrom"
    % "Position"
    % "A1"
    % "A2"
    % "A1Frq"
    % "A1Effect"
    % "SE"
    % "PIP";
    if (makeWindows) out << boost::format("%8s") % "Window";
    out << endl;
    for (unsigned i=0, idx=0; i<numSnps; ++i) {
        SnpInfo *snp = snpInfoVec[i];
        if(!fullSnpFlag[i]) continue;
        //        if(snp->isQTL) continue;)
        double sqrt2pq = sqrt(2.0*snp->af*(1.0-snp->af));
        double effect = (snp->flipped ? -posteriorMean[idx] : posteriorMean[idx]);
        double se = sqrt(posteriorSqrMean[idx]-posteriorMean[idx]*posteriorMean[idx]);
        out << boost::format("%6s %20s %6s %12s %6s %6s %12.4f %12.4e %12.4e %8.4e")
        % (idx+1)
        % snp->rsID
        % snp->chrom
        % snp->physPos
        % (snp->flipped ? snp->a2 : snp->a1)
        % (snp->flipped ? snp->a1 : snp->a2)
        % (snp->flipped ? 1.0-snp->af : snp->af)
        % (noscale ? effect : effect/sqrt2pq)
        % (noscale ? se : se/sqrt2pq)
        % pip[idx];
        if (makeWindows) out << boost::format("%8s") % snp->window;
        out << endl;
        ++idx;
    }
    out.close();
}

void Data::outputSnpResults(const VectorXd &posteriorMean, const VectorXd &posteriorSqrMean, const VectorXd &lastSample, const VectorXd &pip, const bool noscale, const string &filename) const {
    ofstream out(filename.c_str());
    out << boost::format("%6s %20s %6s %12s %6s %6s %12s %15s %15s %15s %15s")
    % "Id"
    % "Name"
    % "Chrom"
    % "Position"
    % "A1"
    % "A2"
    % "A1Frq"
    % "A1Effect"
    % "SE"
    % "PIP"
    % "LastSampleEff";
    if (makeWindows) out << boost::format("%8s") % "Window";
    out << endl;
    for (unsigned i=0, idx=0; i<numSnps; ++i) {
        SnpInfo *snp = snpInfoVec[i];
        if(!fullSnpFlag[i]) continue;
        //        if(snp->isQTL) continue;)
        // double sqrt2pq = sqrt(2.0*snp->af*(1.0-snp->af));
        double sqrtScaleFactor = (snp->scaleFactor);
        double effect = (snp->flipped ? -posteriorMean[idx] : posteriorMean[idx]);
        double lastBeta = (snp->flipped ? -lastSample[idx] : lastSample[idx]);
        double se = sqrt(posteriorSqrMean[idx]-posteriorMean[idx]*posteriorMean[idx]);
        out << boost::format("%6s %20s %6s %12s %6s %6s %12.6f %15.6e %15.6e %15.6e %15.6e")
        % (idx+1)
        % snp->rsID
        % snp->chrom
        % snp->physPos
        % (snp->flipped ? snp->a2 : snp->a1)
        % (snp->flipped ? snp->a1 : snp->a2)
        % (snp->flipped ? 1.0-snp->af : snp->af)
        % (noscale ? effect : effect/sqrtScaleFactor)
        // % ( effect )
        % (noscale ? se : se/sqrtScaleFactor)
        % pip[idx]
        % (noscale ? lastBeta : lastBeta/sqrtScaleFactor);
        if (makeWindows) out << boost::format("%8s") % snp->window;
        out << endl;
        ++idx;
    }
    out.close();
}

void Data::outputGeneEffectResults(const VectorXd &posteriorMean, const VectorXd &posteriorSqrMean,const VectorXd &lastSample, const VectorXd &pip,const string &mcmcType,const string &filename) const {
    ofstream out(filename.c_str());
    out << boost::format("%6s %20s %6s %12s %15s %15s %15s %15s")
    % "Id"
    % "Name"
    % "Chrom"
    % "Position"
    % "Effect"
    % "SE"
    % "PIP";
    if (makeWindows) out << boost::format("%8s") % "Window";
    out << endl;
    for (unsigned i=0, idx=0; i<numKeptGenes; ++i) {
        GeneInfo *gene = keptGeneInfoVec[i];
        // double effect = (snp->flipped ? -posteriorMean[idx] : posteriorMean[idx]);
        // double se = sqrt(posteriorSqrMean[idx]-posteriorMean[idx]*posteriorMean[idx]);
        // out << boost::format("%6s %20s %6s %12s %15.6e %15.6e %15.6e %15.6e")
        // % (idx+1)
        // % gene->ensemblID
        // % gene->chrom
        // % gene->midPhyPos
        // % effect 
        // % se 
        // % 0 // pip[idx]
        // % lastEff;
        // out << endl;
        // ++idx;
    }
    out.close();
}



void Data::inputSnpResults(const string &snpResFile){
    ifstream in(snpResFile.c_str());
    if (!in) LOGGER.e(0,"can not open the SNP result file [" + snpResFile + "] to read.");
    LOGGER << "Reading SNP results from [" + snpResFile + "]." << endl;
    
    SnpInfo *snp;
    map<string, SnpInfo*>::iterator it;
    string name;
    int id, chrom, pos, window;
    double freq, effect, se, pip;
    unsigned line=0, match=0;
    string header;
    getline(in, header);
    while (in >> id >> name >> chrom >> pos >> freq >> effect >> se >> pip >> window) {
        ++line;
        it = snpInfoMap.find(name);
        if (it == snpInfoMap.end()) continue;
        snp = it->second;
        if (snp->included) {
            snp->effect = effect;
            ++match;
        }
    }
    in.close();
    
    LOGGER << match << " matched SNPs in the SNP result file (in total " << line << " SNPs)." << endl;
}

void Data::inputSnpInfoAndResults(const string &snpResFile, const string &bayesType){
    ifstream in(snpResFile.c_str());
    if (!in) LOGGER.e(0,"can not open the SNP result file [" + snpResFile + "] to read.");
    LOGGER << "Reading SNP info and results from [" + snpResFile + "]." << endl;
    
    SnpInfo *snp;
    map<string, SnpInfo*>::iterator it;
    string name;
    int id, chrom, pos, window;
    double freq, effect, se, pip;
    unsigned line=0, match=0;
    string header;
    getline(in, header);
    for (unsigned i=0; i<numSnps; ++i) {
        snp = snpInfoVec[i];
        snp->included = false;
    }
    if (bayesType == "SMix") {
        double piS;
        while (in >> id >> name >> chrom >> pos >> freq >> effect >> se >> pip >> piS) {
            ++line;
            it = snpInfoMap.find(name);
            if (it == snpInfoMap.end()) {
                LOGGER.e(0," SNP " + name + " is not in the LD matrix!");
            }
            snp = it->second;
            snp->included = true;
            snp->gwas_af = freq;
            snp->effect = effect;
        }
    }
    else {
        while (in >> id >> name >> chrom >> pos >> freq >> effect >> se >> pip >> window) {
            ++line;
            //        SnpInfo *snp = new SnpInfo(id-1, name, "NA", "NA", chrom, 0, pos);
            it = snpInfoMap.find(name);
            if (it == snpInfoMap.end()) {
                LOGGER.e(0," SNP " + name + " is not in the LD matrix!");
            }
            snp = it->second;
            snp->included = true;
            snp->gwas_af = freq;
            snp->effect = effect;
            //        snpInfoVec.push_back(snp);
            //        snpInfoMap.insert(pair<string, SnpInfo*>(name, snp));
            //        chromosomes.insert(snp->chrom);
        }
    }
    in.close();
    
//    numSnps = (unsigned) snpInfoVec.size();
    LOGGER << line << " SNPs in the SNP result file." << endl;
}

void Data::summarizeSnpResults(const SpMat &snpEffects, const string &filename) const {
    LOGGER << "SNP results to be summarized in " << filename << endl;
    unsigned nrow = snpEffects.rows();
    VectorXd effectSum(numIncdSnps), effectMean(numIncdSnps);
    VectorXd pipSum(numIncdSnps), pip(numIncdSnps);  // posterior inclusion probability
    for (unsigned i=0; i<numIncdSnps; ++i) {
        effectSum[i] = snpEffects.col(i).sum();
        pipSum[i] = (VectorXd(snpEffects.col(i)).array()!=0).count();
    }
    effectMean = effectSum/(double)nrow;
    pip = pipSum/(double)nrow;
    
    ofstream out(filename.c_str());
    out << boost::format("%6s %20s %6s %12s %6s %6s %12s %12s %8s %8s\n")
    % "Id"
    % "Name"
    % "Chrom"
    % "Position"
    % "A1"
    % "A2"
    % "A1Frq"
    % "A1Effect"
    % "PIP"
    % "Window";
    for (unsigned i=0, idx=0; i<numSnps; ++i) {
        SnpInfo *snp = snpInfoVec[i];
        if(!fullSnpFlag[i]) continue;
        out << boost::format("%6s %20s %6s %12s %6s %6s %12.6f %12.6f %8.3f %8s\n")
        % (idx+1)
        % snp->rsID
        % snp->chrom
        % snp->physPos
        % (snp->flipped ? snp->a2 : snp->a1)
        % (snp->flipped ? snp->a1 : snp->a2)
        % (snp->flipped ? 1.0-snp->af : snp->af)
        % (snp->flipped ? -effectMean[idx] : effectMean[idx])
        % pip[idx]
        % snp->window;
        ++idx;
    }
    out.close();
}


void Data::outputFixedEffects(const MatrixXd &fixedEffects, const string &filename) const {
    ofstream out(filename.c_str());
    long nrow = fixedEffects.rows();
    VectorXd mean = fixedEffects.colwise().mean();
    VectorXd sd = (fixedEffects.rowwise() - mean.transpose()).colwise().squaredNorm().cwiseSqrt()/sqrt(nrow);
    for (unsigned i=0; i<numFixedEffects; ++i) {
        out << boost::format("%20s %12.6f %12.6f\n") % fixedEffectNames[i] %mean[i] %sd[i];
    }
    out.close();
}

void Data::outputRandomEffects(const MatrixXd &randomEffects, const string &filename) const {
    ofstream out(filename.c_str());
    long nrow = randomEffects.rows();
    VectorXd mean = randomEffects.colwise().mean();
    VectorXd sd = (randomEffects.rowwise() - mean.transpose()).colwise().squaredNorm().cwiseSqrt()/sqrt(nrow);
    for (unsigned i=0; i<numRandomEffects; ++i) {
        out << boost::format("%20s %12.6f %12.6f\n") % randomEffectNames[i] %mean[i] %sd[i];
    }
    out.close();
}

void Data::outputWindowResults(const VectorXd &posteriorMean, const string &filename) const {
    ofstream out(filename.c_str());
    out << boost::format("%6s %8s\n") %"Id" %"PIP";
    for (unsigned i=0; i<posteriorMean.size(); ++i) {
        out << boost::format("%6s %8.3f\n")
        % (i+1)
        % posteriorMean[i];
    }
    out.close();
}

void Data::readGwasSummaryFile(const string &gwasFile, const double afDiff, const double mafmin, const double mafmax, const double pValueThreshold, const bool imputeN, const bool removeOutlierN){
    ifstream in(gwasFile.c_str());
    if (!in) LOGGER.e(0,"can not open the GWAS summary data file [" + gwasFile + "] to read.");
    LOGGER  << "............................" << endl;
    LOGGER << "Reading GWAS summary data..." << endl;
    LOGGER  << "............................" << endl;
    LOGGER << "Reading GWAS summary data from [" + gwasFile + "]." << endl;
    
    string header;
    getline(in, header);

    SnpInfo *snp;
    map<string, SnpInfo*>::iterator it;
    string id, allele1, allele2, freq, b, se, pval, n;
    unsigned line=0, match=0;
    unsigned numInconAllele=0, numInconAf=0, numFixed=0, numMafMin=0, numMafMax=0, numOutlierN=0;
    unsigned numPvalPruned=0;
    unsigned numFlip=0;
    bool inconAllele, inconAf, fixed, ismafmin, ismafmax, isPvalPruned;
    double gwas_af;
    while (in >> id >> allele1 >> allele2 >> freq >> b >> se >> pval >> n) {
        it = snpInfoMap.find(id);
        if (it == snpInfoMap.end()) {
            //LOGGER << "cannot find SNP " << id << endl;
            continue;
        }
        snp = it->second;
        if (!snp->included) {
            //LOGGER << "exclude SNP " << id << endl;
            continue;
        }

        inconAllele = inconAf = fixed = ismafmin = ismafmax = isPvalPruned = false;
        if (allele1 == snp->a1 && allele2 == snp->a2) {
            gwas_af = std::stod(freq.c_str());
            snp->gwas_b  = std::stod(b.c_str());
            snp->gwas_af = gwas_af != -999 ? gwas_af : snp->af;  // set -1 in gwas summary file if allele frequencies are not available
            snp->gwas_se = std::stod(se.c_str());
            snp->gwas_n  = std::stod(n.c_str());
            cpp_dec_float_100 highPrecisionValue = cpp_dec_float_100(pval);
            snp->gwas_pvalue = static_cast<double>(highPrecisionValue);
            // snp->gwas_pvalue = std::stod(pval.c_str());
        } else if (allele1 == snp->a2 && allele2 == snp->a1) {
            gwas_af = std::stod(freq.c_str());
            snp->gwas_b  = -std::stod(b.c_str());
            snp->gwas_af = gwas_af != -999 ? 1.0 - gwas_af : 1.0 - snp->af;
            snp->gwas_se = std::stod(se.c_str());
            snp->gwas_n  = std::stod(n.c_str());
            // snp->gwas_pvalue = std::stod(pval.c_str());
            cpp_dec_float_100 highPrecisionValue = cpp_dec_float_100(pval);
            snp->gwas_pvalue = static_cast<double>(highPrecisionValue);
            snp->flipped = true;
            
//            snp->included = false;
            //LOGGER << snp->index << " " << snp->rsID << " " << snp->gwas_af << " " << snp->gwas_b << endl;
            ++numFlip;
        } else {
//            LOGGER << "WARNING: SNP " + id + " has inconsistent allele coding in between the reference and GWAS samples." << endl;
            inconAllele = true;
            ++numInconAllele;
        }
        if (!inconAllele) {
            if (abs(snp->af - snp->gwas_af) > afDiff) {
                inconAf = true;
                ++numInconAf;
            } else if (snp->gwas_af==0 || snp->gwas_af==1) {
                fixed = true;
                ++numFixed;
            } else if (mafmin || mafmax) {
                double maf_ref = snp->af < 0.5 ? snp->af : 1.0 - snp->af;
                double maf_gwas = snp->gwas_af < 0.5 ? snp->gwas_af : 1.0 - snp->gwas_af;
                if (mafmin && (maf_ref < mafmin || maf_gwas < mafmin)) {
                    ismafmin = true;
                    ++numMafMin;
                }
                if (mafmax && (maf_ref > mafmax || maf_gwas > mafmax)) {
                    ismafmax = true;
                    ++numMafMax;
                }
            }
            if (snp->gwas_pvalue > pValueThreshold) {
                isPvalPruned = true;
                ++numPvalPruned;
            }
        }
        if (inconAllele || inconAf || fixed || ismafmin || ismafmax || isPvalPruned) {
            snp->included = false;
            LOGGER << snp->index << " " << snp->rsID << endl;
        } else ++match;
        
        ++line;
    }  // end of while loop
    in.close();
        
    if (removeOutlierN){
        unsigned size = 0;
        for (unsigned i=0; i<numSnps; ++i) {
            snp = snpInfoVec[i];
            if (snp->included && snp->gwas_n != -999) ++size;
        }
        ArrayXd perSnpN(size);
        vector<SnpInfo*> snpvec(size);
        for (unsigned i=0, j=0; i<numSnps; ++i) {
            snp = snpInfoVec[i];
            if (snp->included && snp->gwas_n != -999) {
                perSnpN[j] = snp->gwas_n;
                snpvec[j] = snp;
                ++j;
            }
        }
        
        double n_med = Gadget::findMedian(perSnpN);
        double sd = sqrt(Gadget::calcVariance(perSnpN));
        for (unsigned i=0; i<size; ++i) {
            snp = snpvec[i];
            if (perSnpN[i] < n_med - 3*sd || perSnpN[i] > n_med + 3*sd) {
                snp->included = false;
                ++numOutlierN;
            }
        }
        
        match -= numOutlierN;
    }
    
    numIncdSnps = 0;
    for (unsigned i=0; i<numSnps; ++i) {
        snp = snpInfoVec[i];
        if (!snp->included) continue;
        if (snp->gwas_b == -999) {
            snp->included = false;
        } else {
            ++numIncdSnps;
        }
    }
    numGWASFlip = numFlip;
    if (numFlip) LOGGER << "flipped " << numFlip << " SNPs according to the minor allele in the reference and GWAS samples." << endl;
    if (numInconAllele) LOGGER << "removed " << numInconAllele << " SNPs with inconsistent allele coding in between the reference and GWAS samples." << endl;
    if (numInconAf) LOGGER << "removed " << numInconAf << " SNPs with differences in allele frequency between the reference and GWAS samples > " << afDiff << "." << endl;
    if (numFixed) LOGGER << "removed " << numFixed << " fixed SNPs in the GWAS samples." << endl;
    if (mafmin) LOGGER << "removed " << numMafMin << " SNPs with MAF below " << mafmin << " in either reference and GWAS samples." << endl;
    if (mafmax) LOGGER << "removed " << numMafMax << " SNPs with MAF above " << mafmax << " in either reference and GWAS samples." << endl;
    if (pValueThreshold < 1.0) LOGGER << "removed " << numPvalPruned << " SNPs with GWAS P value greater than " << pValueThreshold << "." << endl;
    if (numOutlierN) LOGGER << "removed " << numOutlierN << " SNPs with per-SNP sample size beyond 3 SD around the median value." << endl;
    if(match == 0){
        LOGGER.e(0,to_string(match) + " matched SNPs in the GWAS summary data (in total "+ to_string(line) + " SNPs).");    
    }
    LOGGER << match << " matched SNPs in the GWAS summary data (in total " << line << " SNPs)." << endl;

    
    if (imputeN) imputePerSnpSampleSize(snpInfoVec, numIncdSnps, 0);
}

void Data::imputePerSnpSampleSize(vector<SnpInfo*> &snpInfoVec, unsigned &numIncdSnps, double sd) {
    // use input allele frequencies, b_hat and se to impute per-snp N
    // then filter SNPs with N > 3 sd apart from the median value
    ArrayXd n(numIncdSnps);
    ArrayXd p(numIncdSnps);
    ArrayXd bsq(numIncdSnps);
    ArrayXd var(numIncdSnps);
    ArrayXd tpq(numIncdSnps);
    ArrayXd ypy(numIncdSnps);
    SnpInfo *snp;
    unsigned j = 0;
    for (unsigned i=0; i<numSnps; ++i) {
        snp = snpInfoVec[i];
        if (!snp->included) continue;
        n[j] = snp->gwas_n;
        p[j] = snp->gwas_af;
        bsq[j] = snp->gwas_b*snp->gwas_b;
        var[j] = snp->gwas_se*snp->gwas_se;
        ++j;
    }
    tpq = 2.0*p*(1.0-p);
    ypy = tpq*n.square()*var + tpq*n*bsq;
    double ypy_med = Gadget::findMedian(ypy);
    // Given ypy and n compute 2pq
    //tpq = ypy / (var*n.square() + bsq*n);
    // Given ypy_med and 2pq compute n
    double n_med = Gadget::findMedian(n);
    double vary = ypy_med / n_med;
    n = (vary - tpq*bsq) / (tpq*var);
    // compute sd of n
    double sdOld = sd;
    sd = sqrt(Gadget::calcVariance(n));
    double p_new;
    j = 0;
    numIncdSnps = 0;
    for (unsigned i=0; i<numSnps; ++i) {
        snp = snpInfoVec[i];
        if (!snp->included) continue;
        if (n[j] < n_med - 3*sd || n[j] > n_med + 3*sd) {
            snp->included = false;
        }
        else {
            snp->gwas_n = n[j];
            //p_new = 0.5 - 0.5*sqrt(1.0-2.0*tpq[j]);
            //snp->gwas_af = p[j] < 0.5 ? p_new : 1.0-p_new;
            ++numIncdSnps;
        }
        ++j;
    }
//    LOGGER << n.mean() << " " << ypy_med << " " << sdOld << " " << sd << " " << numIncdSnps << " " << n.head(10).transpose() << endl;
    if (abs(sd-sdOld) > 0.01) {
        imputePerSnpSampleSize(snpInfoVec, numIncdSnps, sd);
    } else {
        LOGGER << numIncdSnps << " SNPs with per-SNP sample size within 3 sd around the median value of " << n_med << endl;
        string outfile = title + ".imputedPerSnpN";
        ofstream out(outfile.c_str());
        out << boost::format("%15s %12s\n")
        % "ID" % "Imputed_N";
        for (unsigned i=0; i<numSnps; ++i) {
            snp = snpInfoVec[i];
            if (!snp->included) continue;
            out << boost::format("%15s %12s\n")
            % snp->rsID
            % snp->gwas_n;
        }
        out.close();
        return;
    }
}

/*
 * Divide the big matrix to small part (From Zhili)
 * @Param totalPart: number of parts would like to divide 
 * @Param curPart:  current part would like to run. 1 based
 * @Param nSNPs:  total number of SNPs
 * @Return string: SNP idx start-end, 1 based
*/
string makeSNPRangeString(int totalPart, int curPart, int nSNPs){
    int nPartSNP = (nSNPs + totalPart - 1) / totalPart;
    int start = nPartSNP * (curPart - 1) + 1;
    int end = nPartSNP * curPart;
    if(end > nSNPs) end = nSNPs;
    return(to_string(start) + "-" + to_string(end));
}

string Data::partLDMatrix(const string &partParam, const string &outfilename, const string &LDmatType){
    Gadget::Tokenizer token;
    token.getTokens(partParam, ",");
    string snpRange = "";
    if(token.size() == 2){
        int nTotalPart = atoi(token[0].c_str());
        int nCurPart = atoi(token[1].c_str());
        LOGGER << "Dividing LD matrix by " << nTotalPart << " parts, current running part " << nCurPart << std::endl;
        if(nTotalPart < nCurPart || nCurPart <= 0){
            LOGGER.e(0,"--part usage total,curentPart"); 
        }

        snpRange = makeSNPRangeString(nTotalPart, nCurPart, numIncdSnps);
        LOGGER << "  SNP range " << snpRange << endl;

        if(nCurPart == 1){
            ofstream mldmfile((outfilename + ".mldm").c_str());
            for(int i = 1; i <= nTotalPart; i++){
                string snpRange1 = makeSNPRangeString(nTotalPart, i, numIncdSnps);
                string outfilename2 = outfilename + ".snp" + snpRange1 + ".ldm." + LDmatType;
                mldmfile << outfilename2 << endl;
            }
            mldmfile.close();
            LOGGER << " run gctb --mldm " << outfilename << ".mldm --make-full-ldm --out ...  to combine the parted LDM matrix" << std::endl;  
        }
    }
    return snpRange;
}

void Data::makeLDmatrix(const string &bedFile, const string &LDmatType, const double chisqThreshold, const double LDthreshold, const unsigned windowWidth, const string &snpRange, const string &filename, const bool writeLdmTxt){
    
    Gadget::Tokenizer token;
    token.getTokens(snpRange, "-");
    
    unsigned start = 0;
    unsigned end = numIncdSnps;
    
    if (token.size()) {
        start = atoi(token[0].c_str()) - 1;
        end = atoi(token[1].c_str());
        if (end > numIncdSnps) end = numIncdSnps;
    }
    
    unsigned numSnpInRange = end - start;
    
    if (snpRange.empty())
        LOGGER << "Building " + LDmatType + " LD matrix for all SNPs ..." << endl;
    else
        LOGGER << "Building " + LDmatType + " LD matrix for SNPs " << snpRange << " ..." << endl;
    
    if (numIncdSnps == 0) LOGGER.e(0,"No SNP is retained for analysis.");
    if (numKeptInds == 0) LOGGER.e(0,"No individual is retained for analysis.");
    if (start >= numIncdSnps) LOGGER.e(0,"Specified a SNP range of " + snpRange + " but " + to_string(static_cast<long long>(numIncdSnps)) + " SNPs are included.");
    
    Gadget::Timer timer;
    timer.setTime();
    
    unsigned firstWindStart = 0;
    unsigned lastWindEnd = 0;
    if (windowWidth) {
        getWindowInfo(incdSnpInfoVec, windowWidth, windStart, windSize);
    }
    
    // first read in the genotypes of SNPs in the given range
    
    const int bedToGeno[4] = {2, -9, 1, 0};
    unsigned size = (numInds+3)>>2;
    
    MatrixXd ZP(numSnpInRange, numKeptInds);  // SNP x Ind
    D.setZero(numSnpInRange);

    if (numKeptInds < 2) LOGGER.e(0," Cannot calculate LD matrix with number of individuals < 2.");
    
    FILE *in1 = fopen(bedFile.c_str(), "rb");
    if (!in1) LOGGER.e(0,"can not open the file [" + bedFile + "] to read.");
    LOGGER << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
    char header[3];
    fread(header, sizeof(header), 1, in1);
    if (!in1 || header[0] != 0x6c || header[1] != 0x1b || header[2] != 0x01) {
        cerr << "Error: Incorrect first three bytes of bed file: " << bedFile << endl;
        exit(1);
    }

    
    IndInfo *indi = NULL;
    SnpInfo *snpj = NULL;
    SnpInfo *snpk = NULL;
    
    int genoValue;
    unsigned i, j, k;
    unsigned incj, inck; // index of included SNP
    unsigned long long skipj = 0;
    unsigned nmiss;
    double mean;
    
    set<int> chromInRange;
    
    for (j = 0, incj = 0; j < numSnps; j++) {
        snpj = snpInfoVec[j];
        
        if (snpj->index < start || !snpj->included) {
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
        mean /= double(snpj->sampleSize);
        if (nmiss) {
            for (i=0; i<numKeptInds; ++i) {
                if (ZP(incj, i) == -9) ZP(incj, i) = mean;
            }
        }
        
        // compute allele frequency
        snpj->af = 0.5f*mean;
        snp2pq[incj] = snpj->twopq = 2.0f*snpj->af*(1.0-snpj->af);
        
        if (snp2pq[incj]==0) LOGGER.e(0,"" + snpj->rsID + " is a fixed SNP (MAF=0)!");
        
        // standardize genotypes
        //D[incj] = snp2pq[incj]*snpj->sampleSize;
        D[incj] = Gadget::calcVariance(ZP.row(incj))*numKeptInds;
        
        ZP.row(incj) = (ZP.row(incj).array() - ZP.row(incj).mean())/sqrt(D[incj]);
//        ZP.row(incj) = (ZP.row(incj).array() - mean)/sqrt(D[incj]);
        
        if (windowWidth) {
            if (incj == 0) firstWindStart = snpj->windStart;
            if (incj == numSnpInRange-1) lastWindEnd = snpj->windStart + snpj->windSize;
        }
        
        if (++incj == numSnpInRange) break;
    }
    
    fclose(in1);
    
//    ZP = ZP.colwise() - ZP.rowwise().mean();
//    ZP = ZP.array().colwise() / D.cwiseSqrt().array();
    ZPZdiag = ZP.rowwise().squaredNorm();
    
//    LOGGER << ZP.rowwise().mean() << endl << endl;
//    LOGGER << ZP.block(0, 0, 10, 10) << endl;
    
    // then read in the bed file again to compute Z'Z
    

    MatrixXd denseZPZ;
    denseZPZ.setZero(numSnpInRange, numIncdSnps);
    VectorXd Zk(numKeptInds);
    D.setZero(numIncdSnps);
    
//    double numJackknife = numKeptInds-1;
//    MatrixXd ZPZkCwise;
//    MatrixXd samplVarEmp;
//    samplVarEmp.setZero(numSnpInRange, numIncdSnps);
    
    FILE *in2 = fopen(bedFile.c_str(), "rb");
    fseek(in2, 3, SEEK_SET);
    unsigned long long skipk = 0;
    
    set<int>::iterator setend = chromInRange.end();

    if (numSkeletonSnps) {
        for (k = 0, inck = 0; k < numSnps; k++) {
            snpk = snpInfoVec[k];
            
            if (!snpk->included) {
                skipk += size;
                continue;
            }

            if(!(inck%1000)) LOGGER << " read snp " << inck << "\r" << flush;

            if (chromInRange.find(snpk->chrom) == setend && !snpk->skeleton) {
                skipk += size;
                ++inck;       // ensure the index is correct
                continue;
            }
            
            if (windowWidth) {
                if (inck < firstWindStart) {
                    skipk += size;
                    continue;
                } else if (inck > lastWindEnd) {
                    break;
                }
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
            mean /= double(snpk->sampleSize);
            if (nmiss) {
                for (i=0; i<numKeptInds; ++i) {
                    if (Zk[i] == -9) Zk[i] = mean;
                }
            }
            
            // compute allele frequency
            snpk->af = 0.5f*mean;
            snp2pq[inck] = snpk->twopq = 2.0f*snpk->af*(1.0-snpk->af);
            
            if (snp2pq[inck]==0) LOGGER.e(0,"" + snpk->rsID + " is a fixed SNP (MAF=0)!");
            
            // standardize genotypes
            //D[inck] = snp2pq[inck]*snpk->sampleSize;
            D[inck] = Gadget::calcVariance(Zk.row(inck))*numKeptInds;

            Zk = (Zk.array() - Zk.mean())/sqrt(D[inck]);
//            Zk = (Zk.array() - mean)/sqrt(D[inck]);
            
            denseZPZ.col(inck) = ZP * Zk;

//            LOGGER << " inck " << inck << " snpk " << k << " chr " << snpk->chrom << " " << ZP*Zk << endl;
            
            ++inck;
        }
    }
    else {
        for (k = 0, inck = 0; k < numSnps; k++) {
            snpk = snpInfoVec[k];
            
            if (!snpk->included) {
                skipk += size;
                continue;
            }
            
            if (windowWidth) {
                if (inck < firstWindStart) {
                    skipk += size;
                    continue;
                } else if (inck > lastWindEnd) {
                    break;
                }
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
            mean /= double(snpk->sampleSize);
            if (nmiss) {
                for (i=0; i<numKeptInds; ++i) {
                    if (Zk[i] == -9) Zk[i] = mean;
                }
            }
            
            // compute allele frequency
            snpk->af = 0.5f*mean;
            snp2pq[inck] = snpk->twopq = 2.0f*snpk->af*(1.0-snpk->af);
            
            if (snp2pq[inck]==0) LOGGER.e(0,"" + snpk->rsID + " is a fixed SNP (MAF=0)!");
            
            // standardize genotypes
            //D[inck] = snp2pq[inck]*snpk->sampleSize;
            D[inck] = Gadget::calcVariance(Zk)*numKeptInds;

            Zk = (Zk.array() - Zk.mean())/sqrt(D[inck]);
//            Zk = (Zk.array() - mean)/sqrt(D[inck]);
            
            denseZPZ.col(inck) = ZP * Zk;
            
//            // Jackknife estimate of correlation and sampling variance
//            ZPZkCwise = ZP.array().rowwise() * Zk.transpose().array();
//            denseZPZ.col(inck) = ZPZkCwise.rowwise().sum();
//            ZPZkCwise = - (ZPZkCwise.colwise() - denseZPZ.col(inck));
//            ZPZkCwise *= numKeptInds/numJackknife;
//            samplVarEmp.col(inck) = (ZPZkCwise.colwise() - ZPZkCwise.rowwise().mean()).rowwise().squaredNorm() * numJackknife/numKeptInds;
//            // Jackknife end
            
            if(!(inck%1000)) LOGGER << " read snp " << inck << "\r" << flush;
            
            ++inck;
        }
    }

    fclose(in2);
    
    //LOGGER << denseZPZ.block(0, 0, 10, 10) << endl;
    

    // find out per-SNP window position
    
    if (LDmatType == "full") {
        ZPZ.resize(numSnpInRange);
        windStart.setZero(numSnpInRange);
        windSize.setConstant(numSnpInRange, numIncdSnps);
        for (unsigned i=0; i<numSnpInRange; ++i) {
            SnpInfo *snp = incdSnpInfoVec[start+i];
            snp->windStart = 0;
            snp->windSize  = numIncdSnps;
            snp->windEnd   = numIncdSnps-1;
            ZPZ[i] = denseZPZ.row(i);
            snp->ldSamplVar = (1.0 - denseZPZ.row(i).array().square()).square().sum()/numKeptInds;
            snp->ldSum = denseZPZ.row(i).sum();
//            snp->ldSum = samplVarEmp.row(i).sum();
       }
    }
    else if (LDmatType == "band") {
        ZPZ.resize(numSnpInRange);
        if (windowWidth) {  // based on the given window width
            for (unsigned i=0; i<numSnpInRange; ++i) {
                SnpInfo *snp = incdSnpInfoVec[start+i];
                ZPZ[i] = denseZPZ.row(i).segment(snp->windStart, snp->windSize);
                snp->ldSamplVar = (1.0 - ZPZ[i].array().square()).square().sum()/numKeptInds;
                snp->ldSum = ZPZ[i].sum();
            }
        } else {  // based on the given LD threshold
            windStart.setZero(numSnpInRange);
            windSize.setZero(numSnpInRange);
            for (unsigned i=0; i<numSnpInRange; ++i) {
                SnpInfo *snp = incdSnpInfoVec[start+i];
                unsigned windEndi = numIncdSnps;
                for (unsigned j=0; j<numIncdSnps; ++j) {
                    if (abs(denseZPZ(i,j)) > LDthreshold) {
                        windStart[i] = snp->windStart = j;
                        break;
                    }
                }
                for (unsigned j=numIncdSnps; j>0; --j) {
                    if (abs(denseZPZ(i,j-1)) > LDthreshold) {
                        windEndi = j;
                        break;
                    }
                }
                windSize[i] = snp->windSize = windEndi - windStart[i];
                snp->windEnd = windEndi - 1;
                ZPZ[i].resize(windSize[i]);
                VectorXd::Map(&ZPZ[i][0], windSize[i]) = denseZPZ.row(i).segment(windStart[i], windSize[i]);
                snp->ldSamplVar = (1.0 - ZPZ[i].array().square()).square().sum()/numKeptInds;
                snp->ldSum = ZPZ[i].sum();
            }
        }
    }
    else if (LDmatType == "sparse") {
        ZPZsp.resize(numSnpInRange);
        windStart.setZero(numSnpInRange);
        windSize.setZero(numSnpInRange);
        double rsq = 0.0;
        SnpInfo *snpi, *snpj;
        if (numSkeletonSnps) {
            if (LDthreshold) {
                for (unsigned i=0; i<numSnpInRange; ++i) {
                    snpi = incdSnpInfoVec[start+i];
                    snpi->ldSamplVar = 0.0;
                    snpi->ldSum = 0.0;
                    for (unsigned j=0; j<numIncdSnps; ++j) {
                        snpj = incdSnpInfoVec[j];
                        if (snpj->skeleton) {
                            rsq = denseZPZ(i,j)*denseZPZ(i,j);
                            snpi->ldSamplVar += (1.0-rsq)*(1.0-rsq)/numKeptInds;
                            snpi->ldSum += denseZPZ(i,j);
                        }
                        else {
                            if (abs(denseZPZ(i,j)) < LDthreshold) denseZPZ(i,j) = 0;
                            else {
                                rsq = denseZPZ(i,j)*denseZPZ(i,j);
                                snpi->ldSamplVar += (1.0-rsq)*(1.0-rsq)/numKeptInds;
                                snpi->ldSum += denseZPZ(i,j);
                            }
                        }
                    }
                    ZPZsp[i] = denseZPZ.row(i).sparseView();
                    SparseVector <double> ::InnerIterator it(ZPZsp[i]);
                    windStart[i] = snpi->windStart = it.index();
                    windSize[i] = snpi->windSize = ZPZsp[i].nonZeros();
                    for (; it; ++it) snpi->windEnd = it.index();
                }
            } else {
                for (unsigned i=0; i<numSnpInRange; ++i) {
                    snpi = incdSnpInfoVec[start+i];
                    snpi->ldSamplVar = 0.0;
                    snpi->ldSum = 0.0;
//                    unsigned cnt = 0;
                    for (unsigned j=0; j<numIncdSnps; ++j) {
                        snpj = incdSnpInfoVec[j];
                        if (snpj->skeleton) {
                            rsq = denseZPZ(i,j)*denseZPZ(i,j);
                            snpi->ldSamplVar += (1.0-rsq)*(1.0-rsq)/numKeptInds;
                            snpi->ldSum += denseZPZ(i,j);
//                            LOGGER << j << " " << denseZPZ(i,j) << endl;
//                            ++cnt;
                        }
                        else {
                            if (i!=j && denseZPZ(i,j)*denseZPZ(i,j)*snpi->sampleSize < chisqThreshold) denseZPZ(i,j) = 0;
                            else {
                                rsq = denseZPZ(i,j)*denseZPZ(i,j);
                                snpi->ldSamplVar += (1.0-rsq)*(1.0-rsq)/numKeptInds;
                                snpi->ldSum += denseZPZ(i,j);
//                                LOGGER << j << " " << denseZPZ(i,j) << endl;
//                                ++cnt;
                            }
                        }
                    }
//                    LOGGER << "cnt " << cnt << endl;
                    ZPZsp[i] = denseZPZ.row(i).sparseView();
                    SparseVector <double> ::InnerIterator it(ZPZsp[i]);
                    windStart[i] = snpi->windStart = it.index();
                    windSize[i] = snpi->windSize = ZPZsp[i].nonZeros();
                    for (; it; ++it) snpi->windEnd = it.index();
                }
            }
        }
        else {
            if (LDthreshold) {
                for (unsigned i=0; i<numSnpInRange; ++i) {
                    snpi = incdSnpInfoVec[start+i];
                    snpi->ldSamplVar = 0.0;
                    snpi->ldSum = 0.0;
                    for (unsigned j=0; j<numIncdSnps; ++j) {
                        snpj = incdSnpInfoVec[j];
                        if (abs(denseZPZ(i,j)) < LDthreshold) denseZPZ(i,j) = 0;
                        else {
                            rsq = denseZPZ(i,j)*denseZPZ(i,j);
                            snpi->ldSamplVar += (1.0-rsq)*(1.0-rsq)/numKeptInds;
                            snpi->ldSum += denseZPZ(i,j);
                        }
                    }
                    ZPZsp[i] = denseZPZ.row(i).sparseView();
                    SparseVector <double> ::InnerIterator it(ZPZsp[i]);
                    windStart[i] = snpi->windStart = it.index();
                    windSize[i] = snpi->windSize = ZPZsp[i].nonZeros();
                    for (; it; ++it) snpi->windEnd = it.index();
                }
            } else {
                for (unsigned i=0; i<numSnpInRange; ++i) {
                    snpi = incdSnpInfoVec[start+i];
                    snpi->ldSamplVar = 0.0;
                    snpi->ldSum = 0.0;
                    for (unsigned j=0; j<numIncdSnps; ++j) {
                        snpj = incdSnpInfoVec[j];
                        if (i!=j && denseZPZ(i,j)*denseZPZ(i,j)*numKeptInds < chisqThreshold) denseZPZ(i,j) = 0;
                        else {
                            rsq = denseZPZ(i,j)*denseZPZ(i,j);
                            snpi->ldSamplVar += (1.0-rsq)*(1.0-rsq)/numKeptInds;
                            snpi->ldSum += denseZPZ(i,j);
                        }
                    }
                    ZPZsp[i] = denseZPZ.row(i).sparseView();
                    SparseVector <double> ::InnerIterator it(ZPZsp[i]);
                    windStart[i] = snpi->windStart = it.index();
                    windSize[i] = snpi->windSize = ZPZsp[i].nonZeros();
                    for (; it; ++it) snpi->windEnd = it.index();
                }
            }
        }
    }
    
    denseZPZ.resize(0,0);
    
    //    LOGGER << denseZPZ.block(0,0,10,10) << endl;
    //    LOGGER << windStart.transpose() << endl;
    //    LOGGER << windSize.transpose() << endl;
    
    
    timer.getTime();
    
    LOGGER << endl;
    displayAverageWindowSize(windSize);
    LOGGER << "LD matrix diagonal mean " << ZPZdiag.mean() << " variance " << Gadget::calcVariance(ZPZdiag) << "." << endl;
    if (ZPZdiag.mean() < 0.8 || ZPZdiag.mean() > 1.2) LOGGER << "ERROR: The mean of LD matrix diagonal values is expected to be close to one. Something is wrong with the LD matrix!" << endl;
    LOGGER << "Genotype data for " << numKeptInds << " individuals and " << numSnpInRange << " SNPs are included from [" + bedFile + "]." << endl;
    LOGGER << "Build of LD matrix completed (time used: " << timer.format(timer.getElapse()) << ")." << endl;
        
    vector<SnpInfo*> snpVecTmp(numSnpInRange);
    for (unsigned i=0; i<numSnpInRange; ++i) {
        snpVecTmp[i] = incdSnpInfoVec[start+i];
    }
    incdSnpInfoVec = snpVecTmp;
    numIncdSnps = numSnpInRange;
    string outfilename = filename;
    if (!snpRange.empty()) outfilename += ".snp" + snpRange;
    outputLDmatrix(LDmatType, outfilename, writeLdmTxt);
}

void Data::outputLDmatrix(const string &LDmatType, const string &filename, const bool writeLdmTxt) const {
    string outfilename = filename + ".ldm." + LDmatType;
    string outfile1 = outfilename + ".info";
    string outfile2 = outfilename + ".bin";
    ofstream out1(outfile1.c_str());
    FILE *out2 = fopen(outfile2.c_str(), "wb");
    ofstream out3;
    string outfile3;
    if (writeLdmTxt) {
        outfile3 = outfilename + ".txt";
        out3.open(outfile3.c_str());
    }
    out1 << boost::format("%6s %15s %10s %15s %6s %6s %12s %10s %10s %10s %10s %15s %10s %12s %12s\n")
    % "Chrom"
    % "ID"
    % "GenPos"
    % "PhysPos"
    % "A1"
    % "A2"
    % "A1Freq"
    % "Index"
    % "WindStart"
    % "WindEnd"
    % "WindSize"
    % "WindWidth"
    % "N"
    % "SamplVar"
    % "LDsum";
    SnpInfo *snp, *windStart, *windEnd;
    for (unsigned i=0; i<numIncdSnps; ++i) {
        snp = incdSnpInfoVec[i];
        windStart = incdSnpInfoVec[snp->windStart];
        windEnd = incdSnpInfoVec[snp->windEnd];
        out1 << boost::format("%6s %15s %10s %15s %6s %6s %12f %10s %10s %10s %10s %15s %10s %12.6f %12.6f\n")
        % snp->chrom
        % snp->rsID
        % snp->genPos
        % snp->physPos
        % snp->a1
        % snp->a2
        % snp->af
        % snp->index
        % snp->windStart
        % snp->windEnd
        % snp->windSize
        % (windStart->chrom == windEnd->chrom ? windEnd->physPos - windStart->physPos : windStart->chrom-windEnd->chrom)
        % numKeptInds
        % snp->ldSamplVar
        % snp->ldSum;
        if (LDmatType == "sparse") {
            fwrite(ZPZsp[i].innerIndexPtr(), sizeof(unsigned), ZPZsp[i].nonZeros(), out2);
            fwrite(ZPZsp[i].valuePtr(), sizeof(double), ZPZsp[i].nonZeros(), out2);
            if (writeLdmTxt) {
                //out3 << ZPZsp[i].transpose();
                for (SparseVector <double> ::InnerIterator it(ZPZsp[i]); it; ++it) {
                    out3 << boost::format("%-8s %-15s %-8s %-15s %-15s\n")
                    % (i+1)
                    % snp->rsID
                    % (it.index()+1)
                    % incdSnpInfoVec[it.index()]->rsID
                    % it.value();
                }
            }
        } else {
            fwrite(&ZPZ[i][0], sizeof(double), snp->windSize, out2);
            if (writeLdmTxt) out3 << ZPZ[i].transpose() << endl;
        }
    }
    out1.close();
    fclose(out2);
    
    LOGGER << "Written the SNP info into file [" << outfile1 << "]." << endl;
    LOGGER << "Written the LD matrix into file [" << outfile2 << "]." << endl;
    
    if (writeLdmTxt) {
        out3.close();
        LOGGER << "Written the LD matrix into text file [" << outfile3 << "]." << endl;
    }
}


void Data::displayAverageWindowSize(const VectorXi &windSize){
//    double windSizeMean = 0.0;
//    double windSizeSqMean = 0.0;
//    long size = windSize.size();
//    for (unsigned i=0; i<size; ++i) {
//        windSizeMean += (windSize[i] - windSizeMean)/(i+1);
//        windSizeSqMean += (windSize[i]*windSize[i] - windSizeSqMean)/(i+1);
//    }
//    LOGGER << "Per-SNP window size mean " << windSizeMean << " sd " << windSizeSqMean-windSizeMean*windSizeMean << "." << endl;
    VectorXd windSizeFloat = windSize.cast <double> ();
    double windSizeMean = Gadget::calcMean(windSizeFloat);
    double windSizeSD = sqrt(Gadget::calcVariance(windSizeFloat));
    LOGGER << "Per-SNP window size mean " << windSizeMean << " sd " << windSizeSD << "." << endl;
}

void Data::resizeWindow(const vector<SnpInfo *> &incdSnpInfoVec, const VectorXi &windStartOri, const VectorXi &windSizeOri,
                        VectorXi &windStart, VectorXi &windSize){
    bool reindexed = false;
    for (unsigned i=0; i<numSnps; ++i) {
        if (!snpInfoVec[i]->included) {
            reindexed = true;
            break;
        }
    }
    if (reindexed == false) {
        windStart = windStartOri;
        windSize  = windSizeOri;
        return;
    }
    LOGGER << "Resizing per-SNP LD window..." << endl;
    windStart.setZero(numIncdSnps);
    windSize.setZero(numIncdSnps);
    for (unsigned i=0; i<numSnps; ++i) {
        SnpInfo *snpi = snpInfoVec[i];
        if (!snpi->included) continue;
        unsigned windEndOri = windStartOri[i] + windSizeOri[i];
        for (unsigned j=windStartOri[i]; j<windEndOri; ++j) {
            SnpInfo *snpj = snpInfoVec[j];
            if (!snpj->included) continue;
            if (!windSize[snpi->index]) {
                windStart[snpi->index] = snpj->index;
            }
            ++windSize[snpi->index];
        }
        snpi->windStart = windStart[snpi->index];
        snpi->windSize  = windSize[snpi->index];
        snpi->windEnd   = snpi->windStart + snpi->windSize - 1;
    }
}

void Data::readLDmatrixInfoFileOld(const string &ldmatrixFile){   // old format: no allele frequency, no header
    ifstream in(ldmatrixFile.c_str());
    if (!in) LOGGER.e(0,"can not open the file [" + ldmatrixFile + "] to read.");
    //LOGGER << "Reading SNP info from [" + ldmatrixFile + "]." << endl;
    //snpInfoVec.clear();
    //snpInfoMap.clear();
    string header;
    string id, allele1, allele2;
    unsigned chr, physPos;
    double genPos;
    unsigned idx, windStart, windEnd, windSize, windWidth;
    long sampleSize;
    while (in >> chr >> id >> genPos >> physPos >> allele1 >> allele2 >> idx >> windStart >> windEnd >> windSize >> windWidth >> sampleSize) {
        SnpInfo *snp = new SnpInfo(idx, id, allele1, allele2, chr, genPos, physPos);
        snp->windStart = snp->windStartOri = windStart;
        snp->windEnd = snp->windEndOri = windEnd;
        snp->windSize = snp->windSizeOri = windSize;
        snp->sampleSize = sampleSize;
        snpInfoVec.push_back(snp);
        if (snpInfoMap.insert(pair<string, SnpInfo*>(id, snp)).second == false) {
            LOGGER.e(0,"Duplicate SNP ID found: \"" + id + "\".");
        }
    }
    in.close();
    numSnps = (unsigned) snpInfoVec.size();
    LOGGER << numSnps << " SNPs to be included from [" + ldmatrixFile + "]." << endl;
}

void Data::readLDmatrixInfoFile(const string &ldmatrixFile){
    ifstream in(ldmatrixFile.c_str());
    if (!in) LOGGER.e(0,"can not open the file [" + ldmatrixFile + "] to read.");
    LOGGER << "Reading SNP info from [" + ldmatrixFile + "]." << endl;
    //snpInfoVec.clear();
    //snpInfoMap.clear();
    string header;
    string id, allele1, allele2;
    unsigned chr, physPos;
    double genPos, af, ldSamplVar, ldSum;
    unsigned idx, windStart, windEnd, windSize;
    int windWidth;
    long sampleSize;
    bool skeleton;
    getline(in, header);
    Gadget::Tokenizer token;
    token.getTokens(header, " ");
    
    if (token.size() == 12) {
        in.close();
        readLDmatrixInfoFileOld(ldmatrixFile);
        return;
    }

    if (token.back() == "Skeleton") {
        while (in >> chr >> id >> genPos >> physPos >> allele1 >> allele2 >> af >> idx >> windStart >> windEnd >> windSize >> windWidth >> sampleSize >> ldSamplVar >> ldSum >> skeleton) {
            SnpInfo *snp = new SnpInfo(idx, id, allele1, allele2, chr, genPos, physPos);
            snp->af = af;
            snp->twopq = 2.0*af*(1.0-af);
            snp->windStart = snp->windStartOri = windStart;
            snp->windEnd = snp->windEndOri = windEnd;
            snp->windSize = snp->windSizeOri = windSize;
            snp->sampleSize = sampleSize;
            snp->ldSamplVar = ldSamplVar;
            snp->ldSum = ldSum;
            snp->skeleton = skeleton;
            snpInfoVec.push_back(snp);
            if (snpInfoMap.insert(pair<string, SnpInfo*>(id, snp)).second == false) {
                LOGGER.e(0,"Duplicate SNP ID found: \"" + id + "\".");
            }
        }
    }
    else if (token.back() == "LDsum") {
        while (in >> chr >> id >> genPos >> physPos >> allele1 >> allele2 >> af >> idx >> windStart >> windEnd >> windSize >> windWidth >> sampleSize >> ldSamplVar >> ldSum) {
            SnpInfo *snp = new SnpInfo(idx, id, allele1, allele2, chr, genPos, physPos);
            snp->af = af;
            snp->twopq = 2.0*af*(1.0-af);
            snp->windStart = snp->windStartOri = windStart;
            snp->windEnd = snp->windEndOri = windEnd;
            snp->windSize = snp->windSizeOri = windSize;
            snp->sampleSize = sampleSize;
            snp->ldSamplVar = ldSamplVar;
            snp->ldSum = ldSum;
            snpInfoVec.push_back(snp);
            if (snpInfoMap.insert(pair<string, SnpInfo*>(id, snp)).second == false) {
                LOGGER.e(0,"Duplicate SNP ID found: \"" + id + "\".");
            }
        }
    }
    else {
        while (in >> chr >> id >> genPos >> physPos >> allele1 >> allele2 >> af >> idx >> windStart >> windEnd >> windSize >> windWidth >> sampleSize >> ldSamplVar) {
            SnpInfo *snp = new SnpInfo(idx, id, allele1, allele2, chr, genPos, physPos);
            snp->af = af;
            snp->twopq = 2.0*af*(1.0-af);
            snp->windStart = snp->windStartOri = windStart;
            snp->windEnd = snp->windEndOri = windEnd;
            snp->windSize = snp->windSizeOri = windSize;
            snp->sampleSize = sampleSize;
            snp->ldSamplVar = ldSamplVar;
            snpInfoVec.push_back(snp);
            if (snpInfoMap.insert(pair<string, SnpInfo*>(id, snp)).second == false) {
                LOGGER.e(0,"Duplicate SNP ID found: \"" + id + "\".");
            }
        }
    }
    in.close();
    numSnps = (unsigned) snpInfoVec.size();
    LOGGER << numSnps << " SNPs to be included from [" + ldmatrixFile + "]." << endl;
}

void Data::readLDmatrixBinFile(const string &ldmatrixFile){

    Gadget::Timer timer;
    timer.setTime();
    
    Gadget::Tokenizer token;
    token.getTokens(ldmatrixFile, ".");
    string ldmType = token[token.size()-2];
    sparseLDM = ldmType == "sparse" ? true : false;

    // New part for shrunk part ldm 
    shrunkLDM = ldmType == "shrunk" ? true : false;
    
    VectorXi windStartLDM(numSnps);
    VectorXi windSizeLDM(numSnps);
    
    windStart.resize(numIncdSnps);
    windSize.resize(numIncdSnps);
    
    SnpInfo *snpi, *snpj;
    
    for (unsigned i=0; i<numSnps; ++i) {
        SnpInfo *snpi = snpInfoVec[i];
        windStartLDM[i] = snpi->windStart;
        windSizeLDM[i]  = snpi->windSize;
    }

    FILE *in = fopen(ldmatrixFile.c_str(), "rb");
    if (!in) {
        LOGGER.e(0," cannot open LD matrix file " + ldmatrixFile);
    }
    
    if (!sparseLDM) resizeWindow(incdSnpInfoVec, windStartLDM, windSizeLDM, windStart, windSize);
    
    if (numIncdSnps == 0) LOGGER.e(0,"No SNP is retained for analysis.");
    
    LOGGER << "Reading " + ldmType + " LD matrix from [" + ldmatrixFile + "]..." << endl;
    
    double rsq = 0.0;
    
    if (sparseLDM) {
        ZPZsp.resize(numIncdSnps);
        ZPZdiag.resize(numIncdSnps);
       
        for (unsigned i=0, inci=0; i<numSnps; i++) {
            snpi = snpInfoVec[i];
                        
            unsigned d[windSizeLDM[i]];
            double v[windSizeLDM[i]];
            
            if (!snpi->included) {
                fseek(in, sizeof(d), SEEK_CUR);
                fseek(in, sizeof(v), SEEK_CUR);
                continue;
            }
            
            fread(d, sizeof(d), 1, in);
            fread(v, sizeof(v), 1, in);
            
            ZPZsp[inci].resize(windSizeLDM[i]);
            snpi->ldSamplVar = 0.0;
            snpi->ldSum = 0.0;
            if (!readLDscore) snpi->ldsc = 0.0;
            
            for (unsigned j=0; j<windSizeLDM[i]; ++j) {
                snpj = snpInfoVec[d[j]];
        
                if (snpj->included) {
                    ZPZsp[inci].insertBack(snpj->index) = v[j];
                    rsq = v[j]*v[j];
                    snpi->ldSamplVar += (1.0-rsq)*(1.0-rsq)/snpi->sampleSize;
                    snpi->ldSum += v[j];
                    if (!readLDscore) snpi->ldsc += rsq;
                    if (snpj == snpi)
                        ZPZdiag[inci] = v[j];
                }
            }
            SparseVector <double> ::InnerIterator it(ZPZsp[inci]);
            windStart[inci] = snpi->windStart = it.index();
            windSize[inci] = snpi->windSize = ZPZsp[inci].nonZeros();
            for (; it; ++it) snpi->windEnd = it.index();
            snpi->numNonZeroLD = snpi->windSize;
            
            if (++inci == numIncdSnps) break;
        }
    }
    else {
        ZPZ.resize(numIncdSnps);
        ZPZdiag.resize(numIncdSnps);
        
        for (unsigned i=0, inci=0; i<numSnps; i++) {
            snpi = snpInfoVec[i];
            
            double v[windSizeLDM[i]];
            
            if (!snpi->included) {
                fseek(in, sizeof(v), SEEK_CUR);
                continue;
            }
            
            fread(v, sizeof(v), 1, in);
            
            ZPZ[inci].resize(windSize[inci]);
            snpi->ldSamplVar = 0.0;
            snpi->ldSum = 0.0;
            if (!readLDscore) snpi->ldsc = 0.0;
            snpi->numNonZeroLD = snpi->windSize;
            
            for (unsigned j=0, incj=0; j<windSizeLDM[i]; ++j) {
                snpj = snpInfoVec[windStartLDM[i]+j];
                if (snpj->included) {
                    ZPZ[inci][incj++] = v[j];
                    rsq = v[j]*v[j];
                    snpi->ldSamplVar += (1.0-rsq)*(1.0-rsq)/snpi->sampleSize;
                    snpi->ldSum += v[j];
                    if (!readLDscore) snpi->ldsc += rsq;
                    if (snpj == snpi)
                        ZPZdiag[inci] = v[j];
                }
            }
            
            if (++inci == numIncdSnps) break;
        }
    }
    
    fclose(in);
    
    timer.getTime();
    
//    LOGGER << "Window width " << windowWidth << " Mb." << endl;
    displayAverageWindowSize(windSize);
    LOGGER << "LD matrix diagonal mean " << ZPZdiag.mean() << " sd " << sqrt(Gadget::calcVariance(ZPZdiag)) << "." << endl;
    if (ZPZdiag.mean() < 0.8 || ZPZdiag.mean() > 1.2) LOGGER.e(0," The mean of LD matrix diagonal values is expected to be close to one. Something is wrong with the LD matrix!");
    LOGGER << "Read LD matrix for " << numIncdSnps << " SNPs (time used: " << timer.format(timer.getElapse()) << ")." << endl;
}

void Data::readLDmatrixBinFileAndShrink(const string &ldmatrixFile){
    
    Gadget::Timer timer;
    timer.setTime();
    
    Gadget::Tokenizer token;
    token.getTokens(ldmatrixFile, ".");
    string ldmType = token[token.size()-2];
    sparseLDM = ldmType == "sparse" ? true : false;
    
    VectorXi windStartLDM(numSnps);
    VectorXi windSizeLDM(numSnps);
    
    windStart.resize(numIncdSnps);
    windSize.resize(numIncdSnps);
    
    SnpInfo *snpi, *snpj;
    
    for (unsigned i=0; i<numSnps; ++i) {
        SnpInfo *snpi = snpInfoVec[i];
        windStartLDM[i] = snpi->windStart;
        windSizeLDM[i]  = snpi->windSize;
    }
    
    FILE *in = fopen(ldmatrixFile.c_str(), "rb");
    if (!in) {
        LOGGER.e(0," cannot open LD matrix file " + ldmatrixFile);
    }
    
    if (!sparseLDM) resizeWindow(incdSnpInfoVec, windStartLDM, windSizeLDM, windStart, windSize);
    
    if (numIncdSnps == 0) LOGGER.e(0,"No SNP is retained for analysis.");
    
    LOGGER << "Reading and shrinking " + ldmType + " LD matrix from [" + ldmatrixFile + "]..." << endl;
    
    double rsq = 0.0;
    
//    double nref = incdSnpInfoVec[0]->sampleSize;
    double nref = 183;  // this is the sample size for genetic map
    double n = 2.0*nref-1.0;
    // Approximation to the harmonic series
    double nsum = log(n) + 0.5772156649 + 1.0 / (2.0 * n) - 1.0 / (12.0 * pow(n, 2.0)) + 1.0 / (120.0 * pow(n, 4.0));
    double theta = (1.0 / nsum) / (2.0 * (nref) + 1.0 / nsum);
    double off  = (1.0 - theta)*(1.0 - theta);
    double diag = 0.5*theta*(1.0 - 0.5*theta);
    
    double shrunkLD;
    double rho;
    double Ne = 11490.672741;
    double shrinkage;
    double sdprod;
    
    if (sparseLDM) {
        ZPZsp.resize(numIncdSnps);
        ZPZdiag.resize(numIncdSnps);
        
        for (unsigned i=0, inci=0; i<numSnps; i++) {
            snpi = snpInfoVec[i];
            
            unsigned d[windSizeLDM[i]];
            double v[windSizeLDM[i]];
            
            if (!snpi->included) {
                fseek(in, sizeof(d), SEEK_CUR);
                fseek(in, sizeof(v), SEEK_CUR);
                continue;
            }
            
            fread(d, sizeof(d), 1, in);
            fread(v, sizeof(v), 1, in);
            
            ZPZsp[inci].resize(windSizeLDM[i]);
            snpi->ldSamplVar = 0.0;
            snpi->ldSum = 0.0;
            if (!readLDscore) snpi->ldsc = 0.0;
            
            for (unsigned j=0; j<windSizeLDM[i]; ++j) {
                snpj = snpInfoVec[d[j]];
                if (snpj->included) {
                    
                    sdprod = sqrt(snpi->twopq*snpj->twopq);
                    shrunkLD = v[j]*sdprod;
                    rho = 4.0 * Ne * abs(snpi->genPos - snpj->genPos)/100.0;
                    shrinkage = exp(-rho / (2.0*nref)) * off;
                    if (shrinkage <= 1e-5) shrinkage = 0.0;
                    shrunkLD *= shrinkage;
                    shrunkLD /= sdprod;
                    if (snpj == snpi) {
                        shrunkLD += diag/sdprod;
                        ZPZdiag[inci] = shrunkLD;
                    }
                    
                    ZPZsp[inci].insertBack(snpj->index) = shrunkLD;
                    rsq = shrunkLD * shrunkLD;
                    snpi->ldSamplVar += shrinkage*shrinkage*(1.0-rsq)*(1.0-rsq)/snpi->sampleSize;
                    snpi->ldSum += shrunkLD;
                    if (!readLDscore) snpi->ldsc += rsq;
                }
            }
            SparseVector <double> ::InnerIterator it(ZPZsp[inci]);
            windStart[inci] = snpi->windStart = it.index();
            windSize[inci] = snpi->windSize = ZPZsp[inci].nonZeros();
            for (; it; ++it) snpi->windEnd = it.index();
            snpi->numNonZeroLD = snpi->windSize;
            
            if (++inci == numIncdSnps) break;
        }
    }
    else {
        ZPZ.resize(numIncdSnps);
        ZPZdiag.resize(numIncdSnps);
        
        for (unsigned i=0, inci=0; i<numSnps; i++) {
            snpi = snpInfoVec[i];
            
            double v[windSizeLDM[i]];
            
            if (!snpi->included) {
                fseek(in, sizeof(v), SEEK_CUR);
                continue;
            }
            
            fread(v, sizeof(v), 1, in);
            
            ZPZ[inci].resize(windSize[inci]);
            snpi->ldSamplVar = 0.0;
            snpi->ldSum = 0.0;
            if (!readLDscore) snpi->ldsc = 0.0;
            snpi->numNonZeroLD = 0;
            
            for (unsigned j=0, incj=0; j<windSizeLDM[i]; ++j) {
                snpj = snpInfoVec[windStartLDM[i]+j];
                if (snpj->included) {
                    
                    sdprod = sqrt(snpi->twopq*snpj->twopq);
                    shrunkLD = v[j]*sdprod;
                    rho = 4.0 * Ne * abs(snpi->genPos - snpj->genPos)/100.0;
                    shrinkage = exp(-rho / (2.0*nref)) * off;
                    if (shrinkage <= 1e-5) shrinkage = 0.0;
                    else snpi->numNonZeroLD++;
                    shrunkLD *= shrinkage;
                    shrunkLD /= sdprod;
                    if (snpj == snpi) {
                        shrunkLD += diag/sdprod;
                        ZPZdiag[inci] = shrunkLD;
                    }
                    
                    ZPZ[inci][incj++] = shrunkLD;
                    rsq = shrunkLD * shrunkLD;
                    snpi->ldSamplVar += shrinkage*shrinkage*(1.0-rsq)*(1.0-rsq)/snpi->sampleSize;
                    snpi->ldSum += shrunkLD;
                    if (!readLDscore) snpi->ldsc += rsq;
                }
            }
            
            if (++inci == numIncdSnps) break;
        }
    }
    
    fclose(in);
    
    timer.getTime();
    
    //    LOGGER << "Window width " << windowWidth << " Mb." << endl;
    displayAverageWindowSize(windSize);
    LOGGER << "LD matrix diagonal mean " << ZPZdiag.mean() << " sd " << sqrt(Gadget::calcVariance(ZPZdiag)) << "." << endl;
    if (ZPZdiag.mean() < 0.8 || ZPZdiag.mean() > 1.2) LOGGER.e(0," The mean of LD matrix diagonal values is expected to be close to one. Something is wrong with the LD matrix!");
    LOGGER << "Read LD matrix for " << numIncdSnps << " SNPs (time used: " << timer.format(timer.getElapse()) << ")." << endl;
}

void Data::readMultiLDmatInfoFile(const string &mldmatFile){
    ifstream in(mldmatFile.c_str());
    if (!in) LOGGER.e(0,"can not open the file [" + mldmatFile + "] to read.");
    LOGGER << "Reading SNP info from [" + mldmatFile + "]..." << endl;
    string inputStr;
    numSnpMldVec.clear();
    while (getline(in, inputStr)) {
        readLDmatrixInfoFile(inputStr+".info");
        numSnpMldVec.push_back(numSnps);
    }
    SnpInfo *snp = snpInfoVec[numSnpMldVec[0]];
    if (snp->index == 0) reindexed = true;
    else reindexed = false;
}

void Data::readMultiLDmatBinFile(const string &mldmatFile){
    //LOGGER << "Hi I'm in here " << endl;
    vector<string> filenameVec;
    ifstream in1(mldmatFile.c_str());
    if (!in1) LOGGER.e(0,"can not open the file [" + mldmatFile + "] to read.");
    
    Gadget::Timer timer;
    timer.setTime();
    
    string inputStr;
    string ldmType;
    sparseLDM = true;
    while (getline(in1, inputStr)) {
        filenameVec.push_back(inputStr + ".bin");
        Gadget::Tokenizer token;
        token.getTokens(inputStr, ".");
        ldmType = token[token.size()-1];
        sparseLDM = ldmType == "sparse" ? true : false;
    }

    LOGGER << "Reading " + ldmType + " LD matrices from [" + mldmatFile + "]..." << endl;

    VectorXi windStartLDM(numSnps);
    VectorXi windSizeLDM(numSnps);
    //LOGGER << "I made it here " << endl; 
    for (unsigned j=0, i=0, cnt=0; j<numSnps; ++j) {
        SnpInfo *snp = snpInfoVec[j];
        if (j==numSnpMldVec[i]) {
            if (reindexed)
                cnt = numSnpMldVec[i++];
            else
                cnt = 0;
        }
        snp->windStart += cnt;
        snp->windEnd   += cnt;
        windSizeLDM[j]  = snp->windSize;
        windStartLDM[j] = snp->windStart;
    }
    
    if (sparseLDM) {
        windStart.setZero(numIncdSnps);
        windSize.setZero(numIncdSnps);
        ZPZsp.resize(numIncdSnps);
    }
    else {
        resizeWindow(incdSnpInfoVec, windStartLDM, windSizeLDM, windStart, windSize);
        ZPZ.resize(numIncdSnps);
    }
    ZPZdiag.resize(numIncdSnps);
    
    unsigned starti = 0;
    unsigned incj = 0;
    //LOGGER << "I made it here " << endl;    
    long numFiles = filenameVec.size();
    for (unsigned i=0; i<numFiles; ++i) {
        FILE *in2 = fopen(filenameVec[i].c_str(), "rb");
        if (!in2) {
            LOGGER.e(0," cannot open LD matrix file " + filenameVec[i]);
        }
        
        SnpInfo *snpj = NULL;
        SnpInfo *snpk = NULL;
        
        double rsq = 0.0;
        
        if (sparseLDM) {
            for (unsigned j=starti; j<numSnpMldVec[i]; j++) {
                snpj = snpInfoVec[j];
                
                unsigned d[windSizeLDM[j]];
                double v[windSizeLDM[j]];
                
                if (!snpj->included) {
                    fseek(in2, sizeof(d), SEEK_CUR);
                    fseek(in2, sizeof(v), SEEK_CUR);
                    continue;
                }
                
                fread(d, sizeof(d), 1, in2);
                fread(v, sizeof(v), 1, in2);
                
                ZPZsp[incj].resize(windSizeLDM[j]);
                snpj->ldSamplVar = 0.0;
                snpj->ldSum = 0.0;
                if (!readLDscore) snpj->ldsc = 0.0;

                for (unsigned k=0; k<windSizeLDM[j]; ++k) {
                    snpk = snpInfoVec[windStartLDM[j]+d[k]-d[0]];
                    if (snpk->included) {
                        ZPZsp[incj].insertBack(snpk->index) = v[k];
                        rsq = v[k]*v[k];
                        snpj->ldSamplVar += (1.0-rsq)*(1.0-rsq)/snpj->sampleSize;
                        snpj->ldSum += v[k];
                        if (!readLDscore) snpj->ldsc += rsq;
                        if (snpk == snpj)
                            ZPZdiag[incj] = v[k];
                    }
                }
                SparseVector <double> ::InnerIterator it(ZPZsp[incj]);
                windStart[incj] = snpj->windStart = it.index();
                windSize[incj] = snpj->windSize = ZPZsp[incj].nonZeros();
                snpj->numNonZeroLD = snpj->windSize;
                ++incj;
            }
        }
        else {
            //LOGGER << "I's reading the snps " << endl;
            for (unsigned j=starti; j<numSnpMldVec[i]; j++) {
                snpj = snpInfoVec[j];
                double v[windSizeLDM[j]];
	        //LOGGER << "I'm at SNP " << j << " windStartLDM[j] " << windStartLDM[j] << " numSnpMldVec[i] " << numSnpMldVec[i] << " windSizeLDM[j] " << windSizeLDM[j] << endl;
                if (!snpj->included) {
                    fseek(in2, sizeof(v), SEEK_CUR);
                    continue;
                }
                
                fread(v, sizeof(v), 1, in2);
                //LOGGER << "Here 1 " << endl;                 
                ZPZ[incj].resize(windSize[incj]);
                snpj->ldSamplVar = 0.0;
                snpj->ldSum = 0.0;
                if (!readLDscore) snpj->ldsc = 0.0;
                snpj->numNonZeroLD = snpj->windSize;

                for (unsigned k=0, inck=0; k<windSizeLDM[j]; ++k) {
                    //LOGGER << "I'm at k " << k << endl; 
                    snpk = snpInfoVec[windStartLDM[j]+k];
                    if (snpk->included) {
                        ZPZ[incj][inck++] = v[k];
                        rsq = v[k]*v[k];
                        snpj->ldSamplVar += (1.0-rsq)*(1.0-rsq)/snpj->sampleSize;
                        snpj->ldSum += v[k];
                        if (!readLDscore) snpj->ldsc += rsq;
                        if (snpk == snpj)
                            ZPZdiag[incj] = v[k];
                    }
                }
                 //LOGGER << "Here 2 " << endl;
                ++incj;
            }
        }
        
        fclose(in2);
        LOGGER << "Read " + ldmType + " LD matrix for " << numSnpMldVec[i]-starti << " SNPs from [" << filenameVec[i] << "]." << endl;
        
        starti = numSnpMldVec[i];
    }
    
    timer.getTime();
    
    displayAverageWindowSize(windSize);
    LOGGER << "LD matrix diagonal mean " << ZPZdiag.mean() << " sd " << sqrt(Gadget::calcVariance(ZPZdiag)) << "." << endl;
    if (ZPZdiag.mean() < 0.8 || ZPZdiag.mean() > 1.2) LOGGER.e(0," The mean of LD matrix diagonal values is expected to be close to one. Something is wrong with the LD matrix!");
    LOGGER << "Read LD matrix for " << numIncdSnps << " SNPs (time used: " << timer.format(timer.getElapse()) << ")." << endl;
}

void Data::readMultiLDmatBinFileAndShrink(const string &mldmatFile, const double genMapN){
    //LOGGER << "Hi I'm in here " << endl;
    vector<string> filenameVec;
    ifstream in1(mldmatFile.c_str());
    if (!in1) LOGGER.e(0,"can not open the file [" + mldmatFile + "] to read.");
    
    Gadget::Timer timer;
    timer.setTime();
    
    string inputStr;
    string ldmType;
    sparseLDM = true;
    while (getline(in1, inputStr)) {
        filenameVec.push_back(inputStr + ".bin");
        Gadget::Tokenizer token;
        token.getTokens(inputStr, ".");
        ldmType = token[token.size()-1];
        sparseLDM = ldmType == "sparse" ? true : false;
    }
    
    LOGGER << "Reading and shrinking " + ldmType + " LD matrices from [" + mldmatFile + "]..." << endl;
    
    VectorXi windStartLDM(numSnps);
    VectorXi windSizeLDM(numSnps);
    //LOGGER << "I made it here " << endl;
    for (unsigned j=0, i=0, cnt=0; j<numSnps; ++j) {
        SnpInfo *snp = snpInfoVec[j];
        if (j==numSnpMldVec[i]) {
            if (reindexed)
                cnt = numSnpMldVec[i++];
            else
                cnt = 0;
        }
        snp->windStart += cnt;
        snp->windEnd   += cnt;
        windSizeLDM[j]  = snp->windSize;
        windStartLDM[j] = snp->windStart;
    }
    
    if (sparseLDM) {
        windStart.setZero(numIncdSnps);
        windSize.setZero(numIncdSnps);
        ZPZsp.resize(numIncdSnps);
    }
    else {
        resizeWindow(incdSnpInfoVec, windStartLDM, windSizeLDM, windStart, windSize);
        ZPZ.resize(numIncdSnps);
    }
    ZPZdiag.resize(numIncdSnps);
    
    unsigned starti = 0;
    unsigned incj = 0;
    //LOGGER << "I made it here " << endl;
    long numFiles = filenameVec.size();
    for (unsigned i=0; i<numFiles; ++i) {
        FILE *in2 = fopen(filenameVec[i].c_str(), "rb");
        if (!in2) {
            LOGGER.e(0," cannot open LD matrix file " + filenameVec[i]);
        }
        
        SnpInfo *snpj = NULL;
        SnpInfo *snpk = NULL;
        
        double rsq = 0.0;
        
        double nref = genMapN;
        double n = 2.0*nref-1.0;
        // Approximation to the harmonic series
        double nsum = log(n) + 0.5772156649 + 1.0 / (2.0 * n) - 1.0 / (12.0 * pow(n, 2.0)) + 1.0 / (120.0 * pow(n, 4.0));
        double theta = (1.0 / nsum) / (2.0 * (nref) + 1.0 / nsum);
        double off  = (1.0 - theta)*(1.0 - theta);
        double diag = 0.5*theta*(1.0 - 0.5*theta);
        
        double shrunkLD;
        double rho;
        double Ne = 11490.672741;
        double shrinkage;
        double sdprod;

        if (sparseLDM) {
            for (unsigned j=starti; j<numSnpMldVec[i]; j++) {
                snpj = snpInfoVec[j];
                
                unsigned d[windSizeLDM[j]];
                double v[windSizeLDM[j]];
                
                if (!snpj->included) {
                    fseek(in2, sizeof(d), SEEK_CUR);
                    fseek(in2, sizeof(v), SEEK_CUR);
                    continue;
                }
                
                fread(d, sizeof(d), 1, in2);
                fread(v, sizeof(v), 1, in2);
                
                ZPZsp[incj].resize(windSizeLDM[j]);
                snpj->ldSamplVar = 0.0;
                snpj->ldSum = 0.0;
                if (!readLDscore) snpj->ldsc = 0.0;
                
                for (unsigned k=0; k<windSizeLDM[j]; ++k) {
                    snpk = snpInfoVec[windStartLDM[j]+d[k]-d[0]];
                    if (snpk->included) {
                        
                        sdprod = sqrt(snpj->twopq*snpk->twopq);
                        shrunkLD = v[k]*sdprod;
                        rho = 4.0 * Ne * abs(snpj->genPos - snpk->genPos)/100.0;
                        shrinkage = exp(-rho / (2.0*nref)) * off;
                        if (shrinkage <= 1e-5) shrinkage = 0.0;
                        shrunkLD *= shrinkage;
                        shrunkLD /= sdprod;
                        if (snpk == snpj) {
                            shrunkLD += diag/sdprod;
                            ZPZdiag[incj] = shrunkLD;
                        }
                        
                        ZPZsp[incj].insertBack(snpk->index) = shrunkLD;
                        rsq = shrunkLD * shrunkLD;
                        snpj->ldSamplVar += shrinkage*shrinkage*(1.0-rsq)*(1.0-rsq)/snpj->sampleSize;
                        snpj->ldSum += shrunkLD;
                        if (!readLDscore) snpj->ldsc += rsq;
                    }
                }
                SparseVector <double> ::InnerIterator it(ZPZsp[incj]);
                windStart[incj] = snpj->windStart = it.index();
                windSize[incj] = snpj->windSize = ZPZsp[incj].nonZeros();
                snpj->numNonZeroLD = snpj->windSize;
                ++incj;
            }
        }
        else {
            //LOGGER << "I's reading the snps " << endl;
            for (unsigned j=starti; j<numSnpMldVec[i]; j++) {
                snpj = snpInfoVec[j];
                double v[windSizeLDM[j]];
                //LOGGER << "I'm at SNP " << j << " windStartLDM[j] " << windStartLDM[j] << " numSnpMldVec[i] " << numSnpMldVec[i] << " windSizeLDM[j] " << windSizeLDM[j] << endl;
                if (!snpj->included) {
                    fseek(in2, sizeof(v), SEEK_CUR);
                    continue;
                }
                
                fread(v, sizeof(v), 1, in2);
                //LOGGER << "Here 1 " << endl;
                ZPZ[incj].resize(windSize[incj]);
                snpj->ldSamplVar = 0.0;
                snpj->ldSum = 0.0;
                if (!readLDscore) snpj->ldsc = 0.0;
                snpj->numNonZeroLD = 0;
                
                for (unsigned k=0, inck=0; k<windSizeLDM[j]; ++k) {
                    //LOGGER << "I'm at k " << k << endl;
                    snpk = snpInfoVec[windStartLDM[j]+k];
                    if (snpk->included) {
                        
                        sdprod = sqrt(snpj->twopq*snpk->twopq);
                        shrunkLD = v[k]*sdprod;
                        rho = 4.0 * Ne * abs(snpj->genPos - snpk->genPos)/100.0;
                        shrinkage = exp(-rho / (2.0*nref)) * off;
                        if (shrinkage <= 1e-5) shrinkage = 0.0;
                        else snpj->numNonZeroLD++;
                        shrunkLD *= shrinkage;
                        shrunkLD /= sdprod;
                        if (snpk == snpj) {
                            shrunkLD += diag/sdprod;
                            ZPZdiag[incj] = shrunkLD;
                        }
                        
                        ZPZ[incj][inck++] = shrunkLD;
                        rsq = shrunkLD * shrunkLD;
                        snpj->ldSamplVar += shrinkage*shrinkage*(1.0-rsq)*(1.0-rsq)/snpj->sampleSize;
                        snpj->ldSum += shrunkLD;
                        if (!readLDscore) snpj->ldsc += rsq;
                    }
                }
                //LOGGER << "Here 2 " << endl;
                ++incj;
            }
        }
        
        fclose(in2);
        LOGGER << "Read " + ldmType + " LD matrix for " << numSnpMldVec[i]-starti << " SNPs from [" << filenameVec[i] << "]." << endl;
        
        starti = numSnpMldVec[i];
    }
    
    timer.getTime();
    
    displayAverageWindowSize(windSize);
    LOGGER << "LD matrix diagonal mean " << ZPZdiag.mean() << " sd " << sqrt(Gadget::calcVariance(ZPZdiag)) << "." << endl;
    if (ZPZdiag.mean() < 0.8 || ZPZdiag.mean() > 1.2) LOGGER.e(0," The mean of LD matrix diagonal values is expected to be close to one. Something is wrong with the LD matrix!");
    LOGGER << "Read LD matrix for " << numIncdSnps << " SNPs (time used: " << timer.format(timer.getElapse()) << ")." << endl;
}

void Data::resizeLDmatrix(const string &LDmatType, const double chisqThreshold, const unsigned windowWidth, const double LDthreshold, const double effpopNE, const double cutOff, const double genMapN) {
    if (LDmatType == "full") {
//        for (unsigned i=0; i<numIncdSnps; ++i) {
//            SnpInfo *snp = incdSnpInfoVec[i];
//            VectorXd rsq = ZPZ[i].cwiseProduct(ZPZ[i]);
//            VectorXd rsq_adj = rsq - (VectorXd::Ones(snp->windSize) - rsq)/double(snp->sampleSize-2);
//            snp->ldSamplVar = (VectorXd::Ones(snp->windSize) - rsq_adj).squaredNorm()/snp->sampleSize;
//            snp->ldSum = ZPZ[i].sum();
//        }
        return;
    }
    //if (LDmatType == "shrunk") return;  // TMP; to be removed
    snp2pq.resize(numIncdSnps);
    for (unsigned i=0; i<numIncdSnps; ++i) {
        SnpInfo *snp = incdSnpInfoVec[i];
        snp2pq[i] = snp->twopq = 2.0*snp->af*(1.0-snp->af);
    }
    double rsq = 0.0;    
    if (LDmatType == "sparse") {
        if (ZPZsp.size() == 0) {
            LOGGER << "Making a sparse LD matrix by setting the non-significant LD to be zero..." << endl;
            ZPZsp.resize(numIncdSnps);
            SnpInfo *snpi, *snpj;
            if (LDthreshold) {
                for (unsigned i=0; i<numIncdSnps; ++i) {
                    snpi = incdSnpInfoVec[i];
                    ZPZsp[i].resize(snpi->windSize);
                    snpi->ldSamplVar = 0.0;
                    snpi->ldSum = 0.0;
                    for (unsigned j=0; j<snpi->windSize; ++j) {
                        snpj = incdSnpInfoVec[snpi->windStart + j];
                        if (abs(ZPZ[i][j]) > LDthreshold || snpj->skeleton) {
                            ZPZsp[i].insertBack(snpi->windStart + j) = ZPZ[i][j];
                            rsq = ZPZ[i][j]*ZPZ[i][j];
                            snpi->ldSamplVar += (1.0-rsq)*(1.0-rsq)/snpi->sampleSize;
                            snpi->ldSum += ZPZ[i][j];
                        }
                    }
                    SparseVector <double> ::InnerIterator it(ZPZsp[i]);
                    windStart[i] = snpi->windStart = it.index();
                    windSize[i] = snpi->windSize = ZPZsp[i].nonZeros();
                    for (; it; ++it) snpi->windEnd = it.index();
                    ZPZ[i].resize(0);
                }
            } else {
                if (windowWidth) {
                    for (unsigned i=0; i<numIncdSnps; ++i) {
                        snpi = incdSnpInfoVec[i];
                        ZPZsp[i].resize(snpi->windSize);
                        snpi->ldSamplVar = 0.0;
                        snpi->ldSum = 0.0;
                        for (unsigned j=0; j<snpi->windSize; ++j) {
                            snpj = incdSnpInfoVec[snpi->windStart + j];
                            if (i==j || ZPZ[i][j]*ZPZ[i][j]*snpi->sampleSize > chisqThreshold ||
                                snpi->isProximal(*incdSnpInfoVec[snpi->windStart + j], windowWidth/2)) {
                                ZPZsp[i].insertBack(snpi->windStart + j) = ZPZ[i][j];
                                rsq = ZPZ[i][j]*ZPZ[i][j];
                                snpi->ldSamplVar += (1.0-rsq)*(1.0-rsq)/snpi->sampleSize;
                                snpi->ldSum += ZPZ[i][j];
                            }
                        }
                        //            ZPZsp[i] = ZPZ[i].sparseView();
                        SparseVector <double> ::InnerIterator it(ZPZsp[i]);
                        windStart[i] = snpi->windStart = it.index();
                        windSize[i] = snpi->windSize = ZPZsp[i].nonZeros();
                        for (; it; ++it) snpi->windEnd = it.index();
                        ZPZ[i].resize(0);
                        //LOGGER << i << " windsize " << snp->windSize << " " << ZPZsp[i].size() << endl;
                    }
                } else {
                    LOGGER << "Using a chisq threshold of " << chisqThreshold << endl; 
                    for (unsigned i=0; i<numIncdSnps; ++i) {
                        snpi = incdSnpInfoVec[i];
                        ZPZsp[i].resize(snpi->windSize);
                        snpi->ldSamplVar = 0.0;
                        snpi->ldSum = 0.0;
                        for (unsigned j=0; j<snpi->windSize; ++j) {
                            snpj = incdSnpInfoVec[snpi->windStart + j];
                            if (i==j || ZPZ[i][j]*ZPZ[i][j]*snpi->sampleSize > chisqThreshold) {
                                ZPZsp[i].insertBack(snpi->windStart + j) = ZPZ[i][j];
                                rsq = ZPZ[i][j]*ZPZ[i][j];
                                snpi->ldSamplVar += (1.0-rsq)*(1.0-rsq)/snpi->sampleSize;
                                snpi->ldSum += ZPZ[i][j];
                            }
                        }
                        //            ZPZsp[i] = ZPZ[i].sparseView();
                        SparseVector <double> ::InnerIterator it(ZPZsp[i]);
                        windStart[i] = snpi->windStart = it.index();
                        windSize[i] = snpi->windSize = ZPZsp[i].nonZeros();
                        for (; it; ++it) snpi->windEnd = it.index();
                        ZPZ[i].resize(0);
                        //LOGGER << i << " windsize " << snp->windSize << " " << ZPZsp[i].size() << endl;
                    }
                }
            }
        } else {
            LOGGER << "Pruning a sparse LD matrix by chisq threshold of " << chisqThreshold << endl;
            SnpInfo *snpi, *snpj;
            for (unsigned i=0; i<numIncdSnps; ++i) {
                snpi = incdSnpInfoVec[i];
                snpi->ldSamplVar = 0.0;
                snpi->ldSum = 0.0;
                for (SparseVector <double> ::InnerIterator it(ZPZsp[i]); it; ++it) {
                    snpj = incdSnpInfoVec[it.index()];
                    rsq = it.value()*it.value();
                    if (rsq*snpi->sampleSize <= chisqThreshold) it.valueRef() = 0.0;
                    else {
                        snpi->ldSamplVar += (1.0-rsq)*(1.0-rsq)/snpi->sampleSize;
                        snpi->ldSum += it.value();
                    }
                }
                ZPZsp[i].prune(0.0);
                SparseVector <double> ::InnerIterator it(ZPZsp[i]);
                windStart[i] = snpi->windStart = it.index();
                windSize[i] = snpi->windSize = ZPZsp[i].nonZeros();
                for (; it; ++it) snpi->windEnd = it.index();
            }
        }
    }
    if (LDmatType == "band") {
        VectorXi windStartOri = windStart;
        VectorXi windEndi;
        windEndi.resize(numIncdSnps);
        VectorXi windSizeOri = windSize;
        VectorXd ZPZiTmp;
        if (windowWidth) {
            LOGGER << "Resizing LD matrix based on a window width of " << windowWidth*1e-6 << " Mb..." << endl;
            getWindowInfo(incdSnpInfoVec, windowWidth, windStart, windSize);
            for (unsigned i=0; i<numIncdSnps; ++i) {
                SnpInfo *snp = incdSnpInfoVec[i];
                windStart[i] = snp->windStart = max(windStart[i], windStartOri[i]);
                windSize[i]  = snp->windSize  = min(windSize[i], windSizeOri[i]);
                ZPZiTmp = ZPZ[i].segment(windStart[i]-windStartOri[i], windSize[i]);
                ZPZ[i] = ZPZiTmp;
                snp->ldSamplVar = (1.0 - ZPZ[i].array().square()).square().sum()/snp->sampleSize;
                snp->ldSum = ZPZ[i].sum();
            }
        } else if (LDthreshold) {
            LOGGER << "Resizing LD matrix based on a LD threshold of " << LDthreshold << "..." << endl;
            for (unsigned i=0; i<numIncdSnps; ++i) {
                windEndi[i] = windSizeOri[i];
                // Gather up the window ends
                for (unsigned j=windSizeOri[i]; j>i; --j) {
                    if (abs(ZPZ[i][j-1]) > LDthreshold) {
                        windEndi[i] = j;
                        break;
                    }
                }
                // Fill the lower triangle with the necessary zeroes
                //if (windEndi != numIncdSnps) {
                for (unsigned j=windEndi[i]; j<numIncdSnps; ++j) {
                      ZPZ[j][i] = 0.0;
                }
                //}
            }
            // Cycle again over the vectors and get the window start
            for (unsigned i=0; i<numIncdSnps; ++i) {
                SnpInfo *snp = incdSnpInfoVec[i];
                // unsigned windEndi = windSizeOri[i];
                // Run over the lower triangle and make equal
                for (unsigned j=0; j<(i-1); ++j) {
                    if (abs(ZPZ[i][j]) > 0.0) {
                        windStart[i] = snp->windStart = windStartOri[i] + j;
                        break;
                    }
                }
                // Resize the vectors of vectors accordingly
                // LOGGER << "Wind start, Wind end " << windStart[i] << ", "<< windEndi[i] << endl;
                windSize[i] = snp->windSize = windStartOri[i] + windEndi[i] - windStart[i];
                ZPZiTmp = ZPZ[i].segment(windStart[i] - windStartOri[i], windSize[i]);
                ZPZ[i] = ZPZiTmp;
            }
        } else {
            LOGGER << "Resizing LD matrix based on a chisquare threshold of " << chisqThreshold << "..." << endl;
            for (unsigned i=0; i<numIncdSnps; ++i) {
                windEndi[i] = windSizeOri[i];
                 SnpInfo *snp = incdSnpInfoVec[i];
                // Gather up the window ends
                for (unsigned j=windSizeOri[i]; j>i; --j) {
                    if (ZPZ[i][j]*ZPZ[i][j]*snp->sampleSize > chisqThreshold) {
                        windEndi[i] = j;
                        break;
                    }
                }
                // Fill the lower triangle with the necessary zeroes
                //if (windEndi != numIncdSnps) {
                for (unsigned j=windEndi[i]; j<numIncdSnps; ++j) {
                      ZPZ[j][i] = 0.0;
                }
                //}
            }
            // Cycle again over the vectors and get the window start
            for (unsigned i=0; i<numIncdSnps; ++i) {
                SnpInfo *snp = incdSnpInfoVec[i];
                // unsigned windEndi = windSizeOri[i];
                // Run over the lower triangle and make equal
                for (unsigned j=0; j<(i-1); ++j) {
                    if (abs(ZPZ[i][j]) > 0.0) {
                        windStart[i] = snp->windStart = windStartOri[i] + j;
                        break;
                    }
                }
                // Resize the vectors of vectors accordingly
                // LOGGER << "Wind start, Wind end " << windStart[i] << ", "<< windEndi[i] << endl;
                windSize[i] = snp->windSize = windStartOri[i] + windEndi[i] - windStart[i];
                ZPZiTmp = ZPZ[i].segment(windStart[i] - windStartOri[i], windSize[i]);
                ZPZ[i] = ZPZiTmp;
                snp->ldSamplVar = (1.0 - ZPZ[i].array().square()).square().sum()/snp->sampleSize;
                snp->ldSum = ZPZ[i].sum();
            }
//        } else {
//            LOGGER << "Resizing LD matrix based on a chisq threshold of " << chisqThreshold << "..." << endl;
//            for (unsigned i=0; i<numIncdSnps; ++i) {
//                SnpInfo *snp = incdSnpInfoVec[i];
//                unsigned windEndi = windSizeOri[i];
//                for (unsigned j=0; j<windSizeOri[i]; ++j) {
//                    if (ZPZ[i][j]*ZPZ[i][j]*snp->sampleSize > chisqThreshold) {
//                        windStart[i] = snp->windStart = windStartOri[i] + j;
//                        break;
//                    }
//                }
//                for (unsigned j=windSizeOri[i]; j>0; --j) {
//                    if (ZPZ[i][j]*ZPZ[i][j]*snp->sampleSize > chisqThreshold) {
//                        windEndi = j;
//                        break;
//                    }
//                }
//                windSize[i] = snp->windSize = windStartOri[i] + windEndi - windStart[i];
//                ZPZiTmp = ZPZ[i].segment(windStart[i] - windStartOri[i], windSize[i]);
//                ZPZ[i] = ZPZiTmp;
//                snp->ldSamplVar = (1.0 - ZPZ[i].array().square()).square().sum()/snp->sampleSize;
//                snp->ldSum = ZPZ[i].sum();
//            }
        }
    }
    if (LDmatType == "shrunk") {
        VectorXi windStartOri = windStart;
        VectorXi windSizeOri = windSize;
        LOGGER << "Resizing LD matrix using shrunk matrix properties ..." << endl;
        // ----------------------------------------------------
        // NEW - Calculate the mutation rate components 
        // ----------------------------------------------------
        // m is the number of individuals in the reference panel for each variant
        // Need to compute theta, which is related to the mutation rate
        VectorXd nmsumi;
        VectorXd thetai;
        VectorXd mi;
        VectorXd gmapi;
        VectorXd sdss;
        //LOGGER << "Snps size " << incdSnpInfoVec.size() << endl;
        nmsumi.resize(numIncdSnps);
        thetai.resize(numIncdSnps);
        sdss.resize(numIncdSnps);
        // mi.resize(numIncdSnps);
        LOGGER << "\nUsing genetic map sample size of " << genMapN << " please alter with --genmap-n if inappropriate." << endl;
        double m = genMapN;
        gmapi.resize(numIncdSnps);
        for (unsigned i=0; i<numIncdSnps; ++i) {
            SnpInfo *snp = incdSnpInfoVec[i];
            //mi[i] = (snp->sampleSize);
            int  n = 2 * m - 1;
//            LOGGER << "snp " << i << " sample size " << snp->sampleSize << endl;
            // Approximation to the harmonic series
            nmsumi[i] = log(n) + 0.5772156649 + 1.0/ (2.0 * n) - 1.0 / (12.0 * pow(n, 2)) + 1.0 / (120.0 * pow(n, 4));
            //LOGGER << nmsumi[i] << endl;
            // Calculate theta
            //thetai[i] = (1.0 / nmsumi[i]) / (2.0 * (snp->sampleSize) + 1.0 / nmsumi[i]);
            thetai[i] = (1.0 / nmsumi[i]) / (2.0 * m + 1.0 / nmsumi[i]);
            //LOGGER <<  thetai[i] << endl;
            // Pull out the standard deviation for each variant
            sdss[i] = sqrt(2.0 * (snp->af) * (1.0 - (snp->af)));
            // LOGGER << "snp " << i << " af " << snp->af << endl;
            // 
            gmapi[i] = snp->genPos;
        }
        long int nmsum;
        double theta;
        // LOGGER << "I made it here  for shrunk sparse " << endl;
        
        // Rescale Z to the covariance scale
        VectorXd sdssSub;
        for (unsigned i=0; i<numIncdSnps; ++i) {
            sdssSub = sdss.segment(windStart[i], windSizeOri[i]);
            ZPZ[i] = (sdss[i] / 2.0) *  sdssSub.array() * ZPZ[i].array();
        }
        // LOGGER << "I made it here  for shrunk sparse " << endl;
        // --------------------------------------------------------
        // Compute the shrinkage value and then shrink the elements
        // --------------------------------------------------------
        // LOGGER << "Made it here 2 " << endl;
        double mapdiffi; 
        double rho;
        double shrinkage;
        double Ne = effpopNE;
        LOGGER << "Using European effective population size Ne=" << Ne << " please alter with --ne if inappropriate. " << endl;
        double cutoff = cutOff;
        for (unsigned i=0; i<numIncdSnps; ++i) {
            if (!(i%1000)) LOGGER << i << " SNPs processed\r";
            SnpInfo *snp = incdSnpInfoVec[i];
            for (unsigned j=i; j<=snp->windEnd; ++j) {
                mapdiffi = abs(gmapi[j] - gmapi[i]);
                rho = 4 * Ne * (mapdiffi / 100);
                shrinkage = exp(-rho / (2 * m));
                if (shrinkage <= cutoff)
                {
                    shrinkage = 0.0;
                }
                // // Multiple each covariance matrix element with the shrinkage value
                double value = ZPZ[i][j];
                if (i != (windStart[i] + j)) value = value * shrinkage;
                // // Complete as SigHAat from Li and Stephens 2003
                value =  value * ((1.0 - thetai[i]) * (1.0 - thetai[(windStart[i] + j)]));
                // If it's the diagonal element add the extra term
                if (i == (windStart[i] + j))
                {
                    value = value + 0.5f * thetai[i] * (1.0 - 0.5f * thetai[i]);
                }  
                ZPZ[j][i] = ZPZ[i][j] = value;
            }
            if(!(i%1000)) LOGGER << " Completed snp " << i << "\r" << flush;
        }
        // Now back to correlation
        for (unsigned i=0; i<numIncdSnps; ++i) {
            sdssSub = sdss.segment(windStart[i], windSizeOri[i]);
            ZPZ[i]  = (2.0 / sdss[i]) *  (1.0 / sdssSub.array()) * ZPZ[i].array();
        }
    }
    if (LDmatType == "sparseshrunk") {
        LOGGER << "Resizing sparse LD matrix using shrunk matrix properties ..." << endl;
        // ----------------------------------------------------
        // NEW - Calculate the mutation rate components 
        // ----------------------------------------------------
        // m is the number of individuals in the reference panel for each variant
        // Need to compute theta, which is related to the mutation rate
        VectorXd nmsumi;
        VectorXd thetai;
        VectorXd mi;
        VectorXd gmapi;
        VectorXd sdss;
        //LOGGER << "Snps size " << incdSnpInfoVec.size() << endl;
        nmsumi.resize(numIncdSnps);
        thetai.resize(numIncdSnps);
        sdss.resize(numIncdSnps);
        // mi.resize(numIncdSnps);
        LOGGER << "\nUsing genetic map sample size of " << genMapN << " please alter with --genmap-n if inappropriate." << endl;
        double m = genMapN;
        gmapi.resize(numIncdSnps);
        for (unsigned i=0; i<numIncdSnps; ++i) {
            SnpInfo *snp = incdSnpInfoVec[i];
            // mi[i] = (snp->sampleSize);
            int  n = 2.0 * m - 1.0;
            // LOGGER << "snp " << i << " sample size " << snp->sampleSize << endl;
            // Approximation to the harmonic series
            nmsumi[i] = log(n) + 0.5772156649 + 1.0 / (2.0 * n) - 1.0 / (12.0 * pow(n, 2)) + 1.0 / (120.0 * pow(n, 4));
            //LOGGER << nmsumi[i] << endl;
            // Calculate theta
//            thetai[i] = (1.0 / nmsumi[i]) / (2.0 * (snp->sampleSize) + 1 / nmsumi[i]);
            thetai[i] = (1.0 / nmsumi[i]) / (2.0 * m + 1 / nmsumi[i]);
            //LOGGER <<  thetai[i] << endl;
            // Pull out the standard deviation for each variant
            sdss[i] = sqrt(2.0 * (snp->af) * (1.0 - (snp->af)));
            //
            gmapi[i] = snp->genPos;
        }
        long int nmsum;
        double theta;
        // Rescale Z to the covariance scale
        double rsq = 0.0; 
        for (unsigned i=0; i<numIncdSnps; ++i) {
            // LOGGER << "Before " << ZPZsp[i] << endl;
            for (SparseVector <double> ::InnerIterator it(ZPZsp[i]); it; ++it) {
                // snpj = incdSnpInfoVec[it.index()];
                rsq = (sdss[i] / 2.0) *  sdss[it.index()] * it.value();
                // LOGGER << "sdss[i] " << sdss[i] << " sdss[it.index()] " << sdss[it.index()] << " it.value() " << it.value() << " it.index() " << it.index() << endl;
                // LOGGER << "rsq " << rsq << endl;
                it.valueRef() = rsq;
            }
            // LOGGER << "After " << ZPZsp[i] << endl;
        }
        // --------------------------------------------------------
        // Compute the shrinkage value and then shrink the elements
        // --------------------------------------------------------
        double mapdiffi; 
        double rho;
        double shrinkage;
        double Ne = effpopNE;
        LOGGER << "Using European effective population size Ne=" << Ne << " please alter with --ne if inappropriate. " << endl;
        double cutoff = cutOff;
        for (unsigned i=0; i<numIncdSnps; ++i) {
            // -----------------------------
            // Shrinkage using sparse matrix
            // -----------------------------
            double ZPZij = 0.0;
            for (SparseVector <double> ::InnerIterator it(ZPZsp[i]); it; ++it) {
                // LOGGER << " j " << j << " (windStart[i] + j) " << (windStart[i] + j) << endl;
                mapdiffi = abs(gmapi[it.index()] - gmapi[i]);
                rho = 4.0 * Ne * (mapdiffi / 100.0);
                shrinkage = exp(-rho / (2 * m)); 
                // LOGGER << "Shrinkage " << shrinkage << endl;
                if (shrinkage <= cutoff)
                {
                    shrinkage = 0.0;
                }
                // Multiple each covariance matrix element with the shrinkage value
                ZPZij = it.value() * shrinkage;
                // Complete as SigHAat from Li and Stephens 2003
                ZPZij =  ZPZij * ((1.0 - thetai[i]) * (1.0 - thetai[it.index()]));
                // If it's the diagonal element add the extra term
                if (i == (it.index()))
                {
                    ZPZij = ZPZij + 0.5f * thetai[i] * (1.0 - 0.5f * thetai[i]);
                } 
                it.valueRef()= ZPZij; 
            }
            // LOGGER << "After " << ZPZsp[i] << endl;
            if(!(i%1000)) LOGGER << " Completed snp " << i << "\r" << flush;
        }
        // // Now back to correlation
        for (unsigned i=0; i<numIncdSnps; ++i) {
            for (SparseVector <double> ::InnerIterator it(ZPZsp[i]); it; ++it) {
                // snpj = incdSnpInfoVec[it.index()];
                rsq = (2.0 / sdss[i]) *  (1.0 / sdss[it.index()]) * it.value();
                // LOGGER << "sdss[i] " << sdss[i] << " sdss[it.index()] " << sdss[it.index()] << " it.value() " << it.value() << " it.index() " << it.index() << endl;
                // LOGGER << "rsq " << rsq << endl;
                it.valueRef() = rsq;
            }
        }
    }
    displayAverageWindowSize(windSize);
}


// =============================================================================================
// Make shrunk matrix start
// =============================================================================================

// =============================================================================================
// Function read the genetic map and pass matched SNPs to main SNP inclusion exclusion tools
// =============================================================================================

void Data::readGeneticMapFile(const string &geneticMapFile){
    ifstream in(geneticMapFile.c_str());
    if (!in) LOGGER.e(0,"can not open the file [" + geneticMapFile + "] to read.");
    LOGGER << "Reading genetic map info from [" + geneticMapFile + "]." << endl;
    
    SnpInfo *snp;
    map<string, SnpInfo*>::iterator it;
    string id, gmGenPos, gmPPos;
    unsigned line=0, match=0;
    unsigned incon=0;
    while (in >> id >> gmPPos >> gmGenPos) {
        ++line;
        it = snpInfoMap.find(id);
        if (it == snpInfoMap.end()) continue;
        snp = it->second;
        if (!snp->included) continue;
//        snp->gen_map_ppos = std::stod(gmPPos.c_str());
        snp->genPos = std::stod(gmGenPos.c_str());
        ++match;
    }
    in.close();
    
    for (unsigned i=0; i<numSnps; ++i) {
        snp = snpInfoVec[i];
        if (!snp->included) continue;
        if (snp->genPos == -999 || snp->genPos == 0) {
            //LOGGER << "Who went false snp " << i << endl;
            snp->included = false;
        }
    }
    LOGGER << match << " matched SNPs in the genetic map file (in total " << line << " SNPs)." << endl;
}

void Data::readfreqFile(const string &freqFile){
    ifstream in(freqFile.c_str());
    if (!in) LOGGER.e(0,"can not open the file [" + freqFile + "] to read.");
    LOGGER << "Reading allele frequency file from [" + freqFile + "]." << endl;
    
    SnpInfo *snp;
    map<string, SnpInfo*>::iterator it;
    string id, A1, freq;
    unsigned line=0, match=0;
    unsigned incon=0;
    while (in >> id >> A1 >> freq) {
        ++line;
        it = snpInfoMap.find(id);
        if (it == snpInfoMap.end()) continue;
        snp = it->second;
        if (!snp->included) continue;
        snp->af = std::stod(freq.c_str());
        ++match;
    }
    in.close();
    
    for (unsigned i=0; i<numSnps; ++i) {
        snp = snpInfoVec[i];
        if (!snp->included) continue;
        if (snp->af == -1) {
            //LOGGER << "Who went false snp " << i << endl;
            snp->included = false;
        }
    }
    LOGGER << match << " matched SNPs in the allele frequency file (in total " << line << " SNPs)." << endl;
}


// =============================================================================================
// Function to build the shrunk matrix 
// =============================================================================================

void Data::makeshrunkLDmatrix(const string &bedFile, const string &LDmatType, const string &snpRange, const string &filename, const bool writeLdmTxt, const double effpopNE, const double cutOff, const double genMapN){
    

    Gadget::Tokenizer token;
    token.getTokens(snpRange, "-");
    
    unsigned start = 0;
    unsigned end = numIncdSnps;
    
    if (token.size()) {
        start = atoi(token[0].c_str()) - 1;
        end = atoi(token[1].c_str());
        if (end > numIncdSnps) end = numIncdSnps;
    }

    
    unsigned numSnpInRange = end - start;
    
    if (snpRange.empty())
        LOGGER << "Building shrunk LD matrix for all SNPs ..." << endl;
    else
        LOGGER << "Building shrunk LD matrix for SNPs " << snpRange << " ..." << endl;
    
    if (numIncdSnps == 0) LOGGER.e(0,"No SNP is retained for analysis.");
    if (numKeptInds == 0) LOGGER.e(0,"No individual is retained for analysis.");
    if (start >= numIncdSnps) LOGGER.e(0,"Specified a SNP range of " + snpRange + " but " + to_string(static_cast<long long>(numIncdSnps)) + " SNPs are included.");
    
    Gadget::Timer timer;
    timer.setTime();
    
    unsigned firstWindStart = 0;
    unsigned lastWindEnd = 0;
    
    // first read in the genotypes of SNPs in the given range
    
    const int bedToGeno[4] = {2, -9, 1, 0};
    unsigned size = (numInds+3)>>2;
    
    MatrixXd ZP(numSnpInRange, numKeptInds);  // SNP x Ind
    D.setZero(numSnpInRange);

    if (numKeptInds < 2) LOGGER.e(0," Cannot calculate LD matrix with number of individuals < 2.");

    FILE *in1 = fopen(bedFile.c_str(), "rb");
    if (!in1) LOGGER.e(0,"can not open the file [" + bedFile + "] to read.");
    LOGGER << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
    char header[3];
    fread(header, sizeof(header), 1, in1);
    if (!in1 || header[0] != 0x6c || header[1] != 0x1b || header[2] != 0x01) {
        cerr << "Error: Incorrect first three bytes of bed file: " << bedFile << endl;
        exit(1);
    }

    
    IndInfo *indi = NULL;
    SnpInfo *snpj = NULL;
    SnpInfo *snpk = NULL;
    
    int genoValue;
    unsigned i, j, k;
    unsigned incj, inck; // index of included SNP
    unsigned long long skipj = 0;
    unsigned nmiss;
    double mean;
    
    for (j = 0, incj = 0; j < numSnps; j++) {
        snpj = snpInfoVec[j];
        // LOGGER << "Genetic map position " << snpj->gen_map_pos << endl;
        if (snpj->index < start || !snpj->included) {
            skipj += size;
            continue;
        }
        
        if (skipj) fseek(in1, skipj, SEEK_CUR);
        skipj = 0;
        
        char *bedLineIn = new char[size];
        fread(bedLineIn, sizeof(char), size, in1);
        
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
        mean /= double(snpj->sampleSize);
        if (nmiss) {
            for (i=0; i<numKeptInds; ++i) {
                if (ZP(incj, i) == -9) ZP(incj, i) = mean;
            }
        }
        
        // compute allele frequency
        snpj->af = 0.5f*mean;
        snp2pq[incj] = 2.0f*snpj->af*(1.0-snpj->af);
        
        if (snp2pq[incj]==0) LOGGER.e(0,"" + snpj->rsID + " is a fixed SNP (MAF=0)!");
        
        // standardize genotypes
        //D[incj] = snp2pq[incj]*snpj->sampleSize;
        D[incj] = Gadget::calcVariance(ZP.row(incj))*numKeptInds;
        
        ZP.row(incj) = (ZP.row(incj).array() - mean)/sqrt(D[incj]);
        
        if (++incj == numSnpInRange) break;
    }
    
    fclose(in1);
    

    ZPZdiag = ZP.rowwise().squaredNorm(); 
    
    // then read in the bed file again to compute Z'Z
  
    MatrixXd denseZPZ;
    denseZPZ.setZero(numSnpInRange, numIncdSnps);
    VectorXd Zk(numKeptInds);
    D.setZero(numIncdSnps);

    FILE *in2 = fopen(bedFile.c_str(), "rb");
    fseek(in2, 3, SEEK_SET);
    unsigned long long skipk = 0;

    for (k = 0, inck = 0; k < numSnps; k++) {
        snpk = snpInfoVec[k];
        
        if (!snpk->included) {
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
        mean /= double(snpk->sampleSize);
        if (nmiss) {
            for (i=0; i<numKeptInds; ++i) {
                if (Zk[i] == -9) Zk[i] = mean;
            }
        }
        
        // compute allele frequency
        snpk->af = 0.5f*mean;
        snp2pq[inck] = 2.0f*snpk->af*(1.0-snpk->af);
        
        if (snp2pq[inck]==0) LOGGER.e(0,"" + snpk->rsID + " is a fixed SNP (MAF=0)!");
        
        // standardize genotypes
        //D[inck] = snp2pq[inck]*snpk->sampleSize;
        D[inck] = Gadget::calcVariance(Zk)*numKeptInds;
        
        Zk = (Zk.array() - mean)/sqrt(D[inck]);
        
        denseZPZ.col(inck) = ZP * Zk;
        
        if(!(inck%1000)) LOGGER << " read snp " << inck << "\r" << flush;
        
        ++inck;
    }

    fclose(in2);
    
    // ----------------------------------------------------
    // NEW - Calculate the mutation rate components 
    // ----------------------------------------------------
    // m is the number of individuals in the reference panel for each variant
    // Need to compute theta, which is related to the mutation rate
    VectorXd nmsumi(numIncdSnps);
    VectorXd thetai(numIncdSnps);
    VectorXd mi(numIncdSnps);
    VectorXd gmapi(numIncdSnps);
    VectorXd sdss(numIncdSnps);
    LOGGER << "\nUsing genetic map sample size of " << genMapN << " please alter with --genmap-n if inappropriate." << endl;
    double m = genMapN;
    for (unsigned i=0; i<numIncdSnps; ++i) {
            SnpInfo *snp = incdSnpInfoVec[i];
            //mi[i] = (snp->sampleSize);
            int  n = 2.0 * m - 1.0;
            // Approximation to the harmonic series
            nmsumi[i] = log(n) + 0.5772156649 + 1.0 / (2.0 * n) - 1.0 / (12.0 * pow(n, 2.0)) + 1.0 / (120.0 * pow(n, 4.0));
            //LOGGER << nmsumi[i] << endl;
            // Calculate theta
            //thetai[i] = (1.0 / nmsumi[i]) / (2.0 * (snp->sampleSize) + 1.0 / nmsumi[i]);
            thetai[i] = (1.0 / nmsumi[i]) / (2.0 * m + 1.0 / nmsumi[i]);
            //LOGGER <<  thetai[i] << endl;
            // Pull out the standard deviation for each variant
            sdss[i] = sqrt(2.0 * (snp->af) * (1.0 - (snp->af)));
            // 
            gmapi[i] = snp->genPos;
    }
    long int nmsum;
    double theta;
     

    // Rescale Z to the covariance scale
    for (unsigned i=0; i<numIncdSnps; ++i) {
        denseZPZ.col(i) = (sdss[i] / 2.0) *  sdss.array() * denseZPZ.col(i).array();
    }
    // --------------------------------------------------------
    // Compute the shrinkage value and then shrink the elements
    // --------------------------------------------------------
    double mapdiffi; 
    double rho;
    double shrinkage;
    double Ne = effpopNE;
    LOGGER << "\nUsing European effective population size Ne=" << Ne << " please alter with --ne if inappropriate." << endl;
    double cutoff = cutOff;
    for (unsigned i=0; i<numSnpInRange; ++i) {
        for (unsigned j=0; j<numIncdSnps; ++j) {
            mapdiffi = abs(gmapi[j] - gmapi[start+i]);
            rho = 4.0 * Ne * (mapdiffi / 100.0);
            shrinkage = exp(-rho / (2 * m)); 
            // if (i <=10 && j <= 10)
            // { 
            //     LOGGER << "Snp " << i << " " << j << " Mapdiff " << mapdiffi << " shrinkage " << shrinkage << endl;
            //     LOGGER << "Start position " << start  << " start plus i " << start + i << " gmap start plus i " << gmapi[start +i] << endl;
            //     // LOGGER << gmapi[i] << " " << gmapi[j] << endl;
            //     // LOGGER << mapdiffi << endl;
            //     // LOGGER << shrinkage << endl;
            //     // LOGGER << mi[i] << endl;
            // }
            if (shrinkage <= cutoff) {
                shrinkage = 0.0;
            }
            // Multiple each covariance matrix element with the shrinkage value
            denseZPZ(i, j) *= shrinkage;
            // Complete as SigHAat from Li and Stephens 2003
            denseZPZ(i, j) *= (1.0 - thetai[start+i]) * (1.0 - thetai[j]);
            // If it's the diagonal element add the extra term
            if (start+i == j)
            {
               denseZPZ(i, j) = denseZPZ(i, j) + 0.5f * thetai[start+i] * (1 - 0.5f * thetai[start+i]);
            }
            // Make the upper triangle equal to the lower triangle
            //denseZPZ(j, i) =  denseZPZ(i, j);
        }
    }
    // Now back to correlation
    for (unsigned i=0; i<numIncdSnps; ++i) {
        denseZPZ.col(i) = (2.0 / sdss[i]) *  (1.0 / sdss.array()) * denseZPZ.col(i).array();
    }
    // --------------------------------
    // Copy to vector of vectors format
    // --------------------------------
    ZPZ.resize(numSnpInRange);
    windStart.setZero(numSnpInRange);
    windSize.setConstant(numSnpInRange, numIncdSnps);
    // LOGGER << "Num snp range " << numSnpInRange << endl;
    for (unsigned i=0; i<numSnpInRange; ++i) {
        SnpInfo *snp = incdSnpInfoVec[start+i];
        snp->windStart = 0;
        snp->windSize  = numIncdSnps;
        snp->windEnd   = numIncdSnps-1;
        ZPZ[i] = denseZPZ.row(i);
        snp->ldSum = denseZPZ.row(i).sum();
    }
        // LOGGER << ZPZ[0] << endl;
    //}
    // else if (LDmatType == "band") {
    //     ZPZ.resize(numSnpInRange);
    //     if (windowWidth) {  // based on the given window width
    //         for (unsigned i=0; i<numSnpInRange; ++i) {
    //             SnpInfo *snp = incdSnpInfoVec[start+i];
    //             ZPZ[i] = denseZPZ.row(i).segment(snp->windStart, snp->windSize);
    //         }
    //     } else {  // based on the given LD threshold
    //         windStart.setZero(numSnpInRange);
    //         windSize.setZero(numSnpInRange);
    //         for (unsigned i=0; i<numSnpInRange; ++i) {
    //             SnpInfo *snp = incdSnpInfoVec[start+i];
    //             unsigned windEndi = numIncdSnps;
    //             for (unsigned j=0; j<numIncdSnps; ++j) {
    //                 if (abs(denseZPZ(i,j)) > LDthreshold) {
    //                     windStart[i] = snp->windStart = j;
    //                     break;
    //                 }
    //             }
    //             for (unsigned j=numIncdSnps; j>0; --j) {
    //                 if (abs(denseZPZ(i,j-1)) > LDthreshold) {
    //                     windEndi = j;
    //                     break;
    //                 }
    //             }
    //             windSize[i] = snp->windSize = windEndi - windStart[i];
    //             snp->windEnd = windEndi - 1;
    //             ZPZ[i].resize(windSize[i]);
    //             VectorXd::Map(&ZPZ[i][0], windSize[i]) = denseZPZ.row(i).segment(windStart[i], windSize[i]);
    //         }
    //     }
    // }
    // else if (LDmatType == "sparse") {
    //     ZPZsp.resize(numSnpInRange);
    //     windStart.setZero(numSnpInRange);
    //     windSize.setZero(numSnpInRange);
    //     if (LDthreshold) {
    //         for (unsigned i=0; i<numSnpInRange; ++i) {
    //             SnpInfo *snp = incdSnpInfoVec[start+i];
    //             for (unsigned j=0; j<numIncdSnps; ++j) {
    //                 if (abs(denseZPZ(i,j)) < LDthreshold) denseZPZ(i,j) = 0;
    //             }
    //             ZPZsp[i] = denseZPZ.row(i).sparseView();
    //             SparseVector <double> ::InnerIterator it(ZPZsp[i]);
    //             windStart[i] = snp->windStart = it.index();
    //             windSize[i] = snp->windSize = ZPZsp[i].nonZeros();
    //             for (; it; ++it) snp->windEnd = it.index();
    //         }
    //     } else {
    //         for (unsigned i=0; i<numSnpInRange; ++i) {
    //             SnpInfo *snp = incdSnpInfoVec[start+i];
    //             for (unsigned j=0; j<numIncdSnps; ++j) {
    //                 if (denseZPZ(i,j)*denseZPZ(i,j)*snp->sampleSize < chisqThreshold) denseZPZ(i,j) = 0;
    //             }
    //             ZPZsp[i] = denseZPZ.row(i).sparseView();
    //             SparseVector <double> ::InnerIterator it(ZPZsp[i]);
    //             windStart[i] = snp->windStart = it.index();
    //             windSize[i] = snp->windSize = ZPZsp[i].nonZeros();
    //             for (; it; ++it) snp->windEnd = it.index();
    //         }
    //     }
    // }
    
    denseZPZ.resize(0,0);
    
    //    LOGGER << denseZPZ.block(0,0,10,10) << endl;
    //    LOGGER << windStart.transpose() << endl;
    //    LOGGER << windSize.transpose() << endl;
    
    
    timer.getTime();
    
    
    LOGGER << endl;
    displayAverageWindowSize(windSize);
    LOGGER << "LD matrix diagonal mean " << ZPZdiag.mean() << " variance " << Gadget::calcVariance(ZPZdiag) << "." << endl;
    if (ZPZdiag.mean() < 0.8 || ZPZdiag.mean() > 1.2) LOGGER << "ERROR: The mean of LD matrix diagonal values is expected to be close to one. Something is wrong with the LD matrix!" << endl;
LOGGER << "Genotype data for " << numKeptInds << " individuals and " << numSnpInRange << " SNPs are included from [" + bedFile + "]." << endl;
    LOGGER << "Build of LD matrix completed (time used: " << timer.format(timer.getElapse()) << ")." << endl;
    
    
    vector<SnpInfo*> snpVecTmp(numSnpInRange);
    for (unsigned i=0; i<numSnpInRange; ++i) {
        snpVecTmp[i] = incdSnpInfoVec[start+i];
    }
    incdSnpInfoVec = snpVecTmp;
    numIncdSnps = numSnpInRange;
    string outfilename = filename;
    if (!snpRange.empty()) outfilename += ".snp" + snpRange;
    outputLDmatrix(LDmatType, outfilename, writeLdmTxt);
}

// =============================================================================================
// Make shrunk matrix end
// =============================================================================================

void Data::buildSparseMME(const bool sampleOverlap, const bool noscale){
    VectorXd Dref = snp2pq*numKeptInds;
    snp2pq.resize(numIncdSnps);
    D.resize(numIncdSnps);
//    ZPZdiag.resize(numIncdSnps);
    ZPy.resize(numIncdSnps);
    b.resize(numIncdSnps);
    n.resize(numIncdSnps);
    se.resize(numIncdSnps);
    tss.resize(numIncdSnps);
    SnpInfo *snp;
    for (unsigned i=0; i<numIncdSnps; ++i) {
        snp = incdSnpInfoVec[i];
        snp->af = snp->gwas_af;
        snp2pq[i] = snp->twopq = 2.0f*snp->gwas_af*(1.0-snp->gwas_af);
        if(snp2pq[i]==0) LOGGER << "Error: SNP " << snp->rsID << " af " << snp->af << " has 2pq = 0." << endl;
        D[i] = snp2pq[i]*snp->gwas_n;
        b[i] = snp->gwas_b;
        n[i] = snp->gwas_n;
        se[i]= snp->gwas_se;
        tss[i] = D[i]*(n[i]*se[i]*se[i] + b[i]*b[i]);
        // cout << snp->rsID << " b: " << b[i]  << " se: "  << se[i] << endl;
//        D[i] = 1.0/(se[i]*se[i]+b[i]*b[i]/snp->gwas_n);  // NEW!
//        snp2pq[i] = snp->twopq = D[i]/snp->gwas_n;       // NEW!
    }
    //b.array() -= b.mean();  // DO NOT CENTER b

    // estimate phenotypic variance based on the input allele frequencies in GWAS
    //ypy = (D.array()*(n.array()*se.array().square()+b.array().square())).mean();
    VectorXd ypySrt = D.array()*(n.array()*se.array().square()+b.array().square());
    VectorXd varpSrt = ypySrt.array()/n.array();
    std::sort(ypySrt.data(), ypySrt.data() + ypySrt.size());
    std::sort(varpSrt.data(), varpSrt.data() + varpSrt.size());
    ypy = ypySrt[ypySrt.size()/2];  // median
    varPhenotypic = varpSrt[varpSrt.size()/2];

    // cout << "gwas varPhenotypic : " <<  varPhenotypic << endl;

    //numKeptInds = n.mean();

    VectorXd nSrt = n;
    std::sort(nSrt.data(), nSrt.data() + nSrt.size());
    numKeptInds = nSrt[nSrt.size()/2]; // median

        // NEW
        // compute D and snp2pq based on n, se and b, assuming varp = 1
        // these quantities are used in sbayes, as they are more reliable than input allele frequencies
        for (unsigned i=0; i<numIncdSnps; ++i) {
            snp = incdSnpInfoVec[i];
            D[i] = varPhenotypic/(se[i]*se[i]+b[i]*b[i]/snp->gwas_n);  // NEW!
            snp2pq[i] = snp->twopq = D[i]/snp->gwas_n;       // NEW!
            tss[i] = D[i]*(n[i]*se[i]*se[i] + b[i]*b[i]);
            // Need to adjust R and C models X'X matrix depending scale of genotypes or not
            if (noscale == true) {
                D[i] = snp2pq[i]*snp->gwas_n;
            } else {
                D[i] = snp->gwas_n;
            }
        }
        //ypy = numKeptInds;
        // NEW END

    if (ZPZ.size() || ZPZsp.size()) {
        if (sparseLDM == true) {
            for (unsigned i=0; i<numIncdSnps; ++i) {
                snp = incdSnpInfoVec[i];
                //LOGGER << i << " " << ZPZsp[i].nonZeros() << " " << D.size() << endl;
                for (SparseVector <double> ::InnerIterator it(ZPZsp[i]); it; ++it) {
                    //LOGGER << it.index() << " ";
                    it.valueRef() *= sqrt(D[i]*D[it.index()]);
                }
            }
        } else {
            for (unsigned i=0; i<numIncdSnps; ++i) {
                snp = incdSnpInfoVec[i];
                for (unsigned j=0; j<snp->windSize; ++j) {
                    ZPZ[i][j] *= sqrt(D[i]*D[snp->windStart+j]);
                }
            }
        }

        // sum of sampling variance of LD for each SNP with all other SNPs
        // for significant LD, the sampling variance is proportional to the (ratio of ref and gwas n) + 1
        // for insignificant LD, the sampling variance is 1 over gwas n
        LDsamplVar.resize(numIncdSnps);
        LDscore.resize(numIncdSnps);
        for (unsigned i=0; i<numIncdSnps; ++i) {
            snp = incdSnpInfoVec[i];
            if (sampleOverlap) {
                LDsamplVar[i] = 0;
            } else {
                LDsamplVar[i]  = (snp->gwas_n + snp->sampleSize)/double(numIncdSnps)*snp->ldSamplVar;
            }
            LDsamplVar[i] += (numIncdSnps - snp->numNonZeroLD)/double(numIncdSnps);
            //if (sampleOverlap) LDsamplVar[i] = 0;
            LDscore[i] = snp->ldsc; //*snp->gwas_n;
        }

        ZPZdiag.array() *= D.array();
    }
    else {
        Dratio = D.array()/Dref.array();
        ZPZdiag.array() *= Dratio.array();
        DratioSqrt = Dratio.array().sqrt();
        for (unsigned i=0; i<numIncdSnps; ++i) {
            Z.col(i) *= DratioSqrt[i];
        }
    }


    if (noscale) {
        ZPy = ZPZdiag.cwiseProduct(b);
    } else {
        ZPy = ZPZdiag.cwiseProduct(b).cwiseProduct(snp2pq.array().sqrt().matrix());
    }
    chisq = ZPy.cwiseProduct(b);

    numKeptInds = n.mean();


    numFixedEffects = 1;
    fixedEffectNames.resize(1);
    fixedEffectNames[0] = "Intercept";
    ZPX.resize(numIncdSnps,1);
    if (ZPZ.size()) {
        for (unsigned i=0; i<numIncdSnps; ++i) {
            ZPX(i,0) = ZPZ[i].sum();
        }
    } else if (ZPZsp.size()) {
        for (unsigned i=0; i<numIncdSnps; ++i) {
            ZPX(i,0) = ZPZsp[i].sum();
        }
    } else {
        LOGGER.e(0," either dense or sparse ldm does not exist!");
    }
    XPX.resize(1,1);
    XPX << ZPX.col(0).sum();
    XPXdiag.resize(1);
    XPXdiag << XPX(0,0);
    XPy.resize(1,1);
    XPy << ZPy.sum();


    // data summary
    LOGGER << "\nData summary:" << endl;
    LOGGER << boost::format("%40s %8s %8s\n") %"" %"mean" %"sd";
    LOGGER << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP Phenotypic variance" %Gadget::calcMean(varpSrt) %sqrt(Gadget::calcVariance(varpSrt));
    LOGGER << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP heterozygosity" %Gadget::calcMean(snp2pq) %sqrt(Gadget::calcVariance(snp2pq));
    LOGGER << boost::format("%40s %8.0f %8.0f\n") %"GWAS SNP sample size" %Gadget::calcMean(n) %sqrt(Gadget::calcVariance(n));
    LOGGER << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP effect" %Gadget::calcMean(b) %sqrt(Gadget::calcVariance(b));
    LOGGER << boost::format("%40s %8.3f %8.3f\n") %"GWAS SNP SE" %Gadget::calcMean(se) %sqrt(Gadget::calcVariance(se));
    LOGGER << boost::format("%40s %8.3f %8.3f\n") %"MME left-hand-side diagonals" %Gadget::calcMean(ZPZdiag) %sqrt(Gadget::calcVariance(ZPZdiag));
    LOGGER << boost::format("%40s %8.3f %8.3f\n") %"MME right-hand-side" %Gadget::calcMean(ZPy) %sqrt(Gadget::calcVariance(ZPy));
    LOGGER << boost::format("%40s %8.3f %8.3f\n") %"LD sampling variance" %Gadget::calcMean(LDsamplVar) %sqrt(Gadget::calcVariance(LDsamplVar));
    LOGGER << boost::format("%40s %8.3f %8.3f\n") %"LD score" %Gadget::calcMean(LDscore) %sqrt(Gadget::calcVariance(LDscore));
}

void Data::outputSnpEffectSamples(const SpMat &snpEffects, const unsigned burnin, const unsigned outputFreq, const string&snpResFile, const string &filename) const {
    LOGGER << "writing SNP effect samples into " << filename << endl;
    unsigned nrow = snpEffects.rows();
    vector<string> snpName;
    vector <double>  sample;

    ifstream in(snpResFile.c_str());
    if (!in) LOGGER.e(0,"can not open the snpRes file [" + snpResFile + "] to read.");

    Gadget::Tokenizer colData;
    string inputStr;
    string sep(" \t");
    string id;
    unsigned line=0;
    while (getline(in,inputStr)) {
        ++line;
        if (line==1) continue;
        colData.getTokens(inputStr, sep);
        snpName.push_back(colData[1]);
    }
    in.close();
    long numSnps = snpName.size();
    
    ofstream out(filename.c_str());
    out << boost::format("%6s %20s %8s\n")
    % "Iteration"
    % "Name"
    % "Sample";
    
    LOGGER << "Size of mcmc samples " << snpEffects.rows() << " " << snpEffects.cols() << endl;
    
    unsigned idx=0;
    for (unsigned iter=0; iter<nrow; ++iter) {
        if (iter < burnin) continue;
        if (!(iter % outputFreq)) {
            ++idx;
            for (unsigned j=0; j<numSnps; ++j) {
                if (snpEffects.coeff(iter, j)) {
                    out << boost::format("%6s %20s %8s\n")
                    % idx
                    % snpName[j]
                    % snpEffects.coeff(iter, j);
                }
            }
        }
    }
    
    out.close();
}

void Data::directPruneLDmatrix(const string &ldmatrixFile, const string &outLDmatType, const double chisqThreshold, const string &title, const bool writeLdmTxt){
    readLDmatrixInfoFile(ldmatrixFile + ".info");
    Gadget::Tokenizer token;
    token.getTokens(ldmatrixFile, ".");
    string ldmType = token[token.size()-1];
    if (ldmType != "full") {
        LOGGER.e(0," --direct-prune only prunes a full LD matrix to a sparse matrix!");
    }
    LOGGER << "Direct pruning " + ldmType + " LD matrix from [" + ldmatrixFile + "]..." << endl;

    FILE *in = fopen((ldmatrixFile+".bin").c_str(), "rb");
    if (!in) {
        LOGGER.e(0," cannot open LD matrix file " + ldmatrixFile);
    }
    SnpInfo *snp;
    VectorXd vec;
    double rsq = 0.0;

    string outfilename = title + ".ldm." + outLDmatType;
    string outfile1 = outfilename + ".info";
    string outfile2 = outfilename + ".bin";
    ofstream out1(outfile1.c_str());
    FILE *out2 = fopen(outfile2.c_str(), "wb");
    ofstream out3;
    string outfile3;
    if (writeLdmTxt) {
        outfile3 = outfilename + ".txt";
        out3.open(outfile3.c_str());
    }
    out1 << boost::format("%6s %15s %10s %15s %6s %6s %12s %10s %10s %10s %10s %15s %10s %12s %12s\n")
    % "Chrom"
    % "ID"
    % "GenPos"
    % "PhysPos"
    % "A1"
    % "A2"
    % "A1Freq"
    % "Index"
    % "WindStart"
    % "WindEnd"
    % "WindSize"
    % "WindWidth"
    % "N"
    % "SamplVar"
    % "LDsum";

    SnpInfo *windStart, *windEnd;
    for (unsigned i=0; i<numSnps; ++i) {
        snp = snpInfoVec[i];
        snp->ldSamplVar = 0.0;
        snp->ldSum = 0.0;
        
        vec.resize(snp->windSize);
        fread(vec.data(), sizeof(double), snp->windSize, in);
        
//        if (!(i%100)) LOGGER << " SNP " << i << "\r";
        
        for (unsigned j=0; j<snp->windSize; ++j) {
            rsq = vec[j]*vec[j];
            if (i!=j && rsq*snp->sampleSize < chisqThreshold) {
                vec[j] = 0;
            } else {
                snp->ldSamplVar += (1.0-rsq)*(1.0-rsq)/snp->sampleSize;
                snp->ldSum += vec[j];
            }
        }
        
        SparseVector <double>  spvec = vec.sparseView();
        SparseVector <double> ::InnerIterator it(spvec);
        snp->windStart = it.index();
        snp->windSize = spvec.nonZeros();
        for (; it; ++it) snp->windEnd = it.index();
        
        out1 << boost::format("%6s %15s %10s %15s %6s %6s %12f %10s %10s %10s %10s %15s %10s %12.6f %12.6f\n")
        % snp->chrom
        % snp->rsID
        % snp->genPos
        % snp->physPos
        % snp->a1
        % snp->a2
        % snp->af
        % snp->index
        % snp->windStart
        % snp->windEnd
        % snp->windSize
        % -9
        % snp->sampleSize
        % snp->ldSamplVar
        % snp->ldSum;
        fwrite(spvec.innerIndexPtr(), sizeof(unsigned), snp->windSize, out2);
        fwrite(spvec.valuePtr(), sizeof(double), snp->windSize, out2);
        if (writeLdmTxt) out3 << vec.transpose() << endl;
    }
    out1.close();
    fclose(out2);
    
    LOGGER << "Written the SNP info into file [" << outfile1 << "]." << endl;
    LOGGER << "Written the LD matrix into file [" << outfile2 << "]." << endl;
    
    if (writeLdmTxt) {
        out3.close();
        LOGGER << "Written the LD matrix into text file [" << outfile3 << "]." << endl;
    }
}

void Data::addLDmatrixInfo(const string &ldmatrixFile) {
    ifstream in(ldmatrixFile.c_str());
    if (!in) LOGGER.e(0,"can not open the file [" + ldmatrixFile + "] to read.");
    LOGGER << "Reading SNP info from [" + ldmatrixFile + "]." << endl;
    string header;
    string id, allele1, allele2;
    unsigned chr, physPos;
    double genPos, af, ldSamplVar, ldSum;
    unsigned idx, windStart, windEnd, windSize, windWidth;
    long sampleSize;
    bool skeleton;
    getline(in, header);
    Gadget::Tokenizer token;
    token.getTokens(header, " ");
    map<string, SnpInfo*>::iterator it, end = snpInfoMap.end();
    unsigned count = 0;
    if (token.back() == "Skeleton") {
        while (in >> chr >> id >> genPos >> physPos >> allele1 >> allele2 >> af >> idx >> windStart >> windEnd >> windSize >> windWidth >> sampleSize >> ldSamplVar >> ldSum >> skeleton) {
            it = snpInfoMap.find(id);
            if (it == end) LOGGER.e(0,"SNP " + id + " not found in .bim file.");
            SnpInfo *snp = it->second;
            if (snp->index != idx) LOGGER.e(0,"The index of SNP " + id + " is different from that in .bim file.");
            snp->windStart = windStart;
            snp->windEnd = windEnd;
            snp->windSize = windSize;
            snp->sampleSize = sampleSize;
            snp->ldSamplVar = ldSamplVar;
            snp->ldSum = ldSum;
            snp->skeleton = skeleton;
            ++count;
        }
    }
    else if (token.back() == "LDsum") {
        while (in >> chr >> id >> genPos >> physPos >> allele1 >> allele2 >> af >> idx >> windStart >> windEnd >> windSize >> windWidth >> sampleSize >> ldSamplVar >> ldSum) {
            it = snpInfoMap.find(id);
            if (it == end) LOGGER.e(0,"SNP " + id + " not found in .bim file.");
            SnpInfo *snp = it->second;
            if (snp->index != idx) LOGGER.e(0,"The index of SNP " + id + " is different from that in .bim file.");
            snp->windStart = windStart;
            snp->windEnd = windEnd;
            snp->windSize = windSize;
            snp->sampleSize = sampleSize;
            snp->ldSamplVar = ldSamplVar;
            snp->ldSum = ldSum;
            ++count;
        }
    }
    else {
        while (in >> chr >> id >> genPos >> physPos >> allele1 >> allele2 >> af >> idx >> windStart >> windEnd >> windSize >> windWidth >> sampleSize >> ldSamplVar) {
            it = snpInfoMap.find(id);
            if (it == end) LOGGER.e(0,"SNP " + id + " not found in .bim file.");
            SnpInfo *snp = it->second;
            if (snp->index != idx) LOGGER.e(0,"The index of SNP " + id + " is different from that in .bim file.");
            snp->windStart = windStart;
            snp->windEnd = windEnd;
            snp->windSize = windSize;
            snp->sampleSize = sampleSize;
            snp->ldSamplVar = ldSamplVar;
            ++count;
        }
    }
    in.close();
    LOGGER << count << " SNPs to be included from [" + ldmatrixFile + "]." << endl;

}

void Data::jackknifeLDmatrix(const string &ldmatrixFile, const string &outLDmatType, const string &title, const bool writeLdmTxt) {
    // Jackknife estimate of correlation and sampling variance
    Gadget::Tokenizer token;
    token.getTokens(ldmatrixFile, ".");
    if (token[token.size()-1] != "sparse") {
        LOGGER.e(0," --jackknife method only works for sparse LD matrix at the moment!");
    }

    addLDmatrixInfo(ldmatrixFile + ".info");

    FILE *in = fopen((ldmatrixFile + ".bin").c_str(), "rb");
    if (!in) {
        LOGGER.e(0," cannot open LD matrix file " + ldmatrixFile + ".bin");
    }

    if (numIncdSnps == 0) LOGGER.e(0,"No SNP is retained for analysis.");
    
    LOGGER << "Reading and jackknifing sparse LD matrix from [" + ldmatrixFile + ".bin]..." << endl;
    
    SnpInfo *snpi, *snpj;
    SnpInfo *windStart, *windEnd;
    SparseVector <double>  ZPZspvec;
    VectorXd ZiCwiseZj(numKeptInds);
    VectorXd ones;
    ones.setOnes(numKeptInds);
    double numJackknife = numKeptInds-1;
    double samplVarEmp;
    
    string outfilename = title + ".ldm." + outLDmatType;
    string outfile = outfilename + ".info";
    ofstream out(outfile.c_str());
    out << boost::format("%6s %15s %10s %15s %6s %6s %12s %10s %10s %10s %10s %15s %10s %12s %12s\n")
    % "Chrom"
    % "ID"
    % "GenPos"
    % "PhysPos"
    % "A1"
    % "A2"
    % "A1Freq"
    % "Index"
    % "WindStart"
    % "WindEnd"
    % "WindSize"
    % "WindWidth"
    % "N"
    % "SamplVar"
    % "LDsum";
    
    
    // standardize genotype matrix
    for (unsigned i=0; i<numIncdSnps; ++i) {
        Z.col(i) /= sqrt(numKeptInds*snp2pq[i]);
    }

    
    for (unsigned i=0, inci=0; i<numSnps; i++) {
        snpi = snpInfoVec[i];
        
        unsigned d[snpi->windSize];
        double v[snpi->windSize];
        
        if (!snpi->included) {
            fseek(in, sizeof(d), SEEK_CUR);
            fseek(in, sizeof(v), SEEK_CUR);
            continue;
        }
        
        fread(d, sizeof(d), 1, in);
        fread(v, sizeof(v), 1, in);
        
        ZPZspvec.resize(snpi->windSize);
//        snpi->ldSamplVar = 0.0;
        snpi->ldSum = 0.0;
        
        for (unsigned j=0; j<snpi->windSize; ++j) {
            snpj = snpInfoVec[d[j]];
            if (snpj->included) {
                ZiCwiseZj  = v[j]*ones - Z.col(snpi->index).cwiseProduct(Z.col(snpj->index));
                ZiCwiseZj *= numKeptInds/numJackknife;
                samplVarEmp = (ZiCwiseZj.array() - ZiCwiseZj.mean()).square().sum() * numJackknife/numKeptInds;
                snpi->ldSum += samplVarEmp;
            }
        }

        windStart = incdSnpInfoVec[snpi->windStart];
        windEnd = incdSnpInfoVec[snpi->windEnd];
        out << boost::format("%6s %15s %10s %15s %6s %6s %12f %10s %10s %10s %10s %15s %10s %12.6f %12.6f\n")
        % snpi->chrom
        % snpi->rsID
        % snpi->genPos
        % snpi->physPos
        % snpi->a1
        % snpi->a2
        % snpi->af
        % snpi->index
        % snpi->windStart
        % snpi->windEnd
        % snpi->windSize
        % (windStart->chrom == windEnd->chrom ? windEnd->physPos - windStart->physPos : windStart->chrom-windEnd->chrom)
        % snpi->sampleSize
        % snpi->ldSamplVar
        % snpi->ldSum;

        if (++inci == numIncdSnps) break;
    }

    out.close();
    LOGGER << "Written the SNP info into file [" << outfile << "]." << endl;

}

void Data::readLDscoreFile(const string &ldscoreFile) {
    ifstream in(ldscoreFile.c_str());
    if (!in) LOGGER.e(0,"can not open the LD score file [" + ldscoreFile + "] to read.");
    LOGGER << "Reading SNP LD scores from [" + ldscoreFile + "]." << endl;

    readLDscore = true;
    
    SnpInfo *snp;
    map<string, SnpInfo*>::iterator it;
    string id;
    double ldsc;
    unsigned line = 0, match = 0;
    while(in >> id >> ldsc) {
        ++line;
        it = snpInfoMap.find(id);
        if (it == snpInfoMap.end()) continue;
        snp = it->second;
        if (snp->included) {
            snp->ldsc = ldsc;
            ++match;
        }
    }
    in.close();
    
    for (unsigned i=0; i<numSnps; ++i) {
        snp = snpInfoVec[i];
        if (!snp->included) continue;
        if (!snp->ldsc) {
            snp->included = false;
        }
    }
    
    LOGGER << "Read LD scores for " << match << " matched SNPs (in total " << line << " SNPs)." << endl;
}


void Data::getZPZspmat(){ // get ZPZspmat from ZPZsp which is a vector of sparse vectors
    ZPZspmat.resize(numIncdSnps, numIncdSnps);
    vector<Triplet <double>  > tripletList;
    tripletList.reserve(windSize.cast<double>().sum());
    double val = 0.0;
    for (unsigned i=0; i<numIncdSnps; ++i) {
        SnpInfo *snpi = incdSnpInfoVec[i];
        if (!(i % 100000)) LOGGER << "  making sparse LD matrix for SNP " << i << " " << windSize[i] << " " << snpi->windSize << " " << ZPZsp[i].size() << " " << ZPZsp[i].nonZeros() << endl;
        for (SparseVector <double> ::InnerIterator it(ZPZsp[i]); it; ++it) {
            tripletList.push_back(Triplet <double> (i, it.index(), it.value()));
            //if (it.index() > numIncdSnps - 1) LOGGER << i << " " << it.index() << " " << it.value() << endl;
        }
    }
    ZPZspmat.setFromTriplets(tripletList.begin(), tripletList.end());
    ZPZspmat.makeCompressed();
    tripletList.clear();
}

void Data::getZPZmat(){ // get ZPZmat from ZPZ which is a vector of vectors
    ZPZmat.setZero(numIncdSnps, numIncdSnps);
    for (unsigned i=0; i<numIncdSnps; ++i) {
        ZPZmat.col(i) = ZPZ[i];
    }
}

void Data::binSnpByLDrsq(const double rsqThreshold, const string &title){
// only work for full LD matrix at the moment
    
    // initialise a unique window ID for each SNP
    for (unsigned i=0; i<numIncdSnps; ++i) {
        SnpInfo *snp = incdSnpInfoVec[i];
        snp->window = i;
    }
    
    // merge SNP window ID if two SNPs are in high LD
    for (unsigned i=0; i<numIncdSnps; ++i) {
        SnpInfo *snp = incdSnpInfoVec[i];
        VectorXd rsq = ZPZ[i].cwiseProduct(ZPZ[i]);
        for (unsigned j=0; j<numIncdSnps; ++j) {
            if (rsq[j] > rsqThreshold) incdSnpInfoVec[j]->window = snp->window;
        }
    }

    // reorder the window ID
    map<int, vector<SnpInfo*> > uniqueWinIDmap;
    for (unsigned i=0; i<numIncdSnps; ++i) {
        SnpInfo *snp = incdSnpInfoVec[i];
        uniqueWinIDmap[snp->window].push_back(snp);
    }
    map<int, vector<SnpInfo*> >::iterator it, end = uniqueWinIDmap.end();
    int newWinID = 1;
    for (it = uniqueWinIDmap.begin(); it!=end; ++it) {
        for (unsigned i=0; i<it->second.size(); ++i) {
            it->second[i]->window = newWinID;
        }
        ++newWinID;
    }
    
    // output SNP window ID info
    string outfile = title + ".window";
    ofstream out(outfile.c_str());
    out << boost::format("%18s %8s\n") % "SNP" % "Window";
    for (unsigned i=0; i<numIncdSnps; ++i) {
        SnpInfo *snp = incdSnpInfoVec[i];
        out << boost::format("%18s %8s\n") % snp->rsID % snp->window;
    }
    out.close();
    
    LOGGER << "Based on LD Rsq threshold of " << rsqThreshold << ", " << numIncdSnps << " SNPs are grouped into " << newWinID << " bins." << endl;
}

void Data::readWindowFile(const string &windowFile){
    ifstream in(windowFile.c_str());
    if (!in) LOGGER.e(0,"can not open the SNP window file [" + windowFile + "] to read.");
    LOGGER << "Reading SNP window info from [" + windowFile + "]." << endl;
    
    SnpInfo *snp;
    map<string, SnpInfo*>::iterator it;
    string snpID;
    int windowID;
    set<int> uniqueWinID;
    unsigned line=0, match=0;
    string header;
    getline(in, header);
    while (in >> snpID >> windowID) {
        ++line;
        it = snpInfoMap.find(snpID);
        if (it == snpInfoMap.end()) continue;
        snp = it->second;
        snp->window = windowID;
        ++match;
        uniqueWinID.insert(windowID);
    }
    
    numWindows = uniqueWinID.size();

    LOGGER << match << " matched SNPs in the SNP window file (in total " << line << " SNPs) are binned into " << numWindows << " windows." << endl;
}

void Data::binSnpByWindowID(){
    // reorder the window ID in case the original window ID are not continous
    map<int, int> orderedWinIDmap;
    for (unsigned i=0; i<numIncdSnps; ++i) {
        SnpInfo *snp = incdSnpInfoVec[i];
        orderedWinIDmap[snp->window] = 0;
    }
    map<int, int>::iterator it, end = orderedWinIDmap.end();
    int newWinID = 0;
    for (it = orderedWinIDmap.begin(); it!=end; ++it) {
        it->second = newWinID++;
    }

    // set up window-snp map
    windowSnpIdxVec.resize(numWindows);
    for (unsigned i=0; i<numIncdSnps; ++i) {
        SnpInfo *snp = incdSnpInfoVec[i];
        newWinID = orderedWinIDmap[snp->window];
//        LOGGER << i << " " << snp->window << " " << newWinID << " " << windowSnpIdxVec.size() << endl;
        windowSnpIdxVec[newWinID].push_back(snp->index);
    }
    
    windSize.resize(numWindows);
    for (unsigned i=0; i<numWindows; ++i) {
        windSize[i] = windowSnpIdxVec[i].size();
    }
    
    LOGGER << "Average window size for each SNP is " << windSize.mean() << endl;
}

void Data::filterSnpByLDrsq(const double rsqThreshold){
    LOGGER << "Filtering SNPs by rsq ... " << rsqThreshold << endl;
    unsigned numExcdSnpsBeforeFilter = 0;
    unsigned numExcdSnpsAfterFilter = 0;
    for (unsigned i=0; i<numIncdSnps; ++i) {
        SnpInfo *snpi = incdSnpInfoVec[i];
        if (!snpi->included) ++numExcdSnpsBeforeFilter;
    }
    if (sparseLDM) {
        for (unsigned i=0; i<numIncdSnps; ++i) {
            SnpInfo *snpi = incdSnpInfoVec[i];
            if (!snpi->included) continue;
            double rsq = 0.0;
            for (SparseVector <double> ::InnerIterator it(ZPZsp[i]); it; ++it) {
                SnpInfo *snpj = incdSnpInfoVec[it.index()];
                rsq = it.value()*it.value();
                if (snpi!=snpj && rsq > rsqThreshold) {
                    if (snpi->gwas_pvalue < snpj->gwas_pvalue) snpj->included = false;
                    else snpi->included = false;
                }
            }
        }
    }
    else {
        for (unsigned i=0; i<numIncdSnps; ++i) {
            SnpInfo *snpi = incdSnpInfoVec[i];
            if (!snpi->included) continue;
            VectorXd rsq = ZPZ[i].cwiseProduct(ZPZ[i]);
            for (unsigned j=0; j<numIncdSnps; ++j) {
                SnpInfo *snpj = incdSnpInfoVec[j];
                if (!snpj->included) continue;
                if (snpi!=snpj && rsq[j] > rsqThreshold) {
                    if (snpi->gwas_pvalue < snpj->gwas_pvalue) snpj->included = false;
                    else snpi->included = false;
                }
            }
        }
    }
    for (unsigned i=0; i<numIncdSnps; ++i) {
        SnpInfo *snpi = incdSnpInfoVec[i];
        if (!snpi->included) ++numExcdSnpsAfterFilter;
    }
    
    // restore the LD window information for preparing reading LD matrix again
    for (unsigned i=0; i<numSnps; ++i) {
        SnpInfo *snp = snpInfoVec[i];
        snp->windStart = snp->windStartOri;
        snp->windSize = snp->windSizeOri;
        snp->windEnd = snp->windEndOri;
    }

    LOGGER << "Removed " << numExcdSnpsAfterFilter - numExcdSnpsBeforeFilter << " SNPs with LD R-square above " << rsqThreshold << "." << endl;
}

void Data::readLDmatrixTxtFile(const string &ldmatrixFile) {
    Gadget::Timer timer;
    timer.setTime();
    
    Gadget::Tokenizer token;
    token.getTokens(ldmatrixFile, ".");
    string ldmType = token[token.size()-2];
    if (ldmType != "full") LOGGER.e(0," can only read full matrix from text file.");
    
    VectorXi windStartLDM(numSnps);
    VectorXi windSizeLDM(numSnps);
    
    windStart.resize(numIncdSnps);
    windSize.resize(numIncdSnps);
    
    SnpInfo *snpi, *snpj;
    
    for (unsigned i=0; i<numSnps; ++i) {
        SnpInfo *snpi = snpInfoVec[i];
        windStartLDM[i] = snpi->windStart;
        windSizeLDM[i]  = snpi->windSize;
    }

    ifstream in(ldmatrixFile.c_str());
    if (!in) {
        LOGGER.e(0," cannot open LD matrix file " + ldmatrixFile);
    }

    resizeWindow(incdSnpInfoVec, windStartLDM, windSizeLDM, windStart, windSize);
    
    if (numIncdSnps == 0) LOGGER.e(0,"No SNP is retained for analysis.");

    LOGGER << "Reading " + ldmType + " LD matrix from [" + ldmatrixFile + "]..." << endl;

    double rsq = 0.0;

    ZPZ.resize(numIncdSnps);
    ZPZdiag.resize(numIncdSnps);
    Gadget::Tokenizer colData;
    string inputStr;
    string sep(" \t,");
    
    for (unsigned i=0, inci=0; i<numSnps; i++) {
        snpi = snpInfoVec[i];
        getline(in,inputStr);
        
        if (!snpi->included) continue;
        
        colData.getTokens(inputStr, sep);
        
        ZPZ[inci].resize(windSize[inci]);
        snpi->ldSamplVar = 0.0;
        snpi->ldSum = 0.0;
        if (!readLDscore) snpi->ldsc = 0.0;
        snpi->numNonZeroLD = snpi->windSize;
        
        for (unsigned j=0, incj=0; j<windSizeLDM[i]; ++j) {
            unsigned idx = windStartLDM[i] + j;
            snpj = snpInfoVec[idx];
            if (snpj->included) {
                double ld = std::stod(colData[idx].c_str());
                ZPZ[inci][incj++] = ld;
                rsq = ld*ld;
                snpi->ldSamplVar += (1.0-rsq)*(1.0-rsq)/snpi->sampleSize;
                snpi->ldSum += ld;
                if (!readLDscore) snpi->ldsc += rsq;
                if (snpj == snpi)
                    ZPZdiag[inci] = ld;
            }
        }
        
        if (++inci == numIncdSnps) break;
    }
    
    in.close();
    
    timer.getTime();
    
    //    LOGGER << "Window width " << windowWidth << " Mb." << endl;
    displayAverageWindowSize(windSize);
    LOGGER << "LD matrix diagonal mean " << ZPZdiag.mean() << " sd " << sqrt(Gadget::calcVariance(ZPZdiag)) << "." << endl;
    if (ZPZdiag.mean() < 0.8 || ZPZdiag.mean() > 1.2) LOGGER.e(0," The mean of LD matrix diagonal values is expected to be close to one. Something is wrong with the LD matrix!");
    LOGGER << "Read LD matrix for " << numIncdSnps << " SNPs (time used: " << timer.format(timer.getElapse()) << ")." << endl;
}

