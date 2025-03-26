#include "Data.hpp"


void Data::includexQTLSnp(const string &includeSnpFile){
    ifstream in(includeSnpFile.c_str());
    if (!in) LOGGER.e(0,"Can not open the file [" + includeSnpFile + "] to read.");
    for (unsigned i=0; i<numeQTLs; ++i) {
        eqtlInfoVec[i]->included = false;
    }
    map<string, EqtlInfo*>::iterator it, end = eqtlInfoMap.end();
    string id;
    int line = 0;
    while (in >> id) {
        it = eqtlInfoMap.find(id);
        if (it != end) {
            it->second->included = true;
        }
        line ++;
    }
    in.close();
    LOGGER << line << " xQTL SNPs are read from [ " <<  includeSnpFile << " ]."  << endl;

}
////////// debug mode //////////////
// this function aim to extract snpID
void Data::extractParameterFromIndModel(const string bedFile){
        string outPath = bedFile;
        string outFile; std::ofstream file1, file2;
        // gwas snplist
        outFile = ".gwasSnp";
        file1.open((outPath + outFile).c_str());
        file1 << "snpID" << endl;
        for(unsigned i = 0; i < snpEffectNames.size(); i ++){
            file1 << snpEffectNames[i] << endl;
        }
        file1.close();
        // gwas individual list
        outFile = ".gwasInd";
        file1.open((outPath + outFile).c_str());
        file1 << "indID" << endl;
        for(unsigned i = 0; i < gwasIndNames.size(); i ++){
            file1 << gwasIndNames[i] << endl; // "_" << gwasIndNames[i] << endl;
        }
        file1.close();
        outFile = ".gwasGeno";
        file1.open((outPath + outFile).c_str());
        file1 << Z << endl;
        file1.close();
        // GWAS trait phenotype
        outFile = ".gwasPhe";
        file1.open((outPath + outFile).c_str());
        file1 << wbcorr << endl;
        file1.close();
        ///////////////////////////////////////////////
        ///////// pqtl parameter comparison
        if(numGenes != 0){
        // pqtlSnp Name
        outFile = ".eqtlSnp";
        file1.open((outPath + outFile).c_str());
        file1 << "snpID" << endl;
        for(unsigned i = 0; i < cisSnpIDVec.size(); i ++){
            file1 << cisSnpIDVec[i] << endl;
        }
        file1.close();
        // pqtl individual list
        outFile = ".eqtlInd";
        file1.open((outPath + outFile).c_str());
        file1 << "indID" << endl;
        IndInfo *indInfo = NULL;
        for (unsigned i = 0,indIdx =0; i < numInds; i++) {
            indInfo = indInfoVec[i];
            if (!indInfo->hasEQTL || !indInfo->hasEQTL) continue;
            file1 << indInfo->famID << endl; 
        }
        file1.close();
        // gwas ori genotype Matrix
        outFile = ".eqtlGeno";
        file1.open((outPath + outFile).c_str());
        file1 << ZGene << endl;
        file1.close();
        }
        // out gwas summary statistic 
        // gwas ori genotype Matrix
        outFile = ".sumstatis";
        file1.open((outPath + outFile).c_str());
        file1 << "SNP\tA1\tA2\tfreq\n"; 
        for(unsigned i =0; i < numIncdSnps;i++){
            SnpInfo * snp = incdSnpInfoVec[i];
            file1 << snp->rsID << "\t" << snp->a1 << "\t" << snp->a2 << "\t"
                  << snp->af << endl;
        }
        file1.close();
}

// this function aim to extract summary data
void Data::extractParameterFromSumModel(const string bedFile){
        string outPath = bedFile;
        string outFile; std::ofstream file1, file2;
        
        EqtlInfo *eqtl;
        GeneInfo *gene;
        map<string, EqtlInfo*>::iterator it;
        ////////////////////////////////////////
        // check gwas beta 
        ////////////////////////////////////////
        // gwas snplist
        outFile = ".sum.gwas.snp";
        file1.open((outPath + outFile).c_str());
        file1 << "snpID" << endl;
        for(unsigned i = 0; i < snpEffectNames.size(); i ++){
            file1 << snpEffectNames[i] << endl;
        }
        file1.close();
        for (unsigned i = 0; i < numKeptLDBlocks; i++){
            outFile = ".block." + to_string(i + 1) + ".sum.gwas.bhat";
            file1.open((outPath + outFile).c_str());
            file1 << gwasEffectInBlock[i].format(Eigen::IOFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n"));
            file1.close();
        }
        ////////////////////////////////////////
        // check pqtl alpha 
        ////////////////////////////////////////
        // gene 
        outFile = ".sum.gene";
        file1.open((outPath + outFile).c_str());
        file1 << "geneID" << endl;
        for(unsigned i = 0; i < numKeptGenes; i ++){
            gene = keptGeneInfoVec[i];
            file1 << gene->ensemblID << endl;
        }
        file1.close();
        // snplist in gene
        for(unsigned j = 0; j < numKeptGenes; j++){
            gene = keptGeneInfoVec[j];
            outFile = ".gene." + to_string(j + 1) + ".sum.eqtl.snp";
            file1.open((outPath + outFile).c_str());
            file1 << "snpID" << endl;
            for(unsigned i = 0; i < gene->cisSnpNameVec.size(); i ++){
                file1 << gene->cisSnpNameVec[i] << endl;
            }
            file1.close();
        }

        // check eqtl Ahat
        for(unsigned j = 0; j < numKeptGenes; j++){
            gene = keptGeneInfoVec[j];
            outFile = ".gene." + to_string(j + 1) + ".sum.eqtl.Ahat";
            file1.open((outPath + outFile).c_str());
            file1 << eQTLEffAcrossGenes[j].format(Eigen::IOFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n"));
            file1.close();
        }
        ////////////////////////////////////////
        // check LD 
        ////////////////////////////////////////
        // check gwas U and lambda
        for (unsigned i = 0; i < numKeptLDBlocks; i++){
            outFile = ".block." + to_string(i + 1) + ".sum.gwas.U";
            file1.open((outPath + outFile).c_str());
            file1 << eigenVecLdBlock[i].format(Eigen::IOFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n"));
            file1.close();
        }
        for (unsigned i = 0; i < numKeptLDBlocks; i++){
            outFile = ".block." + to_string(i + 1) + ".sum.gwas.lambda";
            file1.open((outPath + outFile).c_str());
            file1 << eigenValLdBlock[i].format(Eigen::IOFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n"));
            file1.close();
        }
        // check pqtl U and lambda
        for(unsigned j = 0; j < numKeptGenes; j++){
            gene = keptGeneInfoVec[j];
            outFile = ".gene." + to_string(j + 1) + ".sum.eqtl.U";
            file1.open((outPath + outFile).c_str());
            file1 << eigenVecGene[j].format(Eigen::IOFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n"));
            file1.close();
        }
        for(unsigned j = 0; j < numKeptGenes; j++){
            gene = keptGeneInfoVec[j];
            outFile = ".gene." + to_string(j + 1) + ".sum.eqtl.lambda";
            file1.open((outPath + outFile).c_str());
            file1 << eigenValGene[j].format(Eigen::IOFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n"));
            file1.close();
        }
}

////////////////////////////////////////////
void Data::readSNPRSIDMaps(const string snpRsidMapFile, const string grchType, const string title){
    std::ifstream file(snpRsidMapFile, std::ios_base::in | std::ios_base::binary);
    if (!file) LOGGER.e(0, " can not open the file [" + snpRsidMapFile + "] to read.");
    LOGGER << "Reading snp info file from [" + snpRsidMapFile + "]." << endl;
    boost::iostreams::filtering_istream in;
    in.push(boost::iostreams::gzip_decompressor());
    in.push(file);
    map<string, EqtlInfo*>::iterator iterEqtl;
    map<string, EqtlInfo*>::iterator iterEqtlrsID;
    string header, id,effAllele, nonEffAllele,rsid;
    double pos19, pos38, phyPos; 
    int chr;
    getline(in, header);
    // LOGGER << "header: " << header << endl;
    string sep(":");
    Gadget::Tokenizer colData;
    string inputStr;
    unsigned idx = 0;
    unsigned dupIdx = 0;
    eqtlInfoMap.clear();
    eqtlInfoVec.clear();
    eqtlInfoAlterMap.clear();
    chromosomes.clear();
    set<string> dupSnpLists;
    // id ALLELE0 ALLELE1
    while (in >> id >> nonEffAllele >> effAllele >> rsid >> pos19 >> pos38) {
        colData.getTokens(id, sep);
        chr = std::stoi(colData[0]);
        // GrCh38 (hg38) or GrCh37 (hg19)
        if(grchType == "hg38") {
            phyPos = pos38;
        } else {
            phyPos = pos19;
        }

        EqtlInfo *eqtl = new EqtlInfo(idx++, rsid,effAllele,nonEffAllele, chr, 0, phyPos,0);
        eqtl->alterId = id;
        eqtl->included = false; 
        
        if (eqtlInfoAlterMap.insert(pair<string, EqtlInfo*>(id, eqtl)).second == false) {
            LOGGER.e(0, " Duplicate eQTL SNP ID found: : [" + id + "].");
            // eqtl->isInGene = false;
            // eqtl->included = false; // since it's only a reference
        }
        if (eqtlInfoMap.insert(pair<string, EqtlInfo*>(rsid, eqtl)).second == false) {
            // LOGGER.w(0, " Duplicate eQTL SNP ID found: rsid: [" + rsid + "].");
            dupSnpLists.insert(eqtlInfoMap[rsid]->alterId);
            dupSnpLists.insert(id);
            // eqtl->included = false;
        }
        eqtlInfoVec.push_back(eqtl);
        chromosomes.insert(eqtl->chrom);
    }

    file.close();
    dupIdx = dupSnpLists.size();
    numeQTLs = (unsigned) eqtlInfoVec.size();
    LOGGER << numeQTLs - dupIdx << " eQTLs to be included from [" + snpRsidMapFile + "]." << endl;

    if(dupSnpLists.size() != 0){
        string dupFile = title + ".duplicateID";
        std::ofstream file1(dupFile.c_str());
        set<string>::iterator iterSet;
        LOGGER.w(0, std::to_string(dupIdx) + " SNP with duplicate rsid are found and stored in [ " + dupFile + " ]");
        file1 << "ID" << endl;
        
        for(iterSet = dupSnpLists.begin(); iterSet != dupSnpLists.end();iterSet ++){
            file1 << *iterSet << endl;
            iterEqtl = eqtlInfoAlterMap.find(*iterSet);
            iterEqtlrsID = eqtlInfoMap.find(iterEqtl->second->rsID);
            if(iterEqtl != eqtlInfoAlterMap.end())eqtlInfoAlterMap.erase(iterEqtl);
            if(iterEqtlrsID != eqtlInfoMap.end()) eqtlInfoMap.erase(iterEqtlrsID);
        }
        file1.close();
    }

}
void Data::readProMaps(const string proteinMapFile, const double pValueThreshold, const string title){
    ifstream in(proteinMapFile.c_str());
    if (!in) LOGGER.e(0, " can not open the file [" + proteinMapFile + "] to read.");
    LOGGER << "Reading protein map  file from [" + proteinMapFile + "]." << endl;
    geneInfoMap.clear();
    geneInfoMap.clear();

    string header;
    getline(in, header);

    Gadget::Tokenizer colData;
    string inputStr;
    string sep("\t");

    unsigned chr, physPos, sampleSize;
    double geneticDis;
    string probe, geneID, geneOri,proPath;
    unsigned idx=0;
    unsigned nonAutoIdx = 0;
    vector<string> NonChrLists;
    vector<string> cisSnpNameVec;
    vector<double>  eQTLMarginEffect;   // eQTL marginal effect
    vector<double> eQTLMarginEffectSE; // eQTL standard error of marginal effect
    map<string,int> cisSnpID2IdxMapInGene;
    set<string> cisSnpIDSetInGene;

    while (getline(in,inputStr)) {
        if (inputStr.size() == 0) continue;            
        cisSnpNameVec.clear();
        eQTLMarginEffect.clear();
        eQTLMarginEffectSE.clear();
        cisSnpID2IdxMapInGene.clear();
        cisSnpIDSetInGene.clear();
        bool hasEqtl = false;
        //////// start
        colData.getTokens(inputStr, sep);
        // remove non-autosomal chromosome
        // LOGGER << "ID: " << colData[0] << endl;
        geneID = colData[0];
        // LOGGER << "chr: " << colData[1] << endl;
        if(!Gadget::isNumber(colData[1].c_str())){
            NonChrLists.push_back(colData[1]);
            nonAutoIdx++;
            continue;
        }
        chr = std::stoi(colData[1]);
        GeneInfo *gene = new GeneInfo(idx++, geneID, chr);
        gene->start = std::stoi(colData[2]);
        gene->end = std::stoi(colData[3]);
        gene->midPhyPos =  (gene->start + gene->end)/ 2.0;  // here we ued middle position as cis-region
        // cout << "gene: " << gene->ensemblID << " start: " << gene->start << " end: "  << gene->end << " pos: " << gene->midPhyPos << endl;
        proPath = colData[4];
        if (geneInfoMap.insert(pair<string, GeneInfo*>(geneID, gene)).second == false) {
            LOGGER.e(0, " Duplicate gene ID found: \"" + geneID + "\".");
        }
        gene->hasEqtl = false;
        readGWASFromProtein(gene->ensemblID,cisSnpNameVec,eQTLMarginEffect,eQTLMarginEffectSE, cisSnpID2IdxMapInGene,cisSnpIDSetInGene, proPath, pValueThreshold, hasEqtl);
        gene->hasEqtl = hasEqtl;
        gene->cisSnpNameVec = cisSnpNameVec;
        gene->cisSnpID2IdxMapInGene = cisSnpID2IdxMapInGene;
        gene->eQTLMarginEffect = Eigen::Map<VectorXd>(&eQTLMarginEffect[0], eQTLMarginEffect.size());
        gene->eQTLMarginEffectSE = Eigen::Map<VectorXd>(&eQTLMarginEffectSE[0], eQTLMarginEffectSE.size());
        gene->cisSnpIDSetInGene = cisSnpIDSetInGene;
        gene->numSnpInGene = cisSnpNameVec.size();
        if(gene->numSnpInGene == 0) {
            gene->kept = false;
            LOGGER.w(0,"Could not find gene [ " +gene->ensemblID + " ].");
        }
        if(gene->cisSnpNameVec.size()!= gene->cisSnpIDSetInGene.size()){
            LOGGER.e(0,"Error, the size of cisSnpNameVec should be equal to size of cisSnpIDSetInGene for gene: [ " +gene->ensemblID + " ].");
        }
        geneInfoVec.push_back(gene);
    } // end of while
    in.close();
    numGenes = (unsigned) geneInfoVec.size();
    LOGGER << numGenes << " genes to be included from [" + proteinMapFile + "]." << endl;

    if(NonChrLists.size() != 0){
        string nonAutoFile = title + ".nonAutosome";
        std::ofstream file1(nonAutoFile.c_str());
        LOGGER.w(0, std::to_string(nonAutoIdx) + " non-autosome are found and stored in [ " + nonAutoFile + " ]");
        file1 << "Chromosome" << endl;
        for(unsigned i=0; i < NonChrLists.size();i++){
            file1 << NonChrLists[i] << endl;
        }
        file1.close();
    }
}
void Data::readGWASFromProtein(const string geneID, vector<string> & cisSnpNameVec, vector<double>  &eQTLMarginEffect, vector<double> &eQTLMarginEffectSE, map<string,int> &cisSnpID2IdxMapInGene, set<string> & cisSnpIDSetInGene, const string proPath,const double pValueThreshold, bool &hasEqtl){
    /// Collect the 
    /////
    namespace fs = boost::filesystem; 
    fs::path targetDir(proPath); 
    map<string, EqtlInfo*>::iterator iterEqtl;
    map<string, EqtlInfo*>::iterator iterEqtlrsID;
    set<int>::iterator iterSet;

    string header, id, allele0, allele1, extra,test;
    double  phyPos,a1freq, info, chisq, log10p; 
    double beta,se;
    int chr,localChr, sampleSize;
    vector<double> sampleSizeVec;
    string sep(":");
    Gadget::Tokenizer colData;
    string inputStr;
    unsigned idx = 0;
    unsigned numFiles,numTotalFiles= 0, numFolders = 0;
    hasEqtl = false;
    bool skipOtherFiles = false;
    bool hasFileInFloder= false;
    string regexStr;
    GeneInfo *gene;
    gene = geneInfoMap.find(geneID)->second;
    if(fs::is_directory(targetDir)) {
        LOGGER << targetDir << " is a directory containing:\n";
        for(fs::directory_entry & singleGZFile : boost::make_iterator_range(fs::directory_iterator(targetDir), {})){
            numFiles = 0;
            fs::path filePath = singleGZFile.path().string();
            // cout << " file: " << singleGZFile.path().string() << endl;
            // regex chromosome
            for(iterSet = chromosomes.begin(); iterSet != chromosomes.end();iterSet ++){
                regexStr = ".*_chr" + to_string(*iterSet) + "_.*";
                boost::regex pattern(regexStr);
                if(boost::regex_match(filePath.string(), pattern)){
                    // if(*iterSet == gene->chrom){
                        hasFileInFloder = true;  // here we only select snps on same genes;
                    // }
                }
            }

            if(!hasFileInFloder) continue;
            
            if (filePath.extension() == ".gz") // Heed the dot.
            {
                // if(!hasFileInFloder) hasFileInFloder = true;
                std::ifstream file(singleGZFile.path().string(), std::ios_base::in | std::ios_base::binary);
                if (!file) LOGGER.e(0, " can not open the file [" + singleGZFile.path().string() + "] to read.");
                boost::iostreams::filtering_istream in;
                in.push(boost::iostreams::gzip_decompressor());
                in.push(file);
                getline(in, header);
                // header: CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTRA
                while (in >> chr >> phyPos >> id >> allele0 >> allele1 >> a1freq >> info >> sampleSize >> test >> beta >> se >> chisq >> log10p >> extra) {
                    // LOGGER << " chr: "  << chr << endl;
                    if(chr == 23 ){continue;} // skip chr for accelerating code
                    if(chromosomes.find(chr) == chromosomes.end()){continue;;} 
                    iterEqtl = eqtlInfoAlterMap.find(id);
                    if(iterEqtl != eqtlInfoAlterMap.end() ){
                        if(!hasEqtl) {hasEqtl = true;}
                        if(!skipOtherFiles){skipOtherFiles = true;}
                        iterEqtl->second->isInGene = true;
                        iterEqtl->second->included = true;
                        iterEqtlrsID = eqtlInfoMap.find(iterEqtl->second->rsID);
                        iterEqtl->second->eqtl_b = beta;
                        iterEqtl->second->eqtl_se = se;
                        iterEqtl->second->eqtl_n = sampleSize;
                        sampleSizeVec.push_back(sampleSize);
                        iterEqtl->second->samVec.push_back(sampleSize);
                        iterEqtl->second->af = a1freq;
                        iterEqtl->second->afVec.push_back(a1freq);
                        iterEqtl->second->eqtl_pvalue = pow(10,log10p);
                        // add to protein
                        cisSnpNameVec.push_back(iterEqtl->second->rsID);
                        cisSnpID2IdxMapInGene.insert(pair<string, int>(iterEqtl->second->rsID,idx));
                        cisSnpIDSetInGene.insert(iterEqtl->second->rsID);
                        eQTLMarginEffect.push_back(beta);
                        eQTLMarginEffectSE.push_back(se);
                        idx ++;
                        // if(iterEqtl->second->rsID == "rs201540471") {
                        //     cout << "af: " << a1freq << " sam: " << sampleSize << endl;
                        // }
                    }
                } // end of while loop
                file.close();
                numFiles ++;
                numTotalFiles ++;
                LOGGER << idx << " cis-SNPs " << numFiles << " files have been found in folder [" <<  singleGZFile.path().string() << "]" << endl;
            } // end of file check of one folder
            if(sampleSizeVec.size() != 0) {
                gene->sampleSize = Gadget::calcMean(sampleSizeVec);
            }
            if(hasFileInFloder) numFolders ++;
            if (skipOtherFiles) {
                break; // exit loop
            }
        } // end of loop in the folder
    } // end of current dir if(fs::is_directory(targetDir))
    // incdEqtlInfoVec = makeIncdEqtlInfoVec(eqtlInfoVec,false);
    // numIncdEqtls = (unsigned) incdEqtlInfoVec.size();
}

void Data::readFlistFile(const string geneInfoFile,const string subGenePath){
    if (geneInfoFile.empty()) return;
    ifstream in(geneInfoFile.c_str());
    if (!in) LOGGER.e(0, " can not open the gene info file [" + geneInfoFile + "] to read.");
    LOGGER << "Reading gene info file from [" + geneInfoFile + "]." << endl;

    geneInfoVec.clear();
    geneInfoMap.clear();

    vector<string> NonChrLists;
    vector<string> cisSnpNameVec;
    vector<double>  eQTLMarginEffect;   // eQTL marginal effect
    vector<double> eQTLMarginEffectSE; // eQTL standard error of marginal effect
    map<string,int> cisSnpID2IdxMapInGene;
    set<string> cisSnpIDSetInGene;

    unsigned chr;
    int start, end;
    string probe, geneID, Orientation,singleGenePheFile;
    unsigned idx = 0;
    double GeneticDistance,ProbeBp,sampleSize;

    Gadget::Tokenizer colData;
    string header;
    string sep(" \t\r");
    getline(in, header);
    int eqtlidx = 0;
    // Step 1. read plist file, where information per gene was stored ( chr, gene start, gene end, ensgid ID, genePath)
    // genePath is absolute path;
    // So when genotypes from gene cis-regions are needed, we should map genotypes to gene based on gene cis-windows. 
    // Chr	ProbeID	GeneticDistance	ProbeBp	Gene	Orientation	PathOfEsd
    while (in >> chr >> probe >> GeneticDistance >> ProbeBp >> geneID >> sampleSize >> Orientation >> singleGenePheFile) {
        bool hasEqtl = false;
        GeneInfo *gene = new GeneInfo(idx++, probe, chr);
        gene->midPhyPos = ProbeBp;
        gene->sampleSize = sampleSize;
        geneInfoVec.push_back(gene);
        if (geneInfoMap.insert(pair<string, GeneInfo*>(geneID, gene)).second == false) {
            LOGGER.e(0, " Duplicate gene ID found: \"" + geneID + "\".");
        }
        // Step 2. read phenotype information from each gene. Here based on genePath from plist file, we located gene file first,
        // and then extract relevant genes. Please note, individuals from each gene should be a subset of individuals from plink fam file,
        // and don't need to be consistent with individuals from complex trait, considering the possibility that gene sample size may be 
        // not equal to sample size of complex trait.
        readSingleEsdFile(gene->ensemblID, singleGenePheFile,eqtlidx,subGenePath);
    }
    // gene info summary
    numGenes = (unsigned) geneInfoVec.size();
    numeQTLs = (unsigned) eqtlInfoVec.size();
    LOGGER << numGenes << " genes with " << numeQTLs << " eQTLs to be included from [" + geneInfoFile + "]." << endl;
}

void Data::readSingleEsdFile(const string geneID,const string singleGenePheFile,int eqtlIdx,const string subGenePath){
    string filename;
    if (singleGenePheFile.empty()) return;
    if(!subGenePath.empty()){
        filename = subGenePath + "/" + singleGenePheFile;
    } else {
        filename = singleGenePheFile;
    }
    ifstream in(filename.c_str());
    if (!in) LOGGER.e(0, " can not open the file [" + filename + "] to read.");
    LOGGER << "Reading eQTLs within gene " +  geneID  + " from [" + filename + "]." << endl;

    EqtlInfo * eqtl;
    map<string, EqtlInfo*>::iterator iterEqtl;
    GeneInfo *gene;
    map<string, GeneInfo*>::iterator iterGene;
    vector<double> eQTLMarginEffect;
    vector<double> eQTLMarginEffectSE;
    string header;
    string id, allele1, allele2;
    int chr, physPos;
    double genPos = 0;
    double allele1Freq,beta,se, eqtl_n;
    int idx = 0;
    int line = 0;
    double p;
    getline(in, header);
    iterGene = geneInfoMap.find(geneID);
    gene = iterGene->second;
    string sep("\t ");
    Gadget::Tokenizer colData;
    string inputStr;
    while (getline(in,inputStr)) {
        colData.getTokens(inputStr, sep);
        if(colData.size() != 9 ){
            LOGGER.e(0,"the number of line "+ to_string(line)+ " is less than 9 in file [" + filename + "].");
        }
        chr = stoi(colData[0]);
        id = colData[1];
        if(id == "rs4920459"){
            cout << "here" << endl;
        }
        physPos = stoi(colData[2]);
        allele1 = colData[3];
        allele2 = colData[4];
        allele1Freq = stod(colData[5]);
        beta = stod(colData[6]);
        se = stod(colData[7]);
        cpp_dec_float_100 highPrecisionValue = cpp_dec_float_100(colData[8]);
        p = static_cast<double>(highPrecisionValue);
  
        genPos = 0;
        /////////////////////
        iterEqtl = eqtlInfoMap.find(id);
        if(iterEqtl != eqtlInfoMap.end()){ //found
            eqtl = iterEqtl->second;
            eqtl->isInGene = true;
            // eqtl->eqtl_n = eqtl_n;
        } else{
            eqtl = new EqtlInfo(idx++, id, allele1, allele2, chr, genPos, physPos,allele1Freq);
            eqtl->eqtl_n = eqtl_n;
            eqtl->isInGene = true;
            eqtlInfoVec.push_back(eqtl);
            if (eqtlInfoMap.insert(pair<string, EqtlInfo*>(id, eqtl)).second == false) {
                LOGGER.w(0, " Duplicate SNP ID found: \"" + id + "\".");
            }
        } 
        if(gene->cisSnpIDSetInGene.find(eqtl->rsID) != gene->cisSnpIDSetInGene.end()){
            LOGGER.w(0, " Duplicate SNP ID found: \"" + id + "\" in the file [ " + filename + " ].");
        }
        eQTLMarginEffect.push_back(beta);
        eQTLMarginEffectSE.push_back(se);
        gene->hasEqtl = true;
        gene->cisSnpNameVec.push_back(eqtl->rsID); 
        gene->cisSnpID2IdxMapInGene.insert(pair<string, int>(eqtl->rsID,line));
        gene->cisSnpIDSetInGene.insert(eqtl->rsID);
        gene->numSnpInGene = gene->cisSnpIDSetInGene.size();
        line ++;
    }
    in.close();

    gene->eQTLMarginEffect = Eigen::Map<VectorXd>(&eQTLMarginEffect[0], eQTLMarginEffect.size());
    gene->eQTLMarginEffectSE = Eigen::Map<VectorXd>(&eQTLMarginEffectSE[0], eQTLMarginEffectSE.size());
    int numeQTLPerGene = eQTLMarginEffect.size();
    LOGGER << numeQTLPerGene << " eQTLs to be included from [" + filename + "]." << endl;
}
//
void Data::readEsiFile(const  string &esiFile, bool &smrBesdBool,const bool haveSnpInfo) {
    // Read bim file: recombination rate is defined between SNP i and SNP i-1
    ifstream in(esiFile.c_str());
    if (!in) LOGGER.e(0, " can not open the file [" + esiFile + "] to read.");
    LOGGER << "Reading BESD format esi file from [" + esiFile + "]." << endl;
    eqtlInfoVec.clear();
    eqtlInfoMap.clear();
    string id, allele1, allele2;
    unsigned chr, physPos;
    double genPos;
    double freq;
    double af;
    int sampleSize;
    unsigned idx = 0;
    unsigned line = 0;
    int testIdx = 0;
    // bool smrBesdBool = false;
    SnpInfo * snp;
    bool warnNullFlag1 =false, warnNullFlag2 =false, warnNullFlag3 =false,warnNullFlag4 =false;
    bool warnNullFlag5 =false, warnNullFlag6 =false, warnNullFlag7 =false,warnNullFlag8 =false;
    Gadget::Tokenizer colData;string inputStr;string sep(" \t\n");
    while (getline(in,inputStr)) {
        colData.getTokens(inputStr, sep);
        if(idx == 0){
            if(colData.size() == 7 ){
                smrBesdBool = true;
                LOGGER << "SMR BESD format is used and per snp sample size is missing." << endl;
            } else if (colData.size() == 8) {
                smrBesdBool = false;
                LOGGER << "A modified SMR BESD format is used and per snp sample size can be read in esi file." << endl;
            } else {
                LOGGER.e(0,"Wrong esi format.");
            }
        }
        line ++;
        id = colData[1];
        // chr
        if(boost::algorithm::to_upper_copy(colData[0])=="NA" || colData[0]=="0"){
            if(!warnNullFlag1)LOGGER.w(0,"One or more chr name is \'NA\'.");
            warnNullFlag1 = true;
            chr = -9;
        } else {
            if(boost::algorithm::to_upper_copy(colData[0]) == "X") chr = 23;
            else if(boost::algorithm::to_upper_copy(colData[0]) == "Y") chr = 24;
            else chr =std::stoi(colData[0]);
        }
        // snp
        if(boost::algorithm::to_upper_copy(colData[1]) == "NA"){
            if(!warnNullFlag2) LOGGER.e(0,"the SNP name is \'NA\' in row " + to_string(line) + ".");
            // warnNullFlag2 = true;
        }
        // cm
        if(boost::algorithm::to_upper_copy(colData[2])=="NA"){
            if(!warnNullFlag3)LOGGER.w(0," One or more genetic distance is \'NA\'.");
            warnNullFlag3 = true;
            genPos = -9;
        } else {
            genPos = std::stod(colData[2]);
        }
        // pos
        if(boost::algorithm::to_upper_copy(colData[3])=="NA"){
            if(!warnNullFlag4)LOGGER.e(0,"Physical position is \'NA\' in row" + to_string(line) +".");
            // warnNullFlag4 = true;
            physPos = -9;
        } else {
            physPos = std::stoi(colData[3]);
        }
        // a1
        if(boost::algorithm::to_upper_copy(colData[4])=="NA"){
            if(!warnNullFlag5)LOGGER.e(0,"Allele 1 is \'NA\' in row" + to_string(line) +".");
            // warnNullFlag5 = true;
            allele1 = "NA";
        } else {
            allele1 = colData[4];
        }
        // a2
        if(boost::algorithm::to_upper_copy(colData[5])=="NA"){
            if(!warnNullFlag6) LOGGER.e(0,"Allele 2 is \'NA\' in row" + to_string(line) +".");
            // warnNullFlag6 = true;
            allele2 = "NA";
        } else {
            allele2 = colData[5];
        }
        // freq
        if(boost::algorithm::to_upper_copy(colData[6])=="NA"){
            if(!warnNullFlag7) LOGGER.w(0,"One or more allele frequency of allele 1 is \'NA\'.");
            warnNullFlag7 = true;
            freq = -999;
        } else {
            freq = std::stod(colData[6]);
        }
        // // N
        if(!smrBesdBool){
            if(boost::algorithm::to_upper_copy(colData[7])=="NA"){
                if(!warnNullFlag8)LOGGER.w(0,"One or more per snp sample size is \'NA\'.");
                warnNullFlag8 = true;
            } else {
            }
        }

        EqtlInfo *eqtl = new EqtlInfo(idx++, 
        id, // id, 
        allele1, 
        allele2, 
        chr, 
        genPos, 
        physPos,
        freq
        );
        if(!smrBesdBool)eqtl->eqtl_n = std::stoi(colData[7]);
        eqtlInfoVec.push_back(eqtl);
        if (eqtlInfoMap.insert(pair<string, EqtlInfo*>(id, eqtl)).second == false) {
            LOGGER.e(0, " Duplicate eQTL SNP ID found: \"" + id + "\".");
        }
        // checkt whether eqtl exists in GWAS snplist
        if(haveSnpInfo){
            eqtl->included = false;
            if(snpInfoMap.find(eqtl->rsID) != snpInfoMap.end()){
                snp = snpInfoMap.find(eqtl->rsID)->second;
                if(snp->included) eqtl->included = true;
            }
        }

    }
    in.close();
    numeQTLs = (unsigned) eqtlInfoVec.size();
    LOGGER << numeQTLs << " eQTLs to be included from [" + esiFile + "]." << endl;
}
void Data::readEpiFile(const  string &epiFile){
    ifstream in(epiFile.c_str());
    if (!in) LOGGER.e(0, " can not open the file [" + epiFile + "] to read.");
    LOGGER << "Reading BESD format epi file from [" + epiFile + "]." << endl;
    geneInfoVec.clear();
    geneInfoMap.clear();
    unsigned chr, sampleSize;
    double geneticDis; //,physPos;
    unsigned physPos;
    string probe, geneID, geneOri;
    unsigned idx = 0;
    // while (in >> chr >> geneID >> geneticDis >> physPos >> probe >> geneOri) {
        // GeneInfo *gene = new GeneInfo(idx++, geneID, chr);
        // gene->probeID = probe;
        // gene->geneticDis = geneticDis;
        // gene->start  = geneticDis;
        // gene->midPhyPos = physPos;
        // gene->geneOri = geneOri;
        // gene->sampleSize = sampleSize;
    Gadget::Tokenizer colData;string inputStr;string sep(" \t\n");
    while (getline(in,inputStr)) {
        colData.getTokens(inputStr, sep);
        if(idx == 0){
            if(colData.size() == 7){
                LOGGER << "sample size are stored in epi file." << endl;
            } else if(colData.size() == 6){
                LOGGER << "sample size are NOT stored in epi file." << endl;
            } else {
                LOGGER.e(0,"Wrong epi format.");
            }
        }
        chr = std::stoi(colData[0]);
        geneID = colData[1];
        GeneInfo *gene = new GeneInfo(idx++, geneID, chr);
        gene->geneticDis = std::stod(colData[2]);
        gene->start  = std::stod(colData[2]);
        gene->midPhyPos = std::stod(colData[3]);
        gene->end = 2 * gene->midPhyPos - gene->start;
        gene->probeID = colData[4];
        if(colData.size() == 6){
            gene->geneOri = colData[5];
        } else if(colData.size() == 7){
            gene->sampleSize = std::stoi(colData[5]);
            gene->geneOri = colData[6];
        }
        geneInfoVec.push_back(gene);
        if (geneInfoMap.insert(pair<string, GeneInfo*>(geneID, gene)).second == false) {
            LOGGER.e(0, " Duplicate gene ID found: \"" + geneID + "\".");
        }
    }
    in.close();
    numGenes = (unsigned) geneInfoVec.size();
    LOGGER << numGenes << " genes to be included from [" + epiFile + "]." << endl;
}
void Data::readBesdFile(const string &besdFile,const bool &hasSnpInfo,const bool &makeGeneLDEigenBool, bool &smrBesdBool){   
    keptGeneInfoVec  = makeKeptGeneInfoVec(geneInfoVec);
    numKeptGenes = (unsigned) keptGeneInfoVec.size();
    incdEqtlInfoVec = makeIncdEqtlInfoVec(eqtlInfoVec);
    numIncdEqtls = (unsigned) incdEqtlInfoVec.size();
    
    if(numIncdEqtls == 0) LOGGER.e(0," No SNP is retained for eQTL analysis.");
    if(numKeptGenes == 0) LOGGER.e(0," No gene is retained for eQTL analysis. Please check gene difference between BESD .epi file and gene LD matrix .ldm.gene.info file.");
    ifstream besd(besdFile.c_str(), ios::in|ios::binary);
    if (!besd) LOGGER.e(0, " can not open the file [" + besdFile + "] to read.");
    LOGGER << "Reading eQTL summary data from BESD format BED file [" + besdFile + "] ..." << endl;
    
    char SIGN[sizeof(uint64_t)+8];
    besd.read(SIGN,4);
    uint32_t header = *(uint32_t *)SIGN;
    if(header == 0x40000000 & header == 0x3f800000){
         LOGGER.e(0, "This is an old BESD format. Please use smr --make-besd command to update the file format." );
         exit(1);
    }
    // Read eQTL summary 
	map<string, SnpInfo *>::iterator iterSnp;
    EqtlInfo *eqtl = NULL;
    GeneInfo *gene = NULL;
    unsigned geneIdx,eqtlIdx;
    float eqtlBeta;
    float eqtlSE;
    int sampleSize;
    unsigned long long skip = 0;

    /////////////////////////////////////////////////////////////
    float totalSteps = 0; // Total number of steps in your process
    /////////////////////////////////////////////////////////////
    // read dense format
    if(header == DENSE_FILE_TYPE_1 || header == DENSE_FILE_TYPE_3 ){
        LOGGER << "Read eqtl info from smr dense besd format." << endl;
        unsigned size = (numeQTLs)<<3;
        // check the completeness of besd
        int length=(RESERVEDUNITS - 1)*sizeof(int);
        char* indicators=new char[length];
        besd.read(indicators,length);
        int* tmp=(int *)indicators;
        sampleSize=*tmp++;
        if(sampleSize!=-9){
            LOGGER << "The eQTL sample size is " << sampleSize << endl;
            numKeptIndseQTL = sampleSize;
        } else {
            // LOGGER.e(0,"ERROR:The sample size is missing. You may use (smr --beqtl-summary myeqtl --add-n 1000 --make-besd --out mybesd) to add it to the BESD file.");
            LOGGER << "Warning:The sample size is missing. You may use (smr --beqtl-summary myeqtl --add-n 1000 --make-besd --out mybesd) to add it to the BESD file." << endl;
        }
        if(*tmp++!=numeQTLs){
            LOGGER.e(0,"ERROR: The SNPs in your .esi file are not in consistency with the one in .besd file.");
        }
        if(*tmp++!=numGenes)
        {
            LOGGER.e(0,"ERROR: The probes in your .epi file are not in consistency with the one in .besd file");
        }
        delete[] indicators;
        // begin to read eqtl beta and se
        for (unsigned j = 0, geneIdx = 0; j < numGenes; j++){
            gene = geneInfoVec[j];
            gene->gene2CisSnpVec.clear();
            if(!gene->kept){
                skip +=size;
                continue;
            }
            if(skip) besd.seekg(skip);
            skip = 0;
            char *besdLineIn = new char[size];
            float * ft;
            float * se_ptr;
            besd.read(besdLineIn,size);
            ft=(float *)besdLineIn;
            VectorXd eqtlEffect(numIncdEqtls), eqtlEffectSE(numIncdEqtls);
            eqtlEffect.setZero();
            eqtlEffectSE.setZero();
            for (unsigned i =0; i < numeQTLs; i++){
                eqtl = eqtlInfoVec[i];
				// if this eqtl can not be found in ld snp reference, this will be removed.
				if(hasSnpInfo){
                    if(snpInfoMap.find(eqtl->rsID) == snpInfoMap.end()) eqtl->included = false;
                } 
                if(!eqtl->included) continue;
                eqtlEffect(i) =  *(ft + eqtl->index);
                se_ptr = ft + numeQTLs;
                eqtlEffectSE(i) = *(se_ptr + eqtl->index);
                if(eqtlEffect(i) == -9 || abs(eqtlEffect(i)) < 5e-6){
                    cout << "eqtl: " << eqtl->rsID << endl;
                }
                if(eqtlEffect(i) != -9 && abs(eqtlEffect(i)) != 0){
                    gene->cisSnpNameVec.push_back(eqtl->rsID);
                    gene->gene2CisSnpVec.push_back(i);
                    gene->hasEqtl = true;
                    eqtl->isInGene = true;
                }
                // For this part, we cannnot use if else statement to remove eqtls since overlapping eQTLs exist.
                // } else {
                //     eqtl->included = false;
                // }
            }
            gene->numSnpInGene = gene->gene2CisSnpVec.size();
            if(!makeGeneLDEigenBool){
                gene->eQTLMarginEffect = eqtlEffect(gene->gene2CisSnpVec);
                gene->eQTLMarginEffectSE = eqtlEffectSE(gene->gene2CisSnpVec);
            }
            delete [] besdLineIn;

            // remove genes if there is no eqtl in gene 
            if (gene->gene2CisSnpVec.size() == 0){ gene->kept = false;}
            if (++geneIdx == numGenes) break;
        }

    }else if (header  == SPARSE_FILE_TYPE_3F || header == SPARSE_FILE_TYPE_3){
        // LOGGER << "Read eqtl info from smr sparse besd format." << endl;
        // clear datastruct for dense befor read sparse
        char* buffer;
        uint64_t colNum=(numGenes<<1)+1;
        uint64_t valNum;
        uint64_t lSize;
        besd.seekg(0,besd.end);
        lSize = besd.tellg();
        besd.seekg(4); //  same as besd.seekg(4, besd.beg);
        
        if(header==SPARSE_FILE_TYPE_3){
            int length=(RESERVEDUNITS-1)*sizeof(int);
            char* indicators=new char[length];
            besd.read(indicators,length);
            int* tmp=(int *)indicators;
            sampleSize=*tmp++;
            if(sampleSize!=-9){
                LOGGER << "The eQTL sample size is " << sampleSize << endl;
                numKeptIndseQTL = sampleSize;
            }else {
                // LOGGER.e(0,"ERROR:The sample size is missing. You may use (smr --beqtl-summary myeqtl --add-n 1000 --make-besd --out mybesd) to add it to the BESD file.");
                LOGGER << "Warning:The sample size is missing. You may use (smr --beqtl-summary myeqtl --add-n 1000 --make-besd --out mybesd) to add it to the BESD file." << endl;
            }
            if(*tmp++!=numeQTLs){
                LOGGER.e(0,"ERROR: The SNPs in your .esi file are not in consistency with the one in .besd file.");
            }
            if(*tmp++!=numGenes)
            {
                LOGGER.e(0,"ERROR: The probes in your .epi file are not in consistency with the one in .besd file");
            }
            delete[] indicators;
        }
        besd.read(SIGN, sizeof(uint64_t));
        valNum=*(uint64_t *)SIGN;

        if(header==SPARSE_FILE_TYPE_3F) {
            if( lSize - (sizeof(uint32_t) + sizeof(uint64_t) + colNum*sizeof(uint64_t) + valNum*sizeof(uint32_t) + valNum*sizeof(float)) != 0){
                // cout << " total size: " << sizeof(uint32_t) + sizeof(uint64_t) + colNum*sizeof(uint64_t) + valNum*sizeof(uint32_t) + valNum*sizeof(float) << endl;
                LOGGER << "The file size is " <<  lSize;
                LOGGER << " " << to_string(sizeof(uint32_t)) << " " <<  to_string(sizeof(uint64_t)) << " " << to_string(colNum*sizeof(uint64_t)) << " "
                     << to_string(valNum*sizeof(uint32_t)) << " " << to_string(valNum*sizeof(float)) << endl;
                LOGGER.e(0,"ERROR: wrong value number. File" + besdFile + " is ruined.\n");
            }
        }else{
            if( lSize - (RESERVEDUNITS*sizeof(int) + sizeof(uint64_t) + colNum*sizeof(uint64_t) + valNum*sizeof(uint32_t) + valNum*sizeof(float)) != 0){
                LOGGER << "The file size is " <<  lSize;
                LOGGER << " " << to_string(RESERVEDUNITS*sizeof(int)) << " " << to_string(sizeof(uint64_t)) << " " << to_string(colNum*sizeof(uint64_t)) << " "
                     << to_string(valNum*sizeof(uint32_t)) << " " << to_string(valNum*sizeof(float)) << endl;
                LOGGER.e(0,"ERROR: wrong value number. File" + besdFile + "is ruined.\n");
            }    
        }
        //for sparse
        vector<uint64_t> _cols;
        vector<uint32_t> _rowid;
        //vector <double>  _val;
        // VectorXd _val;
        vector<float> _val;
        //////////////////////
        // Step 1, read all eqtl data

        buffer = (char*) malloc  (sizeof(char)*(lSize));
        if (buffer == NULL) {fputs ("Memory error",stderr); exit (1);}
        besd.read(buffer,lSize);
        if(header ==SPARSE_FILE_TYPE_3F)
        {
        if (besd.gcount()+sizeof(uint32_t) + sizeof(uint64_t) != lSize){LOGGER.e(0," cannot read besd file: " + besdFile );}
        }else {
            if (besd.gcount()+RESERVEDUNITS*sizeof(int) + sizeof(uint64_t) != lSize){ LOGGER.e(0," cannot read besd file: " + besdFile );}
        }
        
        _cols.resize(colNum);
        _rowid.resize(valNum);
        _val.resize(valNum); 
        uint64_t* ptr;
        ptr=(uint64_t *)buffer;
        for(int i=0;i<colNum;i++) {_cols[i]=*ptr++;} 
        uint32_t* ptr4B=(uint32_t *)ptr;
        for(int i=0;i<valNum;i++){_rowid[i]=*ptr4B++;}
        float* val_ptr=(float*)ptr4B;
        for(int i=0;i<valNum;i++){_val[i] =*val_ptr++;}

        // Step 2. build a map to save eqtls passed QC

        map<int,string> eqtlIndxMap;

        map<int, string>::iterator Intiter;
        map<string, EqtlInfo* >::iterator iterEqtl;
		map<string, SnpInfo *>::iterator iterSnp;

        for(int i = 0; i < numeQTLs;i++){
            eqtl = eqtlInfoVec[i];
            // if(smrBesdBool) eqtl->eqtl_n = sampleSize;
            // eqtlIndxMap.insert(pair<int,string> (eqtl->index,eqtl->rsID));
            // if(eqtl->rsID == "rs547339847") {
            //     cout << "here " << i << endl;
            // }
            eqtlIndxMap.insert(pair<int,string> (i,eqtl->rsID));
        }
        //////////////////////////////////////////////
        totalSteps = _cols.size();
        auto startTime = std::chrono::steady_clock::now();
        //////////////////////////////////////////////
        // Step 2. select subset of eqtl info based one kept eqtl and gene list
        int size = 2; 
        set<string> cisSnpIDSetInGeneTmp;
        for(unsigned j = 0, geneIdx = 0; j < _cols.size() - 1; j += size ){
            gene = geneInfoVec[j/2];
            cisSnpIDSetInGeneTmp = gene->cisSnpIDSetInGene;
            gene->gene2CisSnpVec.clear();
            gene->cisSnpIDSetInGene.clear();
            gene->cisSnpNameVec.clear();
            // gene->eQTLPvalVec.clear();
            // gene->eQTLChisqVec.clear();
            gene->numSnpInGene = 0;
            if(! gene->kept){
                skip += size;
                continue;
            }
            // if(skip) j = j + skip;
            skip = 0;
            int numEqtlInGene = _cols[j + 1] - _cols[j];

            VectorXd eqtlEffect(numEqtlInGene), eqtlEffectSE(numEqtlInGene);
            // VectorXd eqtlSamSize(numEqtlInGene);
            vector<int> eqtlEffectIndxInGene;
            eqtlEffectIndxInGene.clear();
            // remove eqtls not pass QC
            int localIdx,icldIdx;
            eqtlEffect.setZero();
            eqtlEffectSE.setZero();
            // bool haveSamSize = false;
            for(unsigned i  = _cols[j], localIdx  = 0,icldIdx = 0; i < _cols[j + 1] ; i ++){
                // if(eqtlIndxMap.find(_rowid[i]) == "rs547339847"){
                //     cout << "here " << endl;
                // }
                if(eqtlIndxMap.find(_rowid[i]) != eqtlIndxMap.end()){
                    iterEqtl = eqtlInfoMap.find(eqtlIndxMap.at(_rowid[i]));
                    eqtl = iterEqtl->second;
                    // if(eqtl->rsID == "rs547339847") {
                    //     cout << "here" << endl;
                    // }
                    if(!eqtl->included){localIdx ++; continue;}
                    // check gene, used for cis-wind filter
                    if(cisSnpIDSetInGeneTmp.size() !=0){
                        if(cisSnpIDSetInGeneTmp.find(eqtl->rsID) == cisSnpIDSetInGeneTmp.end()){localIdx ++; continue;}
                    }
                    if(hasSnpInfo){ // consider match issue between eqtl and snp
                        iterSnp = snpInfoMap.find(eqtlIndxMap.at(_rowid[i]));
                        if(iterSnp == snpInfoMap.end()){ localIdx ++; continue;}
                    }
                    // beta and se value
                    eqtlEffect(localIdx) = _val[i];
                    eqtlEffectSE(localIdx) = _val[i + numEqtlInGene];
                    /// previously, both eqtl effect and se will be checked, but now I think it's appropriate to check 
                    /// se only
                    // if(( fabs(eqtlEffect(localIdx) + 9) < 1e-6 ) || 
                    //    ( fabs(eqtlEffectSE(localIdx) + 9) < 1e-6 ) || 
                    //      fabs(eqtlEffect(localIdx)      ) < 1e-6 ) {localIdx ++; continue;}
                    if(fabs(eqtlEffectSE(localIdx) + 9) < 1e-6 ) {localIdx ++; continue;}
                    eqtlEffectIndxInGene.push_back(localIdx);
                    gene->cisSnpNameVec.push_back(eqtl->rsID);
                    gene->cisSnpIDSetInGene.insert(eqtl->rsID);
                    gene->gene2CisSnpVec.push_back(_rowid[i]);
                    gene->cisSnpID2IdxMapInGene.insert(pair<string,int>(eqtl->rsID,icldIdx));
                    // gene->eQTLChisqVec.push_back(chisqBuf);
                    // gene->eQTLPvalVec.push_back(pValueBuf);
                    gene->hasEqtl = true;
                    eqtl->isInGene = true;
                    icldIdx ++;
                } //end of if(eqtlIndxMap.find(_rowid[i]) != eqtlIndxMap.end()){
                localIdx ++;
            } // end of loop beta and se
            // if(haveSamSize) {gene->sampleSize = Gadget::findMedian(eqtlSamSize);
            // } else {
            //     gene->sampleSize = estimateSamSize(eqtlEffect,eqtlEffectSE);
            // }
            // cout << "beta: " << eqtlEffect << " se: " << eqtlEffectSE << endl;
            // cout << "estimated sample size for gene " << gene->ensemblID << " is  " << estimateSamSize(eqtlEffect,eqtlEffectSE) << endl;
            gene->numSnpInGene = eqtlEffectIndxInGene.size();
            if(!makeGeneLDEigenBool){
                gene->eQTLMarginEffect = eqtlEffect(eqtlEffectIndxInGene);//.cast<double>();
                gene->eQTLMarginEffectSE = eqtlEffectSE(eqtlEffectIndxInGene);//.cast<double>();
            }
            ///////////////////////////////////////////////////////////////////////
            // float currentStep = j;
            // float progress = (currentStep * 100.0) / totalSteps;
            // LOGGER.p(0, " " + to_string(static_cast<float>(progress)) + "% ","Read eqtl info from sparse besd format:" );
            Gadget::showProgressBar(j, totalSteps, startTime, "Read eqtl info from sparse besd format");
            ///////////////////////////////////////////////////////////////////////
            // remove genes if there is no eqtl in gene 
            if (gene->gene2CisSnpVec.size() == 0){ gene->kept = false;}
            if (++geneIdx == numGenes) break;
        }        
        // terminate
        _cols.clear();
        _rowid.clear();
        _val.clear(); 
        free (buffer);

    }else {
        LOGGER.e(0,"BayesOmics doesn't support this format. Please online manual for details.");
    }
    
    // remove eqtls that do not belong to any gene
    for(unsigned i =0; i < numeQTLs; i++){
        eqtl = eqtlInfoVec[i];
        if(!eqtl->isInGene) {
            eqtl->included = false;
            // cout << "remove eQTL that do no belong gen: " << eqtl->rsID << endl;
        }
    }
    // remove genes that do not have eQTLs
    for(unsigned i =0; i < numGenes ; i++){
        gene = geneInfoVec[i];
        if(!gene->hasEqtl) {
            gene->kept = false;
        }
    }

    includeMatchedEqtl(); // 
    keptGeneInfoVec  = makeKeptGeneInfoVec(geneInfoVec);
    numKeptGenes = (unsigned) keptGeneInfoVec.size();

    LOGGER<<"eQTL summary data of "<< numKeptGenes <<" Probes and "<< numIncdEqtls <<" SNPs to be included from [" + besdFile + "]." << endl;
}

void Data::readQeuryGZFormat(const string &eqtlSummaryQueryFile,vector<GeneInfo *> &geneInfoVecLocal,map<string, GeneInfo *> &geneInfoMapLocal,
                    vector<EqtlInfo *> &eqtlInfoVecLocal,map<string, EqtlInfo *> &eqtlInfoMapLocal, 
                    bool &smrBesdBool,const bool hasSnpInfo,const bool hasGeneLdmInfo){
    std::ifstream file(eqtlSummaryQueryFile, std::ios_base::in | std::ios_base::binary);
    
    if (!file) LOGGER.e(0, " can not open the file [" + eqtlSummaryQueryFile + "] to read.");
    LOGGER << "Reading xQTL information from [" + eqtlSummaryQueryFile + "]." << endl;
    // calculate the totoal rows 
    int totalLines = Gadget::getTotalLineNumFromGzFile(eqtlSummaryQueryFile);
    string header;
    boost::iostreams::filtering_istream in;
    in.push(boost::iostreams::gzip_decompressor());
    in.push(file);
    getline(in, header);
    Gadget::Tokenizer colData;string inputStr;string sep(" \t\n");
    colData.getTokens(header, sep); 
    if(colData.size() != 14){
        LOGGER.e(0,"Wrong format. there should be 14 columns: GeneID\tGeneChr\tGeneLength\tGenePhyPos\tSNPID\tSNPChr\tSNPPhyPos\tA1\tA2\tA1Freq\tBETA\tSE\tPvalue\tN");
    }
    // check each columns
    std::vector<std::string> headderVec = {"GeneID","GeneChr","GeneLength","GenePhyPos","SNPID","SNPChr","SNPPhyPos","A1","A2","A1Freq","BETA","SE","Pvalue","N"};
    for(int i = 0; i < 14; i++){
        if(colData[i] != headderVec[i]){
            LOGGER.e(0,"The header of " + to_string(i)+ "-th column is [" + colData[i] + "], but it should be [" + headderVec[i] + "]");
        }
    }
    /// QC process
    unsigned line=0, match=0;
    unsigned numInconAllele=0, numInconAf=0, numFixed=0, numMafMin=0, numMafMax=0, numOutlierN=0;
    unsigned numPvalPruned=0,numInconsistent2GWAS = 0;
    unsigned numFlip=0;
    bool inconAllele, inconAf, fixed, ismafmin, ismafmax, isPvalPruned;
    //////
    GeneInfo *gene;
    EqtlInfo *eqtl;
    SnpInfo * snp;
    map<string, GeneInfo *>::iterator iterGene;
    map<string, EqtlInfo *>::iterator iterEqtl;
    map<string, SnpInfo*>::iterator iterSnp;
    unsigned geneChr,snpChr, sampleSize;
    int snpPhyPos,geneLength;
    string allele1,allele2;
    double a1freq,beta,se,genePhyPos;
    boost::multiprecision::cpp_dec_float_50 pvalue;
    string  geneID,lastGeneID,snpID;
    double lastGenePos =0, lastSnpPos = 0;
    vector<double>  eQTLMarginEffect;   // eQTL marginal effect
    vector<double> eQTLMarginEffectSE; // eQTL standard error of marginal effect
    set<string> flipSet,removeSNPset;
    unsigned geneIdx = 0, eqtlIdx = 0,eqtlInGene = 0;
    // header: GeneID   GeneChr GeneStart       GeneEnd GenePhyPos      SNPID   SNPChr  PhyPos  A1      A2      A1Freq  BETA    SE      Chisq   LOG10P  N
    auto startTime = std::chrono::steady_clock::now();
    while (getline(in,inputStr)) {
        Gadget::showProgressBar(line, totalLines, startTime,"Read xQTL info from query.gz format");
        if (inputStr.empty()) {continue;}
        colData.getTokens(inputStr, sep);
        geneID = colData[0];
        geneChr = boost::lexical_cast<int>(colData[1]);
        geneLength = boost::lexical_cast<int>(colData[2]);
        genePhyPos = boost::lexical_cast<double>(colData[3]);
        snpID = colData[4];
        snpChr = boost::lexical_cast<int>(colData[5]);
        snpPhyPos = boost::lexical_cast<int>(colData[6]);
        allele1 = colData[7];
        allele2 = colData[8];
        a1freq = boost::lexical_cast<double>(colData[9]);
        beta = boost::lexical_cast<double>(colData[10]);
        se = boost::lexical_cast<double>(colData[11]);

        pvalue = boost::lexical_cast<boost::multiprecision::cpp_dec_float_50>(colData[12]);
        if(!Gadget::isPositiveInteger(colData[13])){
            LOGGER.w(0,"The type of per SNP sample size [" + colData[13] + "] is wrong. Remove this eQTL." );
            continue;
        }
        sampleSize  = boost::lexical_cast<int>(colData[13]);
        //////////////////////////////////////////////////////////
        /// remove snps 
        if(hasSnpInfo) {
            // have GWAS snplist
            if(snpInfoMap.find(snpID) == snpInfoMap.end()) {
                // numInconsistent2GWAS++;
                removeSNPset.insert(snpID);
                continue;
            }
        }
        /// check if gene exist, if not, create one,
        /// or 
        if(geneInfoMapLocal.find(geneID) == geneInfoMapLocal.end()){
            // if the gene is the firt skip storage directly,
            // if not, we need to save them into 
            if(geneInfoMapLocal.size() != 0){
                // now we need to save beta and se into previous last gene
                gene->eQTLMarginEffect = Eigen::Map<VectorXd>(&eQTLMarginEffect[0], eQTLMarginEffect.size());
                gene->eQTLMarginEffectSE = Eigen::Map<VectorXd>(&eQTLMarginEffectSE[0], eQTLMarginEffectSE.size());
                gene->sampleSize = Gadget::calculateMean(gene->cisSnpSampleSizeMap);
                gene->sampleSizeSD = Gadget::calculateStandardDeviation(gene->cisSnpSampleSizeMap,gene->sampleSize);
                eQTLMarginEffect.clear();
                eQTLMarginEffectSE.clear();
                if(gene->eQTLMarginEffect.size() > 1) gene->hasEqtl = true;
            }
            // Not found gene, we need to create a new one
            gene = new GeneInfo(geneIdx++, geneID, geneChr);
            geneInfoVecLocal.push_back(gene);
            if (geneInfoMapLocal.insert(pair<string, GeneInfo*>(gene->ensemblID, gene)).second == false) {
                LOGGER.e(0, " Duplicate  gene ID found: \"" + gene->ensemblID + "\". All snps within the same gene should be put together.");
            }
            // Not found eqtl, we need to create a new one
            if(eqtlInfoMapLocal.find(snpID) == eqtlInfoMapLocal.end()){
                eqtl = new EqtlInfo(eqtlIdx++, snpID,allele1,allele2,snpChr, 0, snpPhyPos,a1freq);
                eqtlInfoVecLocal.push_back(eqtl);
                eqtlInfoMapLocal.insert(pair<string, EqtlInfo*>(eqtl->rsID,eqtl ));
            }
            eqtlInGene = 0;
            eqtl = eqtlInfoMapLocal.find(snpID)->second;
            gene->geneLength = geneLength;
            gene->midPhyPos = genePhyPos;
            gene->start = genePhyPos - geneLength * 0.5;
            gene->end = genePhyPos + geneLength * 0.5;
            gene->cisSnpNameVec.push_back(eqtl->rsID);
            gene->cisSnpIDSetInGene.insert(eqtl->rsID);
            gene->cisSnpID2IdxMapInGene.insert(pair<string,int>(eqtl->rsID,eqtlInGene));
            gene->cisSnpSampleSizeMap.insert(pair<string,int>(eqtl->rsID,sampleSize));
            // gene->sampleSize = 
            eQTLMarginEffect.push_back(beta);
            eQTLMarginEffectSE.push_back(se);
            lastGeneID = geneID;
            eqtl->isInGene = true;
            /////////////////////////////////////
            /// check consistency of allele
            if(hasSnpInfo) {
                iterSnp = snpInfoMap.find(eqtl->rsID);
                if(iterSnp != snpInfoMap.end()){
                    snp = iterSnp->second;
                    if (allele1 == snp->a1 && allele2 == snp->a2) {
                    } else if (allele1 == snp->a2 && allele2 == snp->a1) {
                        eqtl->flipped = true;
                        // ++numFlip;
                        flipSet.insert(eqtl->rsID);
                    } else {
                        inconAllele = true;
                        ++numInconAllele;
                        eqtl->included = false;
                        LOGGER.w(0,"eQTLs "+ eqtl->rsID + " has different allele (A1: " + allele1 + ",A2:" + allele2 + "), however in the GWAs LD reference (A1:" + snp->a1 + ",A2:" + snp->a2 + ")."  );
                    }
                } // end of iterSnp.
            }
            
            ////////////////////////////////////////
        } else {
            // Found gene, now add SNPs to the gene
            gene = geneInfoMapLocal.find(geneID)->second;
            // Not found eqtl, create one
            if(eqtlInfoMapLocal.find(snpID) == eqtlInfoMapLocal.end()){
                eqtl = new EqtlInfo(eqtlIdx++, snpID,allele1,allele2,snpChr, 0, snpPhyPos,a1freq);
                eqtlInfoVecLocal.push_back(eqtl);
                eqtlInfoMapLocal.insert(pair<string, EqtlInfo*>(eqtl->rsID,eqtl ));
            }
            // Found eqtl
            eqtlInGene++;
            eqtl = eqtlInfoMapLocal.find(snpID)->second;
            gene->start = genePhyPos - geneLength * 0.5;
            gene->end = genePhyPos + geneLength * 0.5;
            gene->midPhyPos = genePhyPos;
            gene->cisSnpNameVec.push_back(eqtl->rsID);
            gene->cisSnpIDSetInGene.insert(eqtl->rsID);
            gene->cisSnpID2IdxMapInGene.insert(pair<string,int>(eqtl->rsID,eqtlInGene));
            gene->cisSnpSampleSizeMap.insert(pair<string,int>(eqtl->rsID,sampleSize));
            eQTLMarginEffect.push_back(beta);
            eQTLMarginEffectSE.push_back(se);
            lastGeneID = geneID;
            eqtl->isInGene = true;
            /////////////////////////////////////
            /// check consistency of allele
            if(hasSnpInfo) {
                iterSnp = snpInfoMap.find(eqtl->rsID);
                if(iterSnp != snpInfoMap.end()){
                    snp = iterSnp->second;
                    if (allele1 == snp->a1 && allele2 == snp->a2) {
                    } else if (allele1 == snp->a2 && allele2 == snp->a1) {
                        eqtl->flipped = true;
                        flipSet.insert(eqtl->rsID);
                    } else {
                        LOGGER.e(0,"eQTLs "+ eqtl->rsID + " has different allele (A1: " + allele1 + ",A2:" + allele2 + "), however in the GWAs LD reference (A1:" + snp->a1 + ",A2:" + snp->a2 + ")."  );
                        inconAllele = true;
                        ++numInconAllele;
                        eqtl->included = false;
                    }
                } // end of iterSnp.
            }
            ////////////////////////////////////////
        }
        line++;
    } // end of while
    // for the last gene after ending the loop
    if(geneInfoMapLocal.find(lastGeneID) != geneInfoMapLocal.end()){
        gene = geneInfoMapLocal.find(lastGeneID)->second;
        gene->eQTLMarginEffect = Eigen::Map<VectorXd>(&eQTLMarginEffect[0], eQTLMarginEffect.size());
        gene->eQTLMarginEffectSE = Eigen::Map<VectorXd>(&eQTLMarginEffectSE[0], eQTLMarginEffectSE.size());
        gene->sampleSize = Gadget::calculateMean(gene->cisSnpSampleSizeMap);
        gene->sampleSizeSD = Gadget::calculateStandardDeviation(gene->cisSnpSampleSizeMap,gene->sampleSize);
        if(gene->eQTLMarginEffect.size() > 1) gene->hasEqtl = true;
    }
    // summary data quality
    if(hasSnpInfo) {
        numEqtlFlip = flipSet.size();
        numInconsistent2GWAS = removeSNPset.size();
        if(numInconsistent2GWAS) LOGGER << numInconsistent2GWAS << " SNPs that cannot be found in GWAS LD reference are removed." << endl;
        if (numEqtlFlip) LOGGER << "flipped " << numEqtlFlip << " SNPs according to the minor allele in the reference and GWAS samples." << endl;
        if(numInconAllele) LOGGER << "" << numInconAllele << " inconsistent SNPs will be removed based on GWAS LD reference and samples." << endl;
    }

    // Close the file
    file.close();
    // Now we need to summary gene and eqtl information
    LOGGER << "" << geneInfoVecLocal.size() << " genes with " << eqtlInfoVecLocal.size() <<  " eqtls to be included from [" + eqtlSummaryQueryFile + "]." << endl;
}

void Data::saveBesdFormat(const string title,const bool makeBesdSmrBool){
    //////////////////////////////////////////////
    // update gene and eQTL info;
    //////////////////////////////////////////////
    incdEqtlInfoVec = makeIncdEqtlInfoVec(eqtlInfoVec);
    numIncdEqtls = (unsigned) incdEqtlInfoVec.size();
    keptGeneInfoVec  = makeKeptGeneInfoVec(geneInfoVec);
    numKeptGenes = (unsigned) keptGeneInfoVec.size();

    if(numIncdEqtls == 0) LOGGER.e(0," No SNP is retained in saving besd data process.");
    if(numKeptGenes == 0) LOGGER.e(0," No gene is retained in saving besd data process.");

    LOGGER<<"\nSaving xQTL summary statistics data to besd format..."<<endl;
    //////////////////////////////////////////////
    // save esi file format
    //////////////////////////////////////////////
    string esdfile = title + ".esi";
    ofstream smr(esdfile.c_str());
    if (!smr) LOGGER.e(0,"Can not open the ESI file " + esdfile + " to save!");
    for (int i = 0;i < numIncdEqtls; i++) {
        EqtlInfo *eqtl = incdEqtlInfoVec[i];
        // if(eqtl->rsID == "rs201540471") {
        //     cout <<endl ;
        // }
        if(eqtl->afVec.size() > 1){
            Eigen::VectorXd afVec = Eigen::Map<const Eigen::VectorXd>(eqtl->afVec.data(), eqtl->afVec.size());
            Eigen::VectorXd samVec = Eigen::Map<const Eigen::VectorXd>(eqtl->samVec.data(), eqtl->samVec.size());
            eqtl->af = (afVec.array() * samVec.array()).sum()/samVec.sum();
        }
        smr <<  eqtl->chrom
        <<'\t'<< eqtl->rsID
        <<'\t'<< eqtl->genPos 
        <<'\t'<< eqtl->physPos 
        <<'\t'<< eqtl->a1
        <<'\t'<< eqtl->a2
        <<'\t'<< eqtl->af;
        // if(!makeBesdSmrBool) smr <<'\t'<< eqtl->eqtl_n;
        smr << endl;
    }
    smr.close();
    LOGGER << numIncdEqtls <<" SNPs have been saved in the file [" + esdfile + "]."<<endl;
    //////////////////////////////////////////////
    // save epi file format
    //////////////////////////////////////////////
    esdfile = title + ".epi";
    smr.open(esdfile.c_str());
    if (!smr) LOGGER.e(0," can not open the EPI file " + esdfile + " to save!");

    for (int i = 0;i < numKeptGenes; i++) {
        GeneInfo * gene = keptGeneInfoVec[i];
        smr<< gene->chrom
        <<'\t'<< gene->ensemblID
        <<'\t'<< gene->start
        << std::fixed << std::setprecision(1)
        <<'\t'<< gene->midPhyPos
        <<'\t'<< gene->ensemblID
        << '\t' << gene->sampleSize
        <<'\t'<< gene->geneOri
        << endl;
    }
    smr.close();
    LOGGER << numKeptGenes <<" probes have been saved in the file [" + esdfile + "]."<<endl;
    //////////////////////////////////////////////
    // save besd CSC format
    //////////////////////////////////////////////
    // sparse CSC format, we need construct column index pointers, row indices and data values
    vector<uint64_t> _cols; // column index pointers
    vector<uint32_t> _rowid; // row indices
    vector<float> _val; // data values

    EqtlInfo *eqtl;
    GeneInfo *gene;
    map<string, EqtlInfo*>::iterator iterEqtl;
    vector<float> eqtlBetaVec, eqtlSeVec,eqtlIdxVec;
    int eqtlIdx;
    uint64_t cumsumCols = 0;
    _cols.push_back(cumsumCols);
    // several necessary variables
    // gene->cisSnpIDSetInGene
    // gene->numSnpInGene
    // gene->cisSnpID2IdxMapInGene
    // gene->eQTLMarginEffect
    // gene->eQTLMarginEffectSE
    uint64_t numSnpInGeneCurr = 0;
    auto startTime = std::chrono::steady_clock::now();
    for (int i = 0;i < numKeptGenes; i++) {
        Gadget::showProgressBar(i, numKeptGenes, startTime, "Construct the Column-Wise Sparse Matrix ");
        gene = keptGeneInfoVec[i];
        if(!gene->kept) continue;
        eqtlBetaVec.clear();
        eqtlSeVec.clear();
        eqtlIdxVec.clear();
        numSnpInGeneCurr = 0;
        if(gene->cisSnpNameVec.size()!= gene->cisSnpIDSetInGene.size()){
                LOGGER.e(0,"Error, the size of cisSnpNameVec should be equal to size of cisSnpIDSetInGene for gene: [ " +gene->ensemblID + " ].");
        }
        // for(int j = 0; j < numIncdEqtls; j ++){
        // eqtl = incdEqtlInfoVec[j];
        for(int j =0; j < gene->cisSnpNameVec.size(); j ++){
            iterEqtl = eqtlInfoMap.find(gene->cisSnpNameVec[j]);
            eqtl = iterEqtl->second;
            if(!eqtl->included) continue;
            if(gene->cisSnpIDSetInGene.find(eqtl->rsID) != gene->cisSnpIDSetInGene.end()){
                if(gene->cisSnpID2IdxMapInGene.find(eqtl->rsID) != gene->cisSnpID2IdxMapInGene.end()){
                    eqtlIdx = gene->cisSnpID2IdxMapInGene[eqtl->rsID];
                } else {
                    eqtl->included = false;
                    LOGGER.e(0,"Could not find eqtl " +eqtl->rsID + " in gene " + gene->ensemblID + ".");
                }
                eqtlBetaVec.push_back( (float) gene->eQTLMarginEffect(eqtlIdx));
                eqtlSeVec.push_back( (float) gene->eQTLMarginEffectSE(eqtlIdx));
                eqtlIdxVec.push_back(eqtl->index);
                numSnpInGeneCurr ++;
            }
        }
        if(numSnpInGeneCurr != 0){
            // gene->cisSnpNameVec.clear();
            // gene->eQTLMarginEffect.resize(0); // remove this for save memeroy
            // gene->eQTLMarginEffectSE.resize(0); // remove this for save memeory
            // add beta effect value
            cumsumCols += numSnpInGeneCurr; 
            _cols.push_back(cumsumCols);
            _rowid.insert(_rowid.end(),eqtlIdxVec.begin(),eqtlIdxVec.end());
            _val.insert(_val.end(),eqtlBetaVec.begin(),eqtlBetaVec.end());
            // add se value
            cumsumCols += numSnpInGeneCurr;
            _cols.push_back(cumsumCols);
            _rowid.insert(_rowid.end(),eqtlIdxVec.begin(),eqtlIdxVec.end());
            _val.insert(_val.end(),eqtlSeVec.begin(),eqtlSeVec.end());
        } else {
            LOGGER.e(0,"There is no eqtls in gene " + gene->ensemblID + ".");
        }       
    }

    // // // debug
    // cout << "cols: ";
    // for(unsigned k = 0; k < _cols.size(); k++) cout << _cols[k] << " ";
    // cout << endl;
    // cout << "_rowid: ";
    // for (unsigned k = 0; k < _rowid.size();k ++) cout << _rowid[k] << " ";
    // cout << endl;
    // cout << "_val: ";
    // for (unsigned k = 0; k < _val.size();k ++) cout << _val[k] << " ";
    // cout << endl;

    esdfile = title + ".besd";
    FILE * smrbesd;
    smrbesd = fopen (esdfile.c_str(), "wb");
    if(false)
    {
        // uint64_t bsize=(eqtlinfo->_include.size()*eqtlinfo->_snpNum<<1)+1;
        //     float* buffer=(float*)malloc (sizeof(float)*bsize);
        //     memset(buffer,0,sizeof(float)*bsize);
        //     float* ptr=buffer;
        //     *ptr++=0.0;
        //     uint64_t pro_num=eqtlinfo->_include.size();
        //     uint64_t snp_num=eqtlinfo->_snpNum;
        //     for(int i=0;i<pro_num;i++)
        //     {
        //         memcpy(ptr+(i<<1)*snp_num,&eqtlinfo->_bxz[eqtlinfo->_include[i]][0],sizeof(float)*snp_num);
        //         memcpy(ptr+((i<<1)+1)*snp_num,&eqtlinfo->_sexz[eqtlinfo->_include[i]][0],sizeof(float)*snp_num);
        //     }
        //     fwrite(buffer,sizeof(float), bsize, smrbesd);
        //     free(buffer);
    }
    else
    {
        uint64_t colSize=sizeof(uint64_t)*((numKeptGenes<<1)+1);  // num of probe
        uint64_t rowSize=sizeof(uint32_t)* _val.size();  // total number of non-zero snps
        uint64_t valSize=sizeof(float)* _val.size(); //
        uint64_t valNum= _val.size();
        uint64_t bufsize=sizeof(float)+sizeof(uint64_t)+colSize+rowSize+valSize;
        // cout << "bufsize: " << bufsize << endl;
        // memset function is to set a block of memory to a specific value, typically used for initializing arrays or structures;
        // memcpy is a function in C and C++ used for copying a block of memory from one location to another.
        char* buffer=(char*)malloc (sizeof(char)*bufsize);  //  allocates a block of memory large enough to store bufsize characters and assigns the address of this memory block to buffer. 
        memset(buffer,0,sizeof(char)*bufsize);  // sets the entire memory block pointed to by buffer to zero. This is often done to initialize or clear a buffer before it is used, ensuring that it doesn't contain any garbage values left over from previous use or from the memory allocation process.
        uint32_t ftype=SPARSE_FILE_TYPE_3F;
        memcpy(buffer,&ftype,sizeof(uint32_t)); // Copy 4 bytes of data from the location where ftype is stored into the memory location pointed to by buffer.
        char* wptr=buffer+sizeof(float); // wptr is now a pointer to a location in memory that is sizeof(float) bytes ahead of where buffer points to.
        memcpy(wptr,&valNum,sizeof(uint64_t)); // copies 8 bytes from the location of valNum to the memory location pointed to by wptr.
        wptr+=sizeof(uint64_t); // increments wptr by 8 bytes, effectively moving it forward in memory by the size of a uint64_t.
        uint64_t* uptr=(uint64_t*)wptr; // This line casts wptr to a pointer of type uint64_t*. It is assigning the modified wptr address to a new pointer uptr, which is specifically a pointer to uint64_t.
        *uptr++=0; // 

        uint64_t colNum=(numKeptGenes<<1)+1;// = (2 * numKeptGenes) + 1. shifts the bits of numKeptGenes one position to the left. Bitwise shifting to the left by one position is equivalent to multiplying the number by 2. So, if numKeptGenes is n, after this operation, it effectively becomes 2 * n.
        auto startTime = std::chrono::steady_clock::now();
        for(int i=0;i< numKeptGenes;i++)
        {
            Gadget::showProgressBar(i, numKeptGenes, startTime, "Save the Column-Wise Sparse Matrix");
            gene = keptGeneInfoVec[i];
            *uptr++= _cols[(gene->index<<1)+1];
            *uptr++= _cols[ (gene->index+1)<<1];
        }
        // https://stackoverflow.com/questions/9258432/why-does-memcpy-fail-to-copy-eigen-matrix-data-but-stdcopy-succeed 
        wptr+=colSize;
        memcpy(wptr,& _rowid[0],rowSize);
        wptr+=rowSize;
        memcpy(wptr,& _val[0],valSize);
        fwrite (buffer,sizeof(char), bufsize, smrbesd);
        free(buffer);
    }
    fclose (smrbesd);

    // LOGGER<<"Effect sizes (beta) and SE for "<<eqtlinfo->_include.size()<<" Probes and "<<eqtlinfo->_snpNum<<" SNPs have been saved in the binary file [" + esdfile + "]." <<endl;
    LOGGER<<"Effect sizes (beta) and SE for "<< numKeptGenes <<" Probes and "<< numIncdEqtls <<" SNPs have been saved in the binary file [" + esdfile + "]." <<endl;
}

void Data::saveAnnoPlainMatFormat(const string title, const bool isBinary,const bool hasSnpInfo ){
    string esdfile;
    if (isBinary){
        esdfile = title + ".anno.bin.txt";
    } else {
        esdfile = title + ".anno.con.txt";
    }
    ofstream anno(esdfile.c_str());
    if (!anno) LOGGER.e(0," Can not open the anno file [" + esdfile + "] to save!");
    EqtlInfo * eqtl;
    SnpInfo *snp;
    GeneInfo * gene;
    map<string, EqtlInfo* >::iterator iterEqtl;
    int eqtlIdx;
    double sEqtlValue;

    anno << "SNP" << "\t" << "Intercept" << "\t";
    if(hasSnpInfo){
        anno << "Anno" << endl;
        for(unsigned i = 0; i < numIncdSnps; i ++){
            snp= incdSnpInfoVec[i];
            sEqtlValue = 0;
            anno << snp->rsID << "\t" << "1" << "\t";
            iterEqtl= eqtlInfoMap.find(snp->rsID);
            if(iterEqtl!=eqtlInfoMap.end()){
                if(iterEqtl->second->included){sEqtlValue = 1;};
            }
            // if(!snp->iseQTL) sEqtlValue = 0;
            anno << sEqtlValue << endl;
        }
    } else {
        for(unsigned j = 0; j < numKeptGenes; j++){
            gene = keptGeneInfoVec[j];
            anno << gene->ensemblID << "\t";
        }
        anno << endl;
        for(unsigned i = 0; i < numIncdEqtls; i++){
            eqtl = incdEqtlInfoVec[i];
            anno << eqtl->rsID << "\t" << "1" << "\t";
            for(unsigned j = 0; j < numKeptGenes; j++){
                gene = keptGeneInfoVec[j];
                if(gene->cisSnpID2IdxMapInGene.find(eqtl->rsID) != gene->cisSnpID2IdxMapInGene.end()){
                    eqtlIdx = gene->cisSnpID2IdxMapInGene[eqtl->rsID];
                    sEqtlValue = gene->eQTLMarginEffect(eqtlIdx);
                    if(isBinary) sEqtlValue  = 1;
                } else {
                    sEqtlValue = 0;
                    // anno << "0" << "\t";
                } //
                anno << sEqtlValue << "\t";
            }
            anno << endl;
        }
    }
    anno.close();
    LOGGER<<""<< numKeptGenes <<" genes and "<< numIncdEqtls <<" cis-SNPs have been saved in the annotation file [" + esdfile + "]." <<endl;  
}

void Data::saveQueryMoQTLInfo(const string title){
    incdEqtlInfoVec = makeIncdEqtlInfoVec(eqtlInfoVec);
    numIncdEqtls = (unsigned) incdEqtlInfoVec.size();
    keptGeneInfoVec  = makeKeptGeneInfoVec(geneInfoVec);
    numKeptGenes = (unsigned) keptGeneInfoVec.size();

    if(numIncdEqtls == 0) LOGGER.e(0," No SNP is retained in saving besd data process.");
    if(numKeptGenes == 0) LOGGER.e(0," No gene is retained in saving besd data process.");
    LOGGER << "Save eQTL data into compressed query gz format..." << endl;

    string esdfile = title + ".query.gz";
    // Create an ofstream for output
    std::ofstream file(esdfile, std::ios_base::out | std::ios_base::binary);
    if (!file.is_open()) {LOGGER.e(0, " Can not open the file [" + esdfile + "] to save.");}
    // Create a filtering streambuf and attach the gzip compressor and the ofstream
    boost::iostreams::filtering_streambuf<boost::iostreams::output> outGZ;
    outGZ.push(boost::iostreams::gzip_compressor());
    outGZ.push(file);

    EqtlInfo * eqtl;
    GeneInfo * gene;
    map<string, EqtlInfo *>::iterator iterEqtl;
    int eqtlIdx;
    double sEqtlValue;

    std::stringstream dataStream;
    std::ostringstream genePhyPosRound,beta,se,chisq,pvalue,samSize;
    // string singleLine;
    // header line
    string header = "GeneID\tGeneChr\tGeneLength\tGenePhyPos\tSNPID\tSNPChr\tSNPPhyPos\tA1\tA2\tA1Freq\tBETA\tSE\tPvalue\tN\n";
    dataStream << header;
    double chiBuf, chisqBuf,pValueBuf,numNA = 0;
    auto startTime = std::chrono::steady_clock::now();
    for(unsigned i = 0; i < numKeptGenes; i++){
        Gadget::showProgressBar(i, numKeptGenes, startTime,"Save query.gz format");
        gene = keptGeneInfoVec[i];
        // use on decimal
        genePhyPosRound.str("");
        genePhyPosRound << std::fixed << std::setprecision(1) << gene->midPhyPos;
        if(gene->cisSnpNameVec.size()!= gene->cisSnpIDSetInGene.size()){
            LOGGER.e(0,"Error, the size of cisSnpNameVec should be equal to size of cisSnpIDSetInGene for gene: [ " +gene->ensemblID + " ].");
        }
        // for(int j = 0; j < numIncdEqtls; j ++){
        //     eqtl = incdEqtlInfoVec[j];
        for(int j =0; j < gene->cisSnpNameVec.size(); j ++){
            iterEqtl = eqtlInfoMap.find(gene->cisSnpNameVec[j]);
            eqtl = iterEqtl->second;
            if(gene->cisSnpID2IdxMapInGene.find(eqtl->rsID) != gene->cisSnpID2IdxMapInGene.end()){
                eqtlIdx = gene->cisSnpID2IdxMapInGene[eqtl->rsID];
            } else {
                eqtl->included = false;
            }
            if(!eqtl->included) continue;
            // cout << "eqtlIdx: " << eqtlIdx << endl;
            beta.str("");
            se.str("");
            chisq.str("");
            pvalue.str("");
            samSize.str("");
            // cout << gene->eQTLMarginEffect.size() << endl; 
            if(abs(gene->eQTLMarginEffect(eqtlIdx) ) < 1e-6 ){
                beta << std::scientific << std::setprecision(5) << gene->eQTLMarginEffect(eqtlIdx);
            } else{
                beta << std::fixed << std::setprecision(15) << gene->eQTLMarginEffect(eqtlIdx);
            }
            // se
            if(abs(gene->eQTLMarginEffectSE(eqtlIdx)) < 1e-6 ){
                se << std::scientific << std::setprecision(5) << gene->eQTLMarginEffectSE(eqtlIdx);
            } else{
                se << std::fixed << std::setprecision(15) << gene->eQTLMarginEffectSE(eqtlIdx);
            }
            if(abs(gene->eQTLMarginEffectSE(eqtlIdx)) > 1e-6){
                chiBuf = gene->eQTLMarginEffect(eqtlIdx)/gene->eQTLMarginEffectSE(eqtlIdx);
                chisqBuf = chiBuf * chiBuf;
                pValueBuf = Stat::ChiSq::pchisq(chisqBuf,1);
            }
            // chisq
            if(abs(chisqBuf) < 1e-6 ){
                chisq << std::scientific << std::setprecision(3) << chisqBuf;
            } else{
                chisq << std::fixed << std::setprecision(9) << chisqBuf;
            }
            // pvalue
            if(abs(pValueBuf) < 1e-6 ){
                if(pValueBuf == 0) pValueBuf = std::numeric_limits<double>::min();
                pvalue << std::scientific << std::setprecision(3) << (pValueBuf);
            } else{
                pvalue << std::fixed << std::setprecision(9) << (pValueBuf);
            }

            // sample size 
            samSize.str("");
            bool hasSamSize = false;
            // use the same gene sample size for all snp within the gene
            if(gene->sampleSize != -999){
                samSize.str("");
                samSize << std::fixed << std::setprecision(0) <<  (gene->sampleSize);
                hasSamSize = true;
            }
            // perSNP sample size will be changed if additional sample size provided.
            if(gene->cisSnpSampleSizeMap.find(eqtl->rsID) != gene->cisSnpSampleSizeMap.end()){
                samSize.str("");
                samSize << std::fixed << std::setprecision(0) <<  gene->cisSnpSampleSizeMap.at(eqtl->rsID);
                hasSamSize = true;
            } 
            // if non snp 
            if(!hasSamSize){
                samSize << "NA";
                numNA ++;
            }

            dataStream << gene->ensemblID + "\t" + to_string(gene->chrom) + "\t" + to_string(gene->geneLength) + "\t" + genePhyPosRound.str() + "\t" +
                        eqtl->rsID + "\t" + to_string(eqtl->chrom) + "\t" + to_string(eqtl->physPos) + "\t" + 
                        eqtl->a1 + "\t" + eqtl->a2 + "\t" + to_string(eqtl->af) + "\t" +
                        beta.str() + "\t" + se.str() + "\t" + 
                        pvalue.str() + "\t" + samSize.str() + 
                        "\n";           
        }
    }
    boost::iostreams::copy(dataStream, outGZ);
    // Close the filtering streambuf
    outGZ.pop();  // Pop the gzip compressor
    // Now close the file stream
    file.close();
    if(numNA > 0)LOGGER.w(0, to_string(numNA) + " eQTLs have NA per SNP sample size." );
    LOGGER<<""<< numKeptGenes <<" genes and "<< numIncdEqtls <<" cis-SNPs have been saved in the annotation file [" + esdfile + "]." <<endl;  
}

void Data::saveSumstatsFormat(const string title){
    string esdfile = title + ".sumstats";
    // Create an ofstream for output
    ofstream anno(esdfile.c_str());
    if (!anno) LOGGER.e(0," Can not open the anno file [" + esdfile + "] to save!");
   
    EqtlInfo * eqtl;
    GeneInfo * gene;
    map<string, EqtlInfo *>::iterator iterEqtl;
    int eqtlIdx;
    double sEqtlValue;

    std::stringstream dataStream;
    std::ostringstream genePhyPosRound,zscore,samSize;
    // string singleLine;
    // header line
    LOGGER << "Save eQTL data into ldsc .sumstats format..." << endl;
    string header = "GENE\tGENE_COORD\tSNP\tCHR\tSNP_COORD\tN\tZ\n";
    anno << header;
    double chiBuf, chisqBuf,pValueBuf,numNA = 0;
    auto startTime = std::chrono::steady_clock::now();
    for(unsigned i = 0; i < numKeptGenes; i++){
        Gadget::showProgressBar(i, numKeptGenes, startTime, "Save MESC moQTL summary statistics");
        gene = keptGeneInfoVec[i];
        // use on decimal
        genePhyPosRound.str("");
        genePhyPosRound << std::fixed << std::setprecision(0) << gene->midPhyPos;
        if(gene->cisSnpNameVec.size()!= gene->cisSnpIDSetInGene.size()){
            LOGGER.e(0,"Error, the size of cisSnpNameVec should be equal to size of cisSnpIDSetInGene for gene: [ " +gene->ensemblID + " ].");
        }
        for(int j =0; j < gene->cisSnpNameVec.size(); j ++){
            iterEqtl = eqtlInfoMap.find(gene->cisSnpNameVec[j]);
            eqtl = iterEqtl->second;
            if(gene->cisSnpID2IdxMapInGene.find(eqtl->rsID) != gene->cisSnpID2IdxMapInGene.end()){
                eqtlIdx = gene->cisSnpID2IdxMapInGene[eqtl->rsID];
            } else {
                eqtl->included = false;
            }
            if(!eqtl->included) continue;
            /// N
            // sample size 
            bool hasSamSize = false;
            // use the same gene sample size for all snp within the gene
            if(gene->sampleSize != -999){
                samSize.str("");
                samSize << std::fixed << std::setprecision(0) <<  (gene->sampleSize);
                hasSamSize = true;
            }
            // perSNP sample size will be changed if additional sample size provided.
            if(gene->cisSnpSampleSizeMap.find(eqtl->rsID) != gene->cisSnpSampleSizeMap.end()){
                samSize.str("");
                samSize << std::fixed << std::setprecision(0) <<  gene->cisSnpSampleSizeMap.at(eqtl->rsID);
                hasSamSize = true;
            } 
            // if non snp 
            if(!hasSamSize){
                samSize << "NA";
                numNA++;
            }
            ///// zscore
            zscore.str("");
            if(abs(gene->eQTLMarginEffectSE(eqtlIdx)) > 1e-6 ){
                // cout << "ef: " << gene->eQTLMarginEffect(eqtlIdx) << " se: " << gene->eQTLMarginEffectSE(eqtlIdx) << endl;
                zscore << gene->eQTLMarginEffect(eqtlIdx)/gene->eQTLMarginEffectSE(eqtlIdx);
            }
            
            anno << gene->ensemblID + "\t" + genePhyPosRound.str() + "\t" +
                        eqtl->rsID + "\t" + to_string(eqtl->chrom) + "\t" + to_string(eqtl->physPos) + "\t" + 
                        samSize.str() + "\t" + zscore.str() + "\n";
        }
    }
    // Now close the file stream
    anno.close();
    if(numNA > 0)LOGGER.w(0, to_string(numNA) + " eQTLs have NA per SNP sample size." );
    LOGGER<<""<< numKeptGenes <<" genes and "<< numIncdEqtls <<" cis-SNPs have been saved in the annotation file [" + esdfile + "]." <<endl;  
}

void Data::reSampleXqtlEffect(const int reSamType,const string &besdFile,const string &geneEigenMatrixFile,const int aimN, const double eigenCutoff,const double diag_mod){        
    // read xQTL LD info
    readEigenMatrixGene(geneEigenMatrixFile,eigenCutoff,"gene",0);
    vector<GeneInfo *> geneInfoVecBESD;
    map<string, GeneInfo *> geneInfoMapBESD;
    vector<EqtlInfo *> eqtlInfoVecBESD;
    map<string, EqtlInfo *> eqtlInfoMapBESD;
    bool epiBool, esiBool, besdBool;
    bool hasSnpInfo = true;
    bool hasGeneLdmInfo = true;
    bool smrBesdBool = false;
    // Step 4. read BESD format xQTL info to extract b and se
    epiBool = readMultiEpiFile(besdFile + ".epi", geneInfoVecBESD,geneInfoMapBESD,hasSnpInfo,hasGeneLdmInfo);
    esiBool = readMultiEsiFile(besdFile + ".esi",eqtlInfoVecBESD, eqtlInfoMapBESD,smrBesdBool,hasSnpInfo,hasGeneLdmInfo);
    if(epiBool && esiBool)
        besdBool = readMultiBesdFile(besdFile + ".besd",geneInfoVecBESD,geneInfoMapBESD,eqtlInfoVecBESD, eqtlInfoMapBESD,smrBesdBool,hasSnpInfo,hasGeneLdmInfo);

    // Summary genes
    alignxQTLGeneAndBESDGeneInfo(geneInfoVecBESD,geneInfoMapBESD,eqtlInfoVecBESD, eqtlInfoMapBESD);

    // Step 5. Now we have gene info from xQTL LD and gene info from BESD. we need to align them and then impute missing xQTL effect if needed.
    LOGGER << "Begin to re-sample xqtl effect sample size.." << endl;
    if(reSamType == 0)
            LOGGER << "Re-sample xQTL beta effect and its se based on imputation method (Zheng, Zhili, et al (2024))." << endl;
    if(reSamType == 1)
            LOGGER << "Re-sample xQTL beta effect and its se based on cross-trait LDSC method (Bulik-Sullivan, Brendan, et al. (2015))." << endl;

    for(unsigned j = 0; j < numKeptGenes; j++){
        GeneInfo * gene = keptGeneInfoVec[j];
        if(gene->impCisSnpEffBool){
            LOGGER.e(0,"Please check wether snplists from xQTL LD and BESD are consistent. If not, please check regenerate based on common snplists from GWAS and BESD data.");
        }
        ////// now 
        VectorXd snp2pqTmp;         // 2pq of SNPs
        VectorXd bTmp;              // beta from eQTL summary data
        VectorXd seTmp;             // se from eQTL summary data
        VectorXd nTmpEst;           // sample size for each SNP in eQTL
        VectorXd betaImp, seImp;
        int numIncdEqtlsInGene = gene->cisSnpNameVec.size();
        snp2pqTmp.resize(numIncdEqtlsInGene);
        bTmp.resize(numIncdEqtlsInGene);
        seTmp.resize(numIncdEqtlsInGene);
        nTmpEst.resize(numIncdEqtlsInGene);
        bool imputeSamBool = false;
        for (unsigned i = 0; i < numIncdEqtlsInGene; i++){
            EqtlInfo * eqtl = eqtlInfoMap[gene->cisSnpNameVec[i]];
            snp2pqTmp[i] = 2.0f * eqtl->af *(1.0 - eqtl->af);
            
            if(eqtl->af <= 0) LOGGER.e(0,"Allele frequency of eQTL SNP " + eqtl->rsID + " is negative ( af= " + std::to_string(eqtl->af) + ") ." );
            // flipped
            if(eqtl->flipped){
                gene->eQTLMarginEffect(i) = - gene->eQTLMarginEffect(i);
            } else {

            }
            /////////////////////////////////////////////
            // here we have various ways to get eqtl_n
            if(eqtl->eqtl_n != -999 ) {
                nTmpEst[i]  = eqtl->eqtl_n;
            } else if(gene->sampleSize != -999) {
                eqtl->eqtl_n = gene->sampleSize;
                nTmpEst[i]  = eqtl->eqtl_n;
            } else {
                if(!imputeSamBool){
                    imputeSamBool = true;
                    LOGGER << "Impute per snp sample size..." << endl;
                }
            }
            bTmp[i]  = gene->eQTLMarginEffect(i);
            seTmp[i] = gene->eQTLMarginEffectSE(i);
            // resign gene sample size
            gene->sampleSize = aimN;
            eqtl->eqtl_n = aimN;
        }

        MatrixXd localLD = eigenVecGene[j] * eigenValGene[j].asDiagonal() * eigenVecGene[j].transpose();
        localLD.diagonal().array() += (double)diag_mod;
        if(reSamType == 0) {
            reSampleEffectByImputation(betaImp,seImp,aimN,nTmpEst,snp2pqTmp,bTmp,seTmp,localLD);
        }
        else if(reSamType == 1) {
            reSampleEffectByBivariateLDSC(betaImp,seImp,aimN,nTmpEst,snp2pqTmp,bTmp,seTmp,localLD);
        }
        // cout  << gene->eQTLMarginEffect.size() << " " << gene->eQTLMarginEffect.head(10) << endl;
        gene->eQTLMarginEffect = betaImp;
        gene->eQTLMarginEffectSE = seImp;
    } // End of gene loop

}

void Data::reSampleEffectByImputation( VectorXd &betaImp, VectorXd &seImp,const int aimN, const VectorXd currN,const VectorXd snp2pq, const VectorXd beta, const VectorXd se, MatrixXd localLD){
    /// Step 2. Construct the LD correlation matrix among the typed SNPs(LDtt) and the LD correlation matrix among the missing SNPs and typed SNPs (LDit).
    // Step 2.1 divide SNPs into typed and untyped SNPs
    // VectorXi typedSnpIdx(gene->typSnpIdxGeneLD.size());
    // VectorXi untypedSnpIdx(gene->impSnpIdxGeneLD.size());
    // here we assume only one snp is missing and impute this based on other values
    Stat::Normal normal;
    unsigned numMarker = beta.size();
    vector<int> typedSnpIdx;
    vector<int> untypedSnpIdx;
    VectorXd zTypSnp(numMarker - 1);
    // VectorXd nTypSnp(numMarker - 1);
    VectorXd varyTypSnp(numMarker - 1);
    betaImp.setZero(numMarker);
    seImp.setZero(numMarker);
    for(unsigned i = 0; i < numMarker; i++){
        typedSnpIdx.clear();
        untypedSnpIdx.clear();
        // snpidx i is missing 
        for(unsigned k = 0; k < numMarker; k++){
            if(k == i) {
                untypedSnpIdx.push_back(k);
                continue;
            }
            typedSnpIdx.push_back(k);
        }
        zTypSnp.array() = beta(typedSnpIdx).array() / se(typedSnpIdx).array();
        // nTypSnp.array() = nTypSnp(typedSnpIdx);
        varyTypSnp.array() = snp2pq(typedSnpIdx).array() *
           (currN(typedSnpIdx).array() * se(typedSnpIdx).array().square() + 
           beta(typedSnpIdx).array().square());
        
        // Step 2.2 construct LDtt and LDit and Ztt.
        MatrixXd LDtt = localLD(typedSnpIdx,typedSnpIdx);
        MatrixXd LDit = localLD(untypedSnpIdx,typedSnpIdx);
        // Step 2.3 //  The Z score for the missing SNPs;
        VectorXd LDi_Z = LDtt.ldlt().solve(zTypSnp);
        // VectorXd zImpSnp = LDit * LDi_Z;
        VectorXd zImpSnp = LDit * LDi_Z;
        // Step 3. re-calcualte beta and se
        // if snp is missing use median to replace N
        // std::sort(nTypSnp.data(), nTypSnp.data() + nTypSnp.size());
        // double nMedian = nTypSnp[nTypSnp.size()/2];  // median
        double nMedian = aimN;
        // calcuate median of phenotypic variance
        std::sort(varyTypSnp.data(), varyTypSnp.data() + varyTypSnp.size());
        double varyMedian = varyTypSnp[varyTypSnp.size()/2];  // median
        // begin impute
        double base = sqrt(snp2pq(i) * (nMedian + zImpSnp(0) * zImpSnp(0)));
        betaImp(i) = zImpSnp(0) * sqrt(varyMedian)/base;
        seImp(i) = sqrt(varyMedian) / base;
    }
}

void Data::reSampleEffectByBivariateLDSC( VectorXd &betaImp, VectorXd &seImp,const int aimN, const VectorXd currN,const VectorXd snp2pq, const VectorXd beta, const VectorXd se, MatrixXd localLD){
    /// Step 2. Construct the LD correlation matrix among the typed SNPs(LDtt) and the LD correlation matrix among the missing SNPs and typed SNPs (LDit).
    // Step 2.1 divide SNPs into typed and untyped SNPs
    // VectorXi typedSnpIdx(gene->typSnpIdxGeneLD.size());
    // VectorXi untypedSnpIdx(gene->impSnpIdxGeneLD.size());
    // here we assume only one snp is missing and impute this based on other values
    Stat::Normal normal;
    unsigned numMarker = beta.size();
    vector<int> typedSnpIdx;
    vector<int> untypedSnpIdx;
    VectorXd zTypSnp(numMarker - 1);
    VectorXd nTypSnp(numMarker - 1);
    VectorXd varyTypSnp(numMarker);
    betaImp.setZero(numMarker);
    seImp.setZero(numMarker);
    
    // Cross trait LDSC reference: An atlas of genetic correlations across human diseases and traits
    // parameters
    double rhoPheCorr = 1.0;
    double rhoGeneticCorr = 1.0;
    for(unsigned i = 0; i < numMarker; i++){
        double N1 = currN(i);
        double N2 = aimN; // need to be resampled.
        double Ns;
        if(N2 <= N1){
            Ns = N2;
        } else {
            Ns = N1;
        }
        double ldscorei = localLD.col(i).array().sum();
        double sqrtN1N2 = sqrt(N1*N2);
        double z1 = beta(i)/se(i);
        double z1z2Mean  = sqrtN1N2 * rhoGeneticCorr * ldscorei / (double) numMarker + rhoPheCorr * Ns/ sqrtN1N2;
        double z2 = z1z2Mean/z2;
        // Step 3. re-calcualte beta and se
        // calculate vary 
        // begin impute
        betaImp(i) = z2/sqrt(snp2pq(i) *(N2 + z2 * z2 ));
        seImp(i) = 1/sqrt(snp2pq(i)*(N2 + z2 * z2));
    }
}

void Data::saveQueryCojoMaMoQTLInfo(const string title, const int slurmArrayLimit){
    // if(numKeptGenes > 1) {
    //     LOGGER.e(0, "There is no gene info in ma format, gene number should be one, but gene number is " + to_string(numKeptGenes) + ".");
    // }
    EqtlInfo * eqtl;
    GeneInfo * gene;
    map<string, EqtlInfo *>::iterator iterEqtl;
    int eqtlIdx;
    double sEqtlValue;

    std::stringstream dataStream;
    std::ostringstream genePhyPosRound,beta,se,chisq,pvalue,samSize;
    string header;
    double chiBuf, chisqBuf,pValueBuf,numNA=0;
    // int slurmArrayLimit = 1000;
    std::ofstream file;

    ///////////// output configure file ////////////////
        for(unsigned i = 0; i < numKeptGenes; i++){
        gene = keptGeneInfoVec[i];
        // use on decimal
        string esdfile = title + "-" + to_string(i) + "-" + gene->ensemblID + ".ma";
        int fileIdx;
        if(i % slurmArrayLimit == 0){
            fileIdx = i/slurmArrayLimit;
            string conFile = title + "-" + to_string(fileIdx) + ".txt" ;
            file.open(conFile);
            if (!file.is_open()) {LOGGER.e(0, " Can not open the file [" + esdfile + "] to save.");}
        }
        int endIdx = slurmArrayLimit * (i / slurmArrayLimit + 1);
        // cout << "i: " << i << " end: " << endIdx << endl;
        if(i < endIdx) {
            file << esdfile << endl;
        }
        if(i == endIdx - 1 || i == numKeptGenes -1){
            // cout << "endl;" << endl;
            LOGGER << " File " << fileIdx << "-" << endIdx -1 << " have been saved into file [ " << title + "-" + to_string(fileIdx) + ".txt"  << " ]." << endl;
            file.close();
        }
    }
    // string singleLine;
    // header line
    LOGGER << "Save eQTL data into plain query format..." << endl;
    header = "SNP\tA1\tA2\tfreq\tb\tse\tp\tN\n";
    bool haveSamSize = false;
    for(unsigned i = 0; i < numKeptGenes; i++){
        gene = keptGeneInfoVec[i];
        // use on decimal
        string esdfile = title + "-" + to_string(i) + "-" + gene->ensemblID + ".ma";
        std::ofstream file(esdfile);
        if (!file.is_open()) {LOGGER.e(0, " Can not open the file [" + esdfile + "] to save.");}
        file << header;
        genePhyPosRound.str("");
        genePhyPosRound << std::fixed << std::setprecision(1) << gene->midPhyPos;
        if(gene->cisSnpNameVec.size()!= gene->cisSnpIDSetInGene.size()){
            LOGGER.e(0,"Error, the size of cisSnpNameVec should be equal to size of cisSnpIDSetInGene for gene: [ " +gene->ensemblID + " ].");
        }
        // for(int j = 0; j < numIncdEqtls; j ++){
        //     eqtl = incdEqtlInfoVec[j];
        for(int j =0; j < gene->cisSnpNameVec.size(); j ++){
            iterEqtl = eqtlInfoMap.find(gene->cisSnpNameVec[j]);
            eqtl = iterEqtl->second;
            if(gene->cisSnpID2IdxMapInGene.find(eqtl->rsID) != gene->cisSnpID2IdxMapInGene.end()){
                eqtlIdx = gene->cisSnpID2IdxMapInGene[eqtl->rsID];
            } else {
                eqtl->included = false;
            }
            if(!eqtl->included) continue;
            // cout << "eqtlIdx: " << eqtlIdx << endl;
            beta.str("");
            se.str("");
            chisq.str("");
            pvalue.str("");
            samSize.str("");
            double snp2pqSqrt = sqrt(2.0f * eqtl->af *(1.0 - eqtl->af));
            // cout  << gene->eQTLMarginEffect.size() << " " << gene->eQTLMarginEffect.head(10) << endl;
            if(abs(gene->eQTLMarginEffect(eqtlIdx) ) < 1e-6 ){
                beta << std::scientific << std::setprecision(3) << gene->eQTLMarginEffect(eqtlIdx);
            } else{
                beta << std::fixed << std::setprecision(9) << gene->eQTLMarginEffect(eqtlIdx);
            }
            // se
            if(abs(gene->eQTLMarginEffectSE(eqtlIdx)) < 1e-6 ){
                se << std::scientific << std::setprecision(3) << gene->eQTLMarginEffectSE(eqtlIdx);
            } else{
                se << std::fixed << std::setprecision(9) << gene->eQTLMarginEffectSE(eqtlIdx);
            }
            if(abs(gene->eQTLMarginEffectSE(eqtlIdx)) > 1e-6){
                chiBuf = gene->eQTLMarginEffect(eqtlIdx)/gene->eQTLMarginEffectSE(eqtlIdx);
                chisqBuf = chiBuf * chiBuf;
                pValueBuf = Stat::ChiSq::pchisq(chisqBuf,1);
            }
            // chisq
            if(abs(chisqBuf) < 1e-6 ){
                chisq << std::scientific << std::setprecision(3) << chisqBuf;
            } else{
                chisq << std::fixed << std::setprecision(9) << chisqBuf;
            }
            // pvalue
            if(abs(pValueBuf) < 1e-6 ){
                pvalue << std::scientific << std::setprecision(3) << (pValueBuf);
            } else{
                pvalue << std::fixed << std::setprecision(9) << (pValueBuf);
            }
            // sample size 
            samSize.str("");
            bool hasSamSize = false;
            // use the same gene sample size for all snp within the gene
            samSize.str("");
            if(gene->sampleSize != -999){
                samSize.str("");
                samSize << std::fixed << std::setprecision(0) <<  (gene->sampleSize);
                hasSamSize = true;
            }
            // perSNP sample size will be changed if additional sample size provided.
            if(gene->cisSnpSampleSizeMap.find(eqtl->rsID) != gene->cisSnpSampleSizeMap.end()){
                samSize.str("");
                samSize << std::fixed << std::setprecision(0) <<  gene->cisSnpSampleSizeMap.at(eqtl->rsID);
                hasSamSize = true;
            } 
            // if non snp 
            if(!hasSamSize){
                samSize << "NA";
                numNA ++;
            }
            /// output the results
            if(eqtl->af + 999 < 1e-6){
                file << eqtl->rsID + "\t" + eqtl->a1 + "\t" + eqtl->a2 + "\t" + "NA" + "\t" +
                        beta.str() + "\t" + se.str() + "\t" + pvalue.str() + "\t" + samSize.str() + "\n";   
            } else {
                file << eqtl->rsID + "\t" + eqtl->a1 + "\t" + eqtl->a2 + "\t" + to_string(eqtl->af) + "\t" +
                        beta.str() + "\t" + se.str() + "\t" + pvalue.str() + "\t" + samSize.str() + "\n";   
            }
 
        }
        file.close();
    }
    if(numNA > 0)LOGGER.w(0, to_string(numNA) + " eQTLs have NA per SNP sample size." );
    // Now close the file stream
    LOGGER<<""<< numKeptGenes <<" genes and "<< numIncdEqtls <<" cis-SNPs have been saved in the annotation file." <<endl;  
}

void Data::filterMoQTLBasedOnCisWindow(const double cisRegionWind){
    // in this function, we will use variable cisSnpIDSetInGene
    // save important information and then select snps for each gene
    keptGeneInfoVec  = makeKeptGeneInfoVec(geneInfoVec);
    numKeptGenes = (unsigned) keptGeneInfoVec.size();
    incdEqtlInfoVec = makeIncdEqtlInfoVec(eqtlInfoVec);
    numIncdEqtls = (unsigned) incdEqtlInfoVec.size();
    int i,j;
    vector<locus_bp> cisSnpVec;
    EqtlInfo * eqtl;
    map<int, string> chrBeginSnp;
    map<int, string>  chrEndSnp;
    for (i = 1; i < numIncdEqtls; i++) {
        // eqtl = eqtlInfoVec[i];
        if(incdEqtlInfoVec[i]->chrom != incdEqtlInfoVec[i-1]->chrom){
            chrBeginSnp.insert(pair<int, string>(incdEqtlInfoVec[i]->chrom, incdEqtlInfoVec[i]->rsID));
            chrEndSnp.insert(pair<int, string>(incdEqtlInfoVec[i - 1]->chrom,incdEqtlInfoVec[i - 1]->rsID ));
        }
    }
    chrEndSnp.insert(pair<int, string>(incdEqtlInfoVec[numIncdEqtls - 1]->chrom,incdEqtlInfoVec[numIncdEqtls - 1]->rsID ));
    /////////////////////////////////////////
    // Step 2. Map snps to genes
    /////////////////////////////////////////
    LOGGER  << "Mapping the physical positions of genes to SNP data (gene boundaries: " << cisRegionWind / 1000 << "kb away from middle position of the gene) ..." << endl;
    vector<string> gene2snp_1(numKeptGenes), gene2snp_2(numKeptGenes);
    vector<locus_bp>::iterator iter;
    map<int, string>::iterator chrIter;
    GeneInfo *gene;
    for (i = 0; i < numIncdEqtls ; i++) {
        eqtl = incdEqtlInfoVec[i];
        // eqtl->included = false;
        eqtl->isInGene = false; //  we need to re-match cis-region 
        cisSnpVec.push_back(locus_bp(eqtl->rsID, eqtl->chrom, eqtl->physPos));
    }
    // cis-region as mid-position +/- 1mb
    #pragma omp parallel for private(iter, chrIter)
    for (i = 0; i < numKeptGenes; i++) {
        // find lowest snp_name in the gene
        gene = keptGeneInfoVec[i];
        gene->geneLength = cisRegionWind * 2;
        iter = find_if(cisSnpVec.begin(), cisSnpVec.end(), locus_bp( gene->ensemblID ,gene->chrom, gene->midPhyPos - cisRegionWind));
        if (iter != cisSnpVec.end()) gene2snp_1[i] = iter->locusName;
        else gene2snp_1[i] = "NA";
    }
    #pragma omp parallel for private(iter, chrIter)
    for (i = 0; i < numKeptGenes; i++) { 
        gene = keptGeneInfoVec[i];               
        if (gene2snp_1[i] == "NA") {
            gene2snp_2[i] = "NA";
            continue;
        }
        iter = find_if(cisSnpVec.begin(), cisSnpVec.end(), locus_bp(gene->ensemblID, gene->chrom, gene->midPhyPos + cisRegionWind));
        if (iter != cisSnpVec.end()){
            if (iter->bp ==  gene->midPhyPos + cisRegionWind){
                gene2snp_2[i] = iter->locusName;
            }else {
                if(iter!=cisSnpVec.begin()){
                    iter--;
                    gene2snp_2[i] = iter->locusName;
                }
                else gene2snp_2[i] = "NA";
            }
        }
        else {
            chrIter = chrEndSnp.find(gene->chrom);
            if (chrIter == chrEndSnp.end()) gene2snp_2[i] = "NA";
            else gene2snp_2[i] = chrIter->second;
        }
    }
    int mapped = 0;
    for (i = 0; i < numKeptGenes; i++) {
        gene = keptGeneInfoVec[i];
        if (gene2snp_1[i] != "NA" && gene2snp_2[i] != "NA") 
        {
            mapped++;
            gene->kept = true;
        }else {
            gene->kept = false;
        }
    }
    if (mapped < 1) LOGGER.e(0, "No gene can be mapped to the SNP data. Please check the input data regarding chromosome and bp.");
    
    map<string, int>::iterator iter1, iter2;
    map<string, int> cisSnpNameMap;
    vector<int> snpNumInGene(numKeptGenes);
    vector<VectorXd> cumsumNonNeg(numKeptGenes);
    map<int, ChromInfo*>::iterator iterChr;
    for (i = 0; i < numIncdEqtls; i++) {
        eqtl = incdEqtlInfoVec[i];
        cisSnpNameMap.insert(pair<string,int>(eqtl->rsID, i));
    }
    unsigned numTotCisSnp = 0;
    unsigned numGeneLeft = 0;
    for (i = 0; i < numKeptGenes; i++) {
        gene = keptGeneInfoVec[i]; 
        gene->cisSnpNameVec.clear();
        gene->cisSnpIDSetInGene.clear();
        gene->windows = cisRegionWind;
        gene->start = gene->midPhyPos - cisRegionWind;
        gene->end = gene->midPhyPos + cisRegionWind;
        iterChr = chromInfoMap.find(gene->chrom);
        if(iterChr != chromInfoMap.end()){
            if(gene->start < iterChr->second->startSnPhyPos) gene->start = iterChr->second->startSnPhyPos;
            if(gene->end > iterChr->second->endSnpPhyPos) gene->end = iterChr->second->endSnpPhyPos;
        }
        iter1 = cisSnpNameMap.find(gene2snp_1[i]);
        iter2 = cisSnpNameMap.find(gene2snp_2[i]);
        bool skip = false;
        if (iter1 == cisSnpNameMap.end() || iter2 == cisSnpNameMap.end() || iter1->second >= iter2->second) gene->kept = false;
        snpNumInGene[i] = iter2->second - iter1->second + 1;
        if(!gene->kept) continue;
        vector<int> cisSnpIdx;
        for (j = iter1->second; j <=  iter2->second; j++) {
            gene->cisSnpNameVec.push_back(incdEqtlInfoVec[j]->rsID);
            gene->cisSnpIDSetInGene.insert(incdEqtlInfoVec[j]->rsID);
            incdEqtlInfoVec[j]->isInGene = true;
        }
        numGeneLeft ++;
    } // end of loop genes

    // remove eqtls that do not belong to any gene
    for(unsigned i =0; i < numIncdEqtls; i++){
        eqtl = incdEqtlInfoVec[i];
        if(!eqtl->isInGene) {
            eqtl->included = false;
        }
    }
    // remove genes that do not have eQTLs
    for(unsigned i =0; i < numKeptGenes ; i++){
        gene = keptGeneInfoVec[i];
        if(!gene->hasEqtl) {
            gene->kept = false;
        }
    }
    // // remove non-cisEqtls;
    for (i = 0; i < numIncdEqtls; i++) {
        eqtl = incdEqtlInfoVec[i];
        if(eqtl->included){numTotCisSnp ++;}
    }
    LOGGER  << numGeneLeft << " genes and " << numTotCisSnp  <<   " cis-SNPs are left after using cis-window filter." << endl;
}

void Data::filterMoQTLBasedOnPvalue(const double pValueThreshold){
    keptGeneInfoVec  = makeKeptGeneInfoVec(geneInfoVec);
    numKeptGenes = (unsigned) keptGeneInfoVec.size();
    incdEqtlInfoVec = makeIncdEqtlInfoVec(eqtlInfoVec);
    numIncdEqtls = (unsigned) incdEqtlInfoVec.size();

    map<string, GeneInfo *>::iterator iterGene;
    map<string, EqtlInfo *>::iterator iterEqtl;
    EqtlInfo *eqtl;
    GeneInfo *gene;
    int eqtlIdx;
    vector<string> cisSnpNameVecLocal; // store information from based line as a bechmark
    set<string> cisSnpIDSetInGeneLocal;
    unsigned numTotCisSnp = 0;
    unsigned numGeneLeft = 0;
    double chiBuf, chisqBuf,pValueBuf;
    int numPvalPruned = 0;
    ///  gene 
    for (int i = 0; i < numIncdEqtls ; i++) {
        eqtl = incdEqtlInfoVec[i];
        eqtl->isInGene = false; //  we need to re-match cis-region 
    }
    for (int i = 0; i < numKeptGenes; i++) {
        gene = keptGeneInfoVec[i]; 
        if(!gene->kept) continue;
        cisSnpNameVecLocal.clear();
        cisSnpIDSetInGeneLocal.clear();
        vector<int> cisSnpIdx;
        // Step 1. loop all xqtl for each gene
        for(int j =0; j < gene->cisSnpNameVec.size(); j ++){
            iterEqtl = eqtlInfoMap.find(gene->cisSnpNameVec[j]);
            eqtl = iterEqtl->second;
            if(!eqtl->included) continue;
            // find xqtl in this gene;
            if(gene->cisSnpIDSetInGene.find(eqtl->rsID) != gene->cisSnpIDSetInGene.end()){
                if(gene->cisSnpID2IdxMapInGene.find(eqtl->rsID) != gene->cisSnpID2IdxMapInGene.end()){
                    eqtlIdx = gene->cisSnpID2IdxMapInGene[eqtl->rsID];
                } else {
                    eqtl->included = false;
                    LOGGER.e(0,"Could not find eqtl " +eqtl->rsID + " in gene " + gene->ensemblID + ".");
                }
                /// here we use pvalue threshold
                chiBuf = gene->eQTLMarginEffect(eqtlIdx)/gene->eQTLMarginEffectSE(eqtlIdx);
                chisqBuf = chiBuf * chiBuf;
                pValueBuf = Stat::ChiSq::pchisq(chisqBuf,1);
                if ( pValueBuf >= pValueThreshold ) {
                    continue;
                } 
                ///  collect significant snps into gene
                cisSnpNameVecLocal.push_back(eqtl->rsID);
                cisSnpIDSetInGeneLocal.insert(eqtl->rsID);

            }
        } // end of loop xqtl for one gene
        // Step 2 remove genes with less than 2 xQTLs
        if(cisSnpIDSetInGeneLocal.size() < 2) {
            gene->kept = false;
            gene->hasEqtl = false;
            gene->cisSnpNameVec.clear();
            gene->cisSnpIDSetInGene.clear();
        } else {
            for(int j =0; j < cisSnpNameVecLocal.size(); j ++){
                iterEqtl = eqtlInfoMap.find(gene->cisSnpNameVec[j]);
                iterEqtl->second->isInGene = true;
            }
            gene->cisSnpNameVec = cisSnpNameVecLocal;
            gene->cisSnpIDSetInGene = cisSnpIDSetInGeneLocal;
        }
    } // end of loop genes
    // remove eqtls that do not belong to any gene
    for(unsigned i =0; i < numIncdEqtls; i++){
        eqtl = incdEqtlInfoVec[i];
        if(!eqtl->isInGene) {
            eqtl->included = false;
            numPvalPruned ++;
            continue;
        }
        if(!eqtl->included) continue;
        numTotCisSnp ++;
    }
    // remove genes that do not have eQTLs
    for(unsigned i =0; i < numKeptGenes ; i++){
        gene = keptGeneInfoVec[i];
        if(!gene->hasEqtl) {
            gene->kept = false;
            continue;
        }
        if(!gene->kept) continue;
        numGeneLeft ++;
    }
    LOGGER << numPvalPruned << " SNPs are removed at the significant level (pvalue="<< pValueThreshold << ")." << endl;
    LOGGER  << numGeneLeft << " genes and " << numTotCisSnp  <<   " cis-SNPs are left at the significant level (pvalue="<< pValueThreshold << ")." << endl;
}

void Data::mergeMultiBesdData(const string besdListFile,const string title){
    ifstream in(besdListFile.c_str());
    if (!in) LOGGER.e(0, " Can not open the file [" + besdListFile + "] to read.");
    LOGGER << "Reading besd merge list file from [" + besdListFile + "]." << endl;

    vector<GeneInfo *> geneInfoVecLocal;
    GeneInfo * gene;
    EqtlInfo * eqtl;
    map<string, GeneInfo *> geneInfoMapLocal;
    vector<EqtlInfo *> eqtlInfoVecLocal;
    map<string, EqtlInfo *> eqtlInfoMapLocal;
    map<string, GeneInfo *>::iterator iterGene;
    map<string, EqtlInfo *>::iterator iterEqtl;
    bool epiBool, esiBool, besdBool;
    geneInfoVec.clear();
    geneInfoMap.clear();
    eqtlInfoVec.clear();
    eqtlInfoMap.clear();

    string oneFile;
    unsigned curSnpLengthInGene = 0;
    unsigned snpLenInGeneLocal;
    bool smrBesdBool = false;
    bool hasSnpInfo = false;
    bool hasGeneLdmInfo = false;
    // baseline for first file
    // getline(in, oneFile);
    // epiBool = readMultiEpiFile(oneFile + ".epi", geneInfoVec,geneInfoMap);
    // esiBool = readMultiEsiFile(oneFile + ".esi",eqtlInfoVec, eqtlInfoMap,smrBesdBool);
    // if (epiBool && esiBool)
    //     besdBool = readMultiBesdFile(  + ".besd",geneInfoVec,geneInfoMap,eqtlInfoVec, eqtlInfoMap,smrBesdBool);

    // other files
    while(getline(in, oneFile)){
        if (oneFile.empty()) {continue;} // Skip the empty line
        geneInfoVecLocal.clear();
        geneInfoMapLocal.clear();
        eqtlInfoVecLocal.clear();
        eqtlInfoMapLocal.clear();
        epiBool = readMultiEpiFile(oneFile + ".epi", geneInfoVecLocal,geneInfoMapLocal,hasSnpInfo,hasGeneLdmInfo);
        esiBool = readMultiEsiFile(oneFile + ".esi",eqtlInfoVecLocal, eqtlInfoMapLocal,smrBesdBool,hasSnpInfo,hasGeneLdmInfo);
        if(epiBool && esiBool)
            besdBool = readMultiBesdFile(oneFile + ".besd",geneInfoVecLocal,geneInfoMapLocal,eqtlInfoVecLocal, eqtlInfoMapLocal,smrBesdBool,hasSnpInfo,hasGeneLdmInfo);
        if(besdBool){
            curSnpLengthInGene = 0;
            for(unsigned j = 0; j < geneInfoVecLocal.size();j++){
                gene = geneInfoVecLocal[j];
                iterGene = geneInfoMap.find(gene->ensemblID);
                if(iterGene != geneInfoMap.end()){
                    // already exists
                    // snpName;
                    iterGene->second->cisSnpNameVec.insert(iterGene->second->cisSnpNameVec.begin(),gene->cisSnpNameVec.begin(),gene->cisSnpNameVec.end());
                    // 1. beta 
                    iterGene->second->eQTLMarginEffect.conservativeResize(iterGene->second->numSnpInGene + gene->numSnpInGene);
                    iterGene->second->eQTLMarginEffect.segment(iterGene->second->numSnpInGene,gene->numSnpInGene) = gene->eQTLMarginEffect;
                    // 2. se 
                    iterGene->second->eQTLMarginEffectSE.conservativeResize(iterGene->second->numSnpInGene + gene->numSnpInGene);
                    iterGene->second->eQTLMarginEffectSE.segment(iterGene->second->numSnpInGene,gene->numSnpInGene) = gene->eQTLMarginEffectSE;
                    // 3. snpID to Idx map
                    for(unsigned k = 0; k < gene->numSnpInGene; k++){
                        iterGene->second->cisSnpID2IdxMapInGene.insert(pair<string,int>(gene->cisSnpNameVec[k],iterGene->second->numSnpInGene + k));
                    }
                    // 4. insert current set
                    iterGene->second->cisSnpIDSetInGene.insert(gene->cisSnpIDSetInGene.begin(), gene->cisSnpIDSetInGene.end());
                    // 5. snps number in gene
                    iterGene->second->numSnpInGene += gene->numSnpInGene;
                    // iterGene->second->eQTLMarginEffect.insert();
                } else {
                    geneInfoVec.push_back(gene);
                    geneInfoMap.insert(pair<string, GeneInfo *>(gene->ensemblID,gene));
                }

            }
            for(unsigned j = 0; j < eqtlInfoVecLocal.size();j++){
                eqtl = eqtlInfoVecLocal[j];
                iterEqtl = eqtlInfoMap.find(eqtl->rsID);
                if(iterEqtl != eqtlInfoMap.end()){
                    // already exists
                    // LOGGER.e(0,"Error: duplicated ID " + eqtl->rsID  + " is found.");
                } else {
                    eqtlInfoVec.push_back(eqtl);
                    eqtlInfoMap.insert(pair<string, EqtlInfo *>(eqtl->rsID,eqtl));
                }
            }
        }
    }  // end of while 
    numeQTLs = eqtlInfoVec.size();
    incdEqtlInfoVec = makeIncdEqtlInfoVec(eqtlInfoVec);
    numIncdEqtls = (unsigned) incdEqtlInfoVec.size();

    numGenes = geneInfoVec.size();
    keptGeneInfoVec  = makeKeptGeneInfoVec(geneInfoVec);
    numKeptGenes = (unsigned) keptGeneInfoVec.size();
    LOGGER<<"Eventually, eQTL summary data of "<< numKeptGenes <<" Probes and "<< numIncdEqtls <<" SNPs to be included from [" + besdListFile + "]." << endl;

    saveBesdFormat(title,smrBesdBool); // 
}
bool Data::readMultiEpiFile(const string &epiFile, vector<GeneInfo *> &geneInfoVecLocal,map<string, GeneInfo *> &geneInfoMapLocal,const bool hasSnpInfo,const bool hasGeneLdmInfo){
    ifstream in(epiFile.c_str());
    if (!in){
        LOGGER.w(0, " can not open the file [" + epiFile + "] to read.");
        return false;
    } 
    LOGGER << "Reading BESD format epi file from [" + epiFile + "]." << endl;

    map<string, GeneInfo*>::iterator iterGene;

    unsigned chr, sampleSize;
    double geneticDis;
    unsigned physPos;
    string probe, geneID, geneOri;
    unsigned idx = 0;
    Gadget::Tokenizer colData;string inputStr;string sep(" \t\n");
    while (getline(in,inputStr)) {
        colData.getTokens(inputStr, sep);
        if(idx == 0){
            if(colData.size() == 7){
            } else {
                LOGGER.e(0,"Wrong epi format.");
            }
        }
        chr = std::stoi(colData[0]);
        geneID = colData[1];
        GeneInfo *gene = new GeneInfo(idx++, geneID, chr);
        gene->probeID = colData[4];
        gene->geneticDis = std::stod(colData[2]);
        gene->start  = std::stod(colData[2]);
        gene->midPhyPos = std::stod(colData[3]);
        gene->sampleSize = std::stoi(colData[5]);
        gene->geneOri = colData[6];
        geneInfoVecLocal.push_back(gene);
        if (geneInfoMapLocal.insert(pair<string, GeneInfo*>(geneID, gene)).second == false) {
            LOGGER.w(0, " Duplicate gene ID found: \"" + geneID + "\".");
            return false;
        }
        if(hasGeneLdmInfo){
            // here we need to keep gene that have xQTL LD info only
            iterGene = geneInfoMap.find(gene->ensemblID); // xQTL LD gene info
            if(iterGene != geneInfoMap.end()){
                // TODO other gene consistent check;
                iterGene->second->sampleSize = gene->sampleSize;
            } else {
                gene->kept = false;
            }
        }
    }
    in.close();
    unsigned numGenesLocal = (unsigned) geneInfoVecLocal.size();
    LOGGER << numGenesLocal << " genes to be included from [" + epiFile + "]." << endl;
    return true;

}
bool Data::readMultiEsiFile(const string &esiFile,vector<EqtlInfo *> &eqtlInfoVecLocal,map<string, EqtlInfo *> &eqtlInfoMapLocal,bool &smrBesdBool,const bool hasSnpInfo,const bool hasGeneLdmInfo){
    ifstream in(esiFile.c_str());
    if (!in) {
        LOGGER.w(0, " can not open the file [" + esiFile + "] to read.");
        return false;
    }
    LOGGER << "Reading BESD format esi file from [" + esiFile + "]." << endl;
    string id, allele1, allele2;
    unsigned chr, physPos;
    double genPos;
    double freq;
    double af;
    int sampleSize;
    unsigned idx = 0;
    // QC variable start
    SnpInfo * snp;
    EqtlInfo *eqtlLD;
    map<string, SnpInfo*>::iterator iterSnp;
    unsigned line=0, match=0;
    unsigned numInconAllele=0, numInconAf=0, numFixed=0, numMafMin=0, numMafMax=0, numOutlierN=0;
    unsigned numPvalPruned=0;
    unsigned numFlip=0;
    bool inconAllele, inconAf, fixed, ismafmin, ismafmax, isPvalPruned;
    // QC variable end
    Gadget::Tokenizer colData;string inputStr;string sep(" \t\n");
    while (getline(in,inputStr)) {
        colData.getTokens(inputStr, sep);
        if(idx == 0){
            if(colData.size() == 7 ){
                smrBesdBool = true;
                LOGGER << "SMR BESD format is used and per snp sample size is missing." << endl;
            } else if (colData.size() == 8) {
                smrBesdBool = false;
                LOGGER << "A modified SMR BESD format is used and per snp sample size can be read in esi file." << endl;
            } else {
                LOGGER.e(0,"Wrong esi format.");
            }
        }
        chr = std::stoi(colData[0]);
        id = colData[1];
        genPos = std::stod(colData[2]);
        physPos = ::stoi(colData[3]);
        allele1 = colData[4];
        allele2 = colData[5];
        freq = std::stod(colData[6]);
        EqtlInfo *eqtl = new EqtlInfo(idx++, 
        id, 
        allele1, 
        allele2,
        chr, 
        genPos, 
        physPos,
        freq
        );
        if(!smrBesdBool)eqtl->eqtl_n = std::stoi(colData[7]);
        eqtlInfoVecLocal.push_back(eqtl);
        if (eqtlInfoMapLocal.insert(pair<string, EqtlInfo*>(id, eqtl)).second == false) {
            LOGGER.w(0, " Duplicate eQTL SNP ID found: \"" + id + "\".");
            return false;
        }
        if(hasSnpInfo){
            // we use this to remove eqtls even if snp exist. e.g. we use --keep-block flag to select snps from 1 LD block.
            iterSnp = snpInfoMap.find(eqtl->rsID);
            if(iterSnp != snpInfoMap.end()){
                // Found in the snplist
                eqtl->included = iterSnp->second->included;
            }else {
                // Not found in the snplist
                eqtl->included = false;
            }
        }

        if(hasGeneLdmInfo){
            // here we need to QC for eQTLs, but here we only check allele now.
            // we use SNP LD as benchmard
            iterSnp = snpInfoMap.find(eqtl->rsID);
            if(iterSnp != snpInfoMap.end()){
                snp = iterSnp->second;
                if (allele1 == snp->a1 && allele2 == snp->a2) {
                } else if (allele1 == snp->a2 && allele2 == snp->a1) {
                    eqtl->flipped = true;
                    ++numFlip;
                } else {
                    inconAllele = true;
                    ++numInconAllele;
                    eqtl->included = false;
                }
                // now we need to use eqtl_n and af from BESD to replace that from eqtl
            } // end of iterSnp.
        } // end of hasGeneLdmInfo
        
    }
    in.close();
    unsigned numeQTLsLocal = (unsigned) eqtlInfoVecLocal.size();
    LOGGER << numeQTLsLocal << " eQTLs to be included from [" + esiFile + "]." << endl;

    // summary
    if(hasGeneLdmInfo) {
        numEqtlFlip = numFlip;
        if (numFlip) LOGGER << "flipped " << numFlip << " SNPs according to the minor allele in the reference and GWAS samples." << endl;
    }

    return true;
}
bool Data::readMultiBesdFile(const string &besdFile,vector<GeneInfo *> &geneInfoVecLocal,map<string, GeneInfo *> &geneInfoMapLocal,
                    vector<EqtlInfo *> &eqtlInfoVecLocal,map<string, EqtlInfo *> &eqtlInfoMapLocal, bool &smrBesdBool,const bool hasSnpInfo,const bool hasGeneLdmInfo){
    unsigned numGenesLocal = geneInfoVecLocal.size();
    unsigned numeQTLsLocal = eqtlInfoVecLocal.size();
    if(numeQTLsLocal == 0) {LOGGER.w(0," No SNP is retained for eQTL analysis."); return false;}
    if(numGenesLocal == 0) {LOGGER.w(0," No gene is retained for eQTL analysis. ");return false;}
    ifstream besd(besdFile.c_str(), ios::in|ios::binary);
    if (!besd) {LOGGER.w(0, " can not open the file [" + besdFile + "] to read.");return false;}
    LOGGER << "Reading eQTL summary data from BESD format BED file [" + besdFile + "] ..." << endl;
    
    char SIGN[sizeof(uint64_t)+8];
    besd.read(SIGN,4);
    uint32_t header = *(uint32_t *)SIGN;
    if(header == 0x40000000 & header == 0x3f800000){
         LOGGER.w(0, "This is an old BESD format. Please use smr --make-besd command to update the file format." );
        //  exit(1);
        return false;
    }
    // Read eQTL summary 
	map<string, SnpInfo *>::iterator iterSnp;
    EqtlInfo *eqtl = NULL;
    GeneInfo *gene = NULL;
    unsigned geneIdx,eqtlIdx;
    float eqtlBeta;
    float eqtlSE;
    int sampleSize;
    unsigned long long skip = 0;

    /////////////////////////////////////////////////////////////
    float totalSteps = 0; // Total number of steps in your process
    /////////////////////////////////////////////////////////////
    // read dense format
    if (header  == SPARSE_FILE_TYPE_3F || header == SPARSE_FILE_TYPE_3){
        LOGGER << "Read eqtl info from smr sparse besd format." << endl;
        char* buffer;
        uint64_t colNum=(numGenesLocal<<1)+1;
        uint64_t valNum;
        uint64_t lSize;
        besd.seekg(0,besd.end);
        lSize = besd.tellg();
        besd.seekg(4); //  same as besd.seekg(4, besd.beg);
        
        if(header==SPARSE_FILE_TYPE_3){
            int length=(RESERVEDUNITS-1)*sizeof(int);
            char* indicators=new char[length];
            besd.read(indicators,length);
            int* tmp=(int *)indicators;
            sampleSize=*tmp++;
            if(sampleSize!=-9){
                LOGGER << "The eQTL sample size is " << sampleSize << endl;
                // numKeptIndseQTL = sampleSize;
            }else {
                // LOGGER.e(0,"ERROR:The sample size is missing. You may use (smr --beqtl-summary myeqtl --add-n 1000 --make-besd --out mybesd) to add it to the BESD file.");
                LOGGER << "Warning:The sample size is missing. You may use (smr --beqtl-summary myeqtl --add-n 1000 --make-besd --out mybesd) to add it to the BESD file." << endl;
            }
            if(*tmp++!=numeQTLsLocal){
                LOGGER.w(0,"ERROR: The SNPs in your .esi file are not in consistency with the one in .besd file.");
                return false;
            }
            if(*tmp++!=numGenesLocal)
            {
                LOGGER.w(0,"ERROR: The probes in your .epi file are not in consistency with the one in .besd file");
                return false;
            }
            delete[] indicators;
        }
        besd.read(SIGN, sizeof(uint64_t));
        valNum=*(uint64_t *)SIGN;

        if(header==SPARSE_FILE_TYPE_3F) {
            if( lSize - (sizeof(uint32_t) + sizeof(uint64_t) + colNum*sizeof(uint64_t) + valNum*sizeof(uint32_t) + valNum*sizeof(float)) != 0){
                // cout << " total size: " << sizeof(uint32_t) + sizeof(uint64_t) + colNum*sizeof(uint64_t) + valNum*sizeof(uint32_t) + valNum*sizeof(float) << endl;
                LOGGER << "The file size is " <<  lSize;
                LOGGER << " " << to_string(sizeof(uint32_t)) << " " <<  to_string(sizeof(uint64_t)) << " " << to_string(colNum*sizeof(uint64_t)) << " "
                     << to_string(valNum*sizeof(uint32_t)) << " " << to_string(valNum*sizeof(float)) << endl;
                LOGGER.w(0,"ERROR: wrong value number. File" + besdFile + " is ruined.\n");
                return false;
            }
        }else{
            if( lSize - (RESERVEDUNITS*sizeof(int) + sizeof(uint64_t) + colNum*sizeof(uint64_t) + valNum*sizeof(uint32_t) + valNum*sizeof(float)) != 0){
                LOGGER << "The file size is " <<  lSize;
                LOGGER << " " << to_string(RESERVEDUNITS*sizeof(int)) << " " << to_string(sizeof(uint64_t)) << " " << to_string(colNum*sizeof(uint64_t)) << " "
                     << to_string(valNum*sizeof(uint32_t)) << " " << to_string(valNum*sizeof(float)) << endl;
                LOGGER.w(0,"ERROR: wrong value number. File" + besdFile + "is ruined.\n");
                return false;
            }    
        }
        //for sparse
        vector<uint64_t> _cols;
        vector<uint32_t> _rowid;
        //vector <double>  _val;
        VectorXd _val;
        //////////////////////
        // Step 1, read all eqtl data
        buffer = (char*) malloc  (sizeof(char)*(lSize));
        if (buffer == NULL) {fputs ("Memory error",stderr); exit (1);}
        besd.read(buffer,lSize);
        if(header ==SPARSE_FILE_TYPE_3F)
        {
        if (besd.gcount()+sizeof(uint32_t) + sizeof(uint64_t) != lSize){LOGGER.w(0," cannot read besd file: " + besdFile );return false;}
        }else {
            if (besd.gcount()+RESERVEDUNITS*sizeof(int) + sizeof(uint64_t) != lSize){ LOGGER.w(0," cannot read besd file: " + besdFile );return false;}
        }
        
        _cols.resize(colNum);
        _rowid.resize(valNum);
        _val.resize(valNum); 
        uint64_t* ptr;
        ptr=(uint64_t *)buffer;
        for(int i=0;i<colNum;i++) {_cols[i]=*ptr++; } 
        uint32_t* ptr4B=(uint32_t *)ptr;
        for(int i=0;i<valNum;i++){_rowid[i]=*ptr4B++;}
        float* val_ptr=(float*)ptr4B;
        for(int i=0;i<valNum;i++){_val(i) =*val_ptr++;}

        // Step 2. build a map to save eqtls passed QC

        map<int,string> eqtlIndxMap;

        map<int, string>::iterator Intiter;
        map<string, EqtlInfo* >::iterator iterEqtl;
		map<string, SnpInfo *>::iterator iterSnp;
        // for(int i = 0; i < numIncdEqtls;i++){
        //     eqtl = incdEqtlInfoVec[i];
        for(int i = 0; i < numeQTLsLocal;i++){
            eqtl = eqtlInfoVecLocal[i];
            // if(smrBesdBool) eqtl->eqtl_n = sampleSize;
            eqtlIndxMap.insert(pair<int,string> (i,eqtl->rsID));
        }
        //////////////////////////////////////////////
        totalSteps = _cols.size();
        //////////////////////////////////////////////
        // Step 2. select subset of eqtl info based one kept eqtl and gene list
        int size = 2; 
        for(unsigned j = 0, geneIdx = 0; j < _cols.size() - 1; j += size ){
            gene = geneInfoVecLocal[j/2];
            // if(!hasGeneLdmInfo){
                // these are global parameter, when gene ld is used, they will be set by info from LD
                gene->cisSnpNameVec.clear();
                gene->gene2CisSnpVec.clear();
                gene->cisSnpIDSetInGene.clear();
                gene->cisSnpID2IdxMapInGene.clear();
            // }
            gene->numSnpInGene = 0;
            if(! gene->kept){
                skip += size;
                continue;
            }
            // if(skip) j = j + skip;
            skip = 0;
            int numEqtlInGene = _cols[j + 1] - _cols[j];

            VectorXf eqtlEffect(numEqtlInGene), eqtlEffectSE(numEqtlInGene);
            vector<int> eqtlEffectIndxInGene;
            eqtlEffectIndxInGene.clear();
            // remove eqtls not pass QC
            int localIdx,icldIdx;
            eqtlEffect.setZero();
            eqtlEffectSE.setZero();
            for(unsigned i  = _cols[j], localIdx  = 0,icldIdx = 0; i < _cols[j + 1] ; i ++){
                if(eqtlIndxMap.find(_rowid[i]) != eqtlIndxMap.end()){
                    iterEqtl = eqtlInfoMapLocal.find(eqtlIndxMap.at(_rowid[i]));
                    eqtl = iterEqtl->second;
                    if(!eqtl->included){localIdx ++; continue;}
                    eqtlEffect(localIdx) = _val[i];
                    eqtlEffectSE(localIdx) = _val[i + numEqtlInGene];
                    if(fabs(eqtlEffectSE(localIdx) + 9) < 1e-6 ) {localIdx ++; continue;}

                    eqtlEffectIndxInGene.push_back(localIdx);
                    // if(!hasGeneLdmInfo){
                        // global parameter
                        gene->cisSnpNameVec.push_back(eqtl->rsID);
                        gene->gene2CisSnpVec.push_back(_rowid[i]);
                        gene->cisSnpIDSetInGene.insert(eqtl->rsID);
                        gene->cisSnpID2IdxMapInGene.insert(pair<string,int>(eqtl->rsID,icldIdx));
                    // } else {
                    //     // besed related parameter
                    //     gene->cisSnpNameVecBESD.push_back(eqtl->rsID);
                    //     gene->cisSnpID2IdxMapInGeneBESD.insert(pair<string,int>(eqtl->rsID,icldIdx));
                    // }
                    gene->hasEqtl = true;
                    eqtl->isInGene = true;
                    icldIdx ++;
                } //end of if(eqtlIndxMap.find(_rowid[i]) != eqtlIndxMap.end()){
                localIdx ++;
            } // end of loop beta and se
            gene->numSnpInGene = eqtlEffectIndxInGene.size();
            gene->eQTLMarginEffect = eqtlEffect(eqtlEffectIndxInGene).cast<double>();
            gene->eQTLMarginEffectSE = eqtlEffectSE(eqtlEffectIndxInGene).cast<double>();

            ///////////////////////////////////////////////////////////////////////
            float currentStep = j;
            float progress = (currentStep * 100.0) / totalSteps;
            LOGGER.p(0, " " + to_string(static_cast<float>(progress)) + "% ","Read eqtl info from sparse besd format:" );
            ///////////////////////////////////////////////////////////////////////
            // remove genes if there is no eqtl in gene 
            if (gene->gene2CisSnpVec.size() == 0){ gene->kept = false;}
            if (++geneIdx == numGenesLocal) break;
        }        
        // terminate     
        // terminate
        _cols.clear();
        _rowid.clear();
        _val.resize(0); 
        free (buffer);

    }else {
        LOGGER.w(0,"GCTB doesn't support this format. Please use OSCA (http://cnsgenomics.com/software/osca) to transform it to SMR format.");
        return false;
    }
    
    // remove eqtls that do not belong to any gene
    for(unsigned i =0; i < numeQTLsLocal; i++){
        eqtl = eqtlInfoVecLocal[i];
        if(!eqtl->isInGene) {
            eqtl->included = false;
        }
    }
    // remove genes that do not have eQTLs
    for(unsigned i =0; i < numGenesLocal ; i++){
        gene = geneInfoVecLocal[i];
        if(!gene->hasEqtl) {
            gene->kept = false;
        }
    }

    geneInfoVecLocal  = makeKeptGeneInfoVec(geneInfoVecLocal);
    numGenesLocal = (unsigned) geneInfoVecLocal.size();
    eqtlInfoVecLocal = makeIncdEqtlInfoVec(eqtlInfoVecLocal);
    numeQTLsLocal = (unsigned) eqtlInfoVecLocal.size();

    LOGGER<<"eQTL summary data of "<< numGenesLocal <<" Probes and "<< numeQTLsLocal <<" SNPs to be included from [" + besdFile + "]." << endl;
    return true;
}

void Data::alignxQTLGeneAndBESDGeneInfoWithoutImp(vector<GeneInfo *> &geneInfoVecBESD,map<string, GeneInfo *> &geneInfoMapBESD,
        vector<EqtlInfo *> &eqtlInfoVecBESD,map<string, EqtlInfo *> &eqtlInfoMapBESD, double diag_mod){
        // align xQTL effect and se.
        GeneInfo * gene, *geneBESD; // gene is info from xQTL LD
        EqtlInfo * eqtl, *eqtlBESD;
        map<string, GeneInfo*>::iterator iterGene;
        map<string, EqtlInfo*>::iterator iterEqtl;
        string eqtlID;
        vector<int> gene2CisSnpVecBESDLocal;
        for(unsigned j = 0; j < numKeptGenes; j++){
            gene = keptGeneInfoVec[j];
            iterGene = geneInfoMapBESD.find(gene->ensemblID);
            if(iterGene == geneInfoMapBESD.end()){
                gene->kept = false; // cannot find LD gene in BESD gene;
                continue;
            }
            geneBESD = iterGene->second;

            // now we need to check xQTL missing situation and reorder xQTL b and se based on gene LD 
            gene->typSnpIdxBesd.clear(); // used to extract eqtl beta and se values;
            gene->typSnpIdxGeneLD.clear(); // used to extract eqtl with beta idx in gene ld 
            gene->impSnpIdxGeneLD.clear(); // used to extract eqtl idx need to impute
            gene2CisSnpVecBESDLocal.clear();
            for(unsigned i = 0; i < gene->cisSnpNameVec.size();i ++)
            {
                eqtlID = gene->cisSnpNameVec[i];
                if(geneBESD->cisSnpID2IdxMapInGene.find(eqtlID) == geneBESD->cisSnpID2IdxMapInGene.end()){
                    // missing;
                    gene->impSnpIdxGeneLD.push_back(i);
                    gene->impCisSnpEffBool = true;
                } else {
                    // found
                    // extract se and b index from besd data 
                    int eqtlIdxBESD = geneBESD->cisSnpID2IdxMapInGene.find(eqtlID)->second;
                    gene2CisSnpVecBESDLocal.push_back(eqtlIdxBESD);
                    gene->typSnpIdxGeneLD.push_back(i);
                    // extract sample size and allele frequency to besd;
                    eqtl = eqtlInfoMap.find(eqtlID)->second;
                    eqtlBESD = eqtlInfoMapBESD.find(eqtlID)->second;
                    if(eqtlBESD->eqtl_n + 999 < 1e-6 ) eqtl->eqtl_n = eqtlBESD->eqtl_n;
                    if(eqtlBESD->af > 0 ) eqtl->af = eqtlBESD->af;
                }
            } // end of loop cisSnpNameVec
            ///

            if(gene->impCisSnpEffBool){
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
            } else {
                gene->eQTLMarginEffect = geneBESD->eQTLMarginEffect(gene2CisSnpVecBESDLocal);
                gene->eQTLMarginEffectSE = geneBESD->eQTLMarginEffectSE(gene2CisSnpVecBESDLocal);
            } 
        } // end of gene->impCisSnpEffBool

    }

void Data::mergeMultiEigenMat(const string title, const string besdListFile,const string LDMatType, const double &eigenCutoff){
    ifstream in(besdListFile.c_str());
    if (!in) LOGGER.e(0, " Can not open the file [" + besdListFile + "] to read.");
    LOGGER << "Reading eigen list file from [" + besdListFile + "]." << endl;

    vector<GeneInfo *> geneInfoVecLocal;
    GeneInfo * gene;
    EqtlInfo * eqtl;
    map<string, GeneInfo *> geneInfoMapLocal;
    vector<EqtlInfo *> eqtlInfoVecLocal;
    map<string, EqtlInfo *> eqtlInfoMapLocal;
    map<string, GeneInfo *>::iterator iterGene;
    map<string, EqtlInfo *>::iterator iterEqtl;
    bool geneBool, snpBool, binBool;
    geneInfoVec.clear();
    geneInfoMap.clear();
    eqtlInfoVec.clear();
    eqtlInfoMap.clear();

    string oneFile;
    unsigned curSnpLengthInGene = 0;
    unsigned snpLenInGeneLocal;
    // baseline for first file
    getline(in, oneFile);
    bool smrBesdBool = false;
    geneBool = readMultiEigenMatInfoFile(oneFile + ".eigen.gene.info", geneInfoVec,geneInfoMap);
    snpBool = readMultiEigenMatSnpInfoFile(oneFile + ".eigen.gene.snp.info",eqtlInfoVec, eqtlInfoMap);
    if (geneBool && snpBool)
        binBool = readMultiEigenMatBinFile(oneFile + ".eigen.gene.bin",geneInfoVec,geneInfoMap,eqtlInfoVec, eqtlInfoMap,eigenCutoff);

    while(getline(in, oneFile)){
        if (oneFile.empty()) {continue;} // Skip the empty line
        geneInfoVecLocal.clear();
        geneInfoMapLocal.clear();
        eqtlInfoVecLocal.clear();
        eqtlInfoMapLocal.clear();
        geneBool = readMultiEigenMatInfoFile(oneFile + ".eigen.gene.info", geneInfoVecLocal,geneInfoMapLocal);
        snpBool = readMultiEigenMatSnpInfoFile(oneFile + ".eigen.gene.snp.info",eqtlInfoVecLocal, eqtlInfoMapLocal);
        if(geneBool && snpBool){
            binBool = readMultiEigenMatBinFile(oneFile + ".eigen.gene.bin",geneInfoVecLocal,geneInfoMapLocal,eqtlInfoVecLocal, eqtlInfoMapLocal,eigenCutoff);
        } else {
            binBool = false;
        }
        if(binBool){
            curSnpLengthInGene = 0;
            for(unsigned j = 0; j < geneInfoVecLocal.size();j++){
                gene = geneInfoVecLocal[j];
                iterGene = geneInfoMap.find(gene->ensemblID);
                if(iterGene != geneInfoMap.end()){
                    // already exists
                    LOGGER.w(0, "Duplicate gene [ " + iterGene->second->ensemblID + " ] is found. Remove it");
                    iterGene->second->kept = false;
                    // iterGene->second->eQTLMarginEffect.insert();
                } else {
                    geneInfoVec.push_back(gene);
                    geneInfoMap.insert(pair<string, GeneInfo *>(gene->ensemblID,gene));
                }

            }
            for(unsigned j = 0; j < eqtlInfoVecLocal.size();j++){
                eqtl = eqtlInfoVecLocal[j];
                iterEqtl = eqtlInfoMap.find(eqtl->rsID);
                if(iterEqtl != eqtlInfoMap.end()){
                    // already exists
                    // LOGGER.e(0,"Error: duplicated ID " + eqtl->rsID  + " is found.");
                } else {
                    eqtlInfoVec.push_back(eqtl);
                    eqtlInfoMap.insert(pair<string, EqtlInfo *>(eqtl->rsID,eqtl));
                }
            }
        }
    }  // end of while 
    numeQTLs = eqtlInfoVec.size();
    incdEqtlInfoVec = makeIncdEqtlInfoVec(eqtlInfoVec);
    numIncdEqtls = (unsigned) incdEqtlInfoVec.size();

    numGenes = geneInfoVec.size();
    keptGeneInfoVec  = makeKeptGeneInfoVec(geneInfoVec);
    numKeptGenes = (unsigned) keptGeneInfoVec.size();
    LOGGER<< "Eventually, eQTL summary data of " << numKeptGenes <<" Probes and "<< numIncdEqtls <<" SNPs to be included from [" + besdListFile + "]." << endl;

    // save multiple eigen gene 
    LOGGER << "Save multiple gene low-rank LD matrice..." << endl;
    MergeMultiBinEigen(title,LDMatType);
    outputEigenMatIndInfoFromLDMat(LDMatType, title,true); 
}

bool Data::readMultiEigenMatInfoFile(const string title,vector<GeneInfo *> &geneInfoVecLD, map<string, GeneInfo *> &geneInfoMapLD){
        // Read bim file: recombination rate is defined between SNP i and SNP i-1
    ifstream in(title.c_str());
    if (!in){
        LOGGER.w(0, " can not open the file [" + title + "] to read.");
        return false;
    }
    LOGGER << "Reading eigen gene info from file [" + title + "]." << endl;
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
            LOGGER.w(0, " Duplicate LDBlock-SNP pair found: \"" + id + "_" + snp + "\".");
            return false;
        } else{
            // Step 1.1 build a new geneInfo firstly
            if(snpCount == 1){
                gene = new GeneInfo(idx++, id, chr);
                gene->start = start;
                gene->end   = end;
                gene->typSnpIdxBesd.clear(); // used to extract eqtl beta and se values;
                gene->typSnpIdxGeneLD.clear(); // used to extract eqtl with beta idx in gene ld 
                gene->impSnpIdxGeneLD.clear(); // used to extract eqtl idx need to impute
                gene->cisSnpNameVec.clear();
                gene->cisSnpID2IdxMapInGene.clear();
            }
            // Step 1.2 continue to add eQTLs into gene
            if(snpCount < snpNum){
                gene->cisSnpNameVec.push_back(snp);
                gene->cisSnpID2IdxMapInGene.insert(pair<string,int>(snp,snpCount -1));
            }
            // Step 1.3 add the last eQTL into gene 
            if(snpCount == snpNum){
                gene->cisSnpNameVec.push_back(snp);
                gene->cisSnpID2IdxMapInGene.insert(pair<string,int>(snp,snpCount -1));
                gene->numSnpInGene = snpNum;
                geneInfoVecLD.push_back(gene);
                if (geneInfoMapLD.insert(pair<string, GeneInfo*>(id, gene)).second == false) {
                    LOGGER.w(0, " Duplicate LD block ID found: \"" + id + "\".");
                    return false;
                }
                snpCount = 1;
                continue;  // snpCount ++ not works
            }
            snpCount++;
        } // end of non-dup
    } //  end of while loop
    in.close();
    LOGGER << geneInfoVecLD.size() << " genes to be included from [" + title + "]." << endl;
    return true;

}
bool Data::readMultiEigenMatSnpInfoFile(const string title,vector<EqtlInfo *> &eqtlInfoVecLD, map<string, EqtlInfo *> &eqtlInfoMapLD){
    ifstream in(title.c_str());
    if (!in){ 
        LOGGER.w(0, " can not open the file [" + title + "] to read.");
        return false;
    }
    LOGGER << "Reading eigen snp info from file [" + title + "]." << endl;
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
        eqtl->ld_n  = ld_n;
        eqtlInfoVecLD.push_back(eqtl);
        if (eqtlInfoMapLD.insert(pair<string, EqtlInfo*>(id, eqtl)).second == false) {
            LOGGER.w(0, " Duplicate SNP ID found: \"" + id + "\".");
            return false;
        }
    }
    in.close();
    LOGGER << eqtlInfoVecLD.size() << " SNPs to be included from [" + title + "]." << endl;
    return true;
}
bool Data::readMultiEigenMatBinFile(const string title, vector<GeneInfo*> &geneInfoVecLD,map<string, GeneInfo *> &geneInfoMapLD,vector<EqtlInfo*> &eqtlInfoVecLD,map<string, EqtlInfo *> eqtlInfoMapLD, float eigenCutoff){
    FILE *fp = fopen(title.c_str(), "rb");
    if(!fp){
        LOGGER.e(0, " can not open the file [" + title + "] to read.");
    }
    LOGGER << "Reading eigen binary info from file [" + title + "]." << endl;
    GeneInfo * gene;  // gene from eqtl summary statistics
    GeneInfo * geneLD;  // gene from per gene ld matrix
    map<string, GeneInfo*>::iterator iterGene;
    map<string, SnpInfo*>::iterator iterSnp;
    map<string, EqtlInfo*>::iterator iterEqtl;
    // map<string,MatrixXd> geneEigenVecMap;
    // map<string,VectorXd> geneEigenValMap;
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
            LOGGER.e(0,"Read " + title + " error (m)");
        }
        if(cur_mknum != numSnpInRegion){
            LOGGER.e(0,"In region  " + to_string(i) + ", inconsistent marker number to marker information in " + title);
        }
        // 2. ncol of eigenVec (number of eigenvalues)
        if(fread(&cur_k, sizeof(int32_t), 1, fp) != 1){
            LOGGER.e(0,"In region  " + to_string(i) + ", error about number of eigenvalues in  " + title);
        }
        // 3. sum of all positive eigenvalues
        if(fread(&sumPosEigVal, sizeof(float), 1, fp) != 1){
            LOGGER.e(0,"In region  " + to_string(i) + ", error about sumPosEigVal in " + title);
        }
        // 4. eigenCutoff
        if(fread(&oldEigenCutoff, sizeof(float), 1, fp) != 1){
            LOGGER.e(0,"In region  " + to_string(i) + ", error about oldEigenCutoff used in " + title);
        }
        // 5. eigenvalues
        VectorXf lambda(cur_k);
        if(fread(lambda.data(), sizeof(float), cur_k, fp) != cur_k){
            LOGGER.e(0,"In region  " + to_string(i) + ",size error about eigenvalues in " + title);
        }
        //6. read eigenvector
        MatrixXf U(cur_mknum, cur_k);
        uint64_t nElements = (uint64_t)cur_mknum * (uint64_t)cur_k;
        if(fread(U.data(), sizeof(float), nElements, fp) != nElements){
            LOGGER << "fread(U.data(), sizeof(double), nElements, fp): " << fread(U.data(), sizeof(float), nElements, fp) << endl;
            LOGGER << "nEle: " << nElements << " U.size: " << U.size() <<  " U.col: " << U.cols() << " row: " << U.rows() << endl;
            LOGGER.e(0,"In region  " + to_string(i) + ",size error about eigenvectors in " + title);
        }
        // 8. select eigenvectors based on proportion
        MatrixXd eigenVector;
        VectorXd eigenValue;
        bool haveValue = false;
        int revIdx = 0;
        if(oldEigenCutoff != eigenCutoff & i == 0){
            LOGGER << "Warning: current proportion of variance in region  is set as " + to_string(eigenCutoff)+ ". But the proportion of variance is set as "<< to_string(oldEigenCutoff) + " in "  + title + ".\n";
            // return true;
        }
        eigenValue = lambda.cast<double>();
        eigenVector = U.cast<double>();
        //////////////////////////////////////////////////////////////////////////////////////////
        //////// Step 1.4 build build map between gene and its eigen  matrix of gene ld matrix  /////
        //////////////////////////////////////////////////////////////////////////////////////////
        geneEigenVecMap.insert(pair<string,MatrixXd>(geneLD->ensemblID,eigenVector));
        geneEigenValMap.insert(pair<string,VectorXd>(geneLD->ensemblID,eigenValue));
        geneEigenSumPosEigValMap.insert(pair<string,float>(geneLD->ensemblID,sumPosEigVal));
        geneEigenCutoffMap.insert(pair<string,float>(geneLD->ensemblID,oldEigenCutoff));
    }
    LOGGER << geneInfoVecLD.size()  << " genes are matched between gene eigen  LD reference and eQTL summary data." << endl;
    return true;

}

void Data::MergeMultiBinEigen(const string filename,const string LDMatType){
        string filenameFull = filename + ".eigen." + LDMatType + ".bin";
        FILE *out3 = fopen(filenameFull.c_str(), "wb");
        map<string,MatrixXd>::iterator iterEigenVec;
        map<string,VectorXd>::iterator iterEigenVal;
        map<string,float>::iterator iterEigenSum;
        map<string,float>::iterator iterEigenCutoff;
        auto startTime = std::chrono::steady_clock::now();
        for (unsigned i = 0; i < numKeptGenes; i++) {
            GeneInfo *gene = keptGeneInfoVec[i]; 
            // LOGGER  << " Generate and save eigen of LD matrix from gene region " << i + 1 << "/" << numKeptGenes << "\r" << flush;  
            string taksName = "Merge LD matrices from multiple region";
            Gadget::showProgressBar(i, numKeptGenes, startTime,taksName);
            MatrixXf eigenVec;
            VectorXf eigenVal;
            float sumPosEigVal = 0;
            float eigenCutoff = 0;
            vector<string> snplists = gene->cisSnpNameVec;

            // float eigenCutoff = geneEigenCutoffMap[gene->ensemblID];
            // sumPosEigVal = geneEigenSumPosEigValMap[gene->ensemblID];
            // eigenVal = geneEigenValMap[gene->ensemblID].cast<float>();
            // eigenVec = geneEigenVecMap[gene->ensemblID].cast<float>();
            if(geneEigenCutoffMap.find(gene->ensemblID) != geneEigenCutoffMap.end()){
                eigenCutoff = geneEigenCutoffMap.find(gene->ensemblID)->second;
            } else {
                LOGGER.e(0, "Can not find gene Eigen cutoff vlaue from gene [" + gene->ensemblID + " ].");
            }
            if(geneEigenSumPosEigValMap.find(gene->ensemblID) != geneEigenSumPosEigValMap.end()){
                sumPosEigVal = geneEigenSumPosEigValMap.find(gene->ensemblID)->second;
            } else {
                LOGGER.e(0, "Can not find gene sumPosEigVal vlaue from gene [" + gene->ensemblID + " ].");
            }

            if(geneEigenValMap.find(gene->ensemblID) != geneEigenValMap.end()){
                eigenVal = geneEigenValMap.find(gene->ensemblID)->second.cast<float>();
            } else {
                LOGGER.e(0, "Can not find gene eigenVal vlaue from gene [" + gene->ensemblID + " ].");
            }

            if(geneEigenVecMap.find(gene->ensemblID) != geneEigenVecMap.end()){
                eigenVec = geneEigenVecMap.find(gene->ensemblID)->second.cast<float>();
            } else {
                LOGGER.e(0, "Can not find gene eigenVec vlaue from gene [" + gene->ensemblID + " ].");
            }


            // save eigen  matrix
            int32_t numEigenValue = eigenVal.size();
            // int32_t numSnpInBlock = ldblock->snpNameVec.size();
            int32_t numSnpInBlock = snplists.size();
            double eigenValueSum = eigenVal.sum();
            // 1, nrow of eigenVecGene[i] 
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
        }
        fclose(out3);
}

void Data::readEigenMatLDBlockInfoFile(const string &infoFile){
       ifstream in(infoFile.c_str());
    if (!in) LOGGER.e(0, " can not open the file [" + infoFile + "] to read.");
    LOGGER << "Reading LD Block LDM block info from file [" + infoFile + "]." << endl;
    map<string,int> snp2snpMap;
    string header;
    string id;
    int  chr, start, end, windows, snpNum;
    int idx = 0;
    int snpCount =  1;
    string snpID;
    LDBlockInfo *ldblock;
    SnpInfo * snp;
    // Step 1. read gene info from ld gene matrix. In this step, if gene in LD file does not belong to 
    // gene set in BESD file. gene->kept will be set as false;
    getline(in, header);
    while (in >>chr>>id >>start>>end >> snpID >> snpNum ) {
        if (snp2snpMap.insert(pair<string, int>(id + "_" + snpID , snpCount)).second == false) {
            LOGGER.e(0, " Duplicate LDBlock-SNP pair found: \"" + id + "_" + snpID + "\".");
        } else{
            // check whether snp can be found in ld snp info 
                if(snpInfoMap.find(snpID) == snpInfoMap.end()){
                    LOGGER.e(0,"Cannot find snp [ "+ snpID + " ].");
                }else {
                    snp = snpInfoMap.find(snpID)->second;
                    snp->block = id;
                }
            // Step 1.1 build a new geneInfo firstly
            if(snpCount == 1){
                ldblock = new LDBlockInfo(idx++, id, chr);
                ldblock->startPos = start;
                ldblock->endPos = end;
            }
            // Step 1.3 continue to add eQTLs into gene
            if(snpCount < snpNum){
                ldblock->snpNameVec.push_back(snpID);  
                ldblock->snpInfoVec.push_back(snp);   
            }
            // Step 1.5 add the last eQTL into gene 
            if(snpCount == snpNum){
                ldblock->snpNameVec.push_back(snpID);
                ldblock->snpInfoVec.push_back(snp);
                ldblock->numSnpInBlock = snpNum;

                ldBlockInfoVec.push_back(ldblock);
                if (ldBlockInfoMap.insert(pair<string, LDBlockInfo*>(id, ldblock)).second == false) {
                    LOGGER.e(0, " Duplicate LD block ID found: \"" + id + "\".");}
                snpCount = 1;
                continue;  // snpCount ++ not works
            }
            snpCount++;
        } // end of non-dup
    } //  end of while loop
    in.close();
    numLDBlocks = ldBlockInfoVec.size();
    LOGGER << numLDBlocks << " LD blocks to be included from [" + infoFile + "]." << endl;
};
void Data::readEigenMatLDBlockSnpInfoFile(const string &snpInfoFile,const string LDMatType){
        ifstream in(snpInfoFile.c_str());
    if (!in) LOGGER.e(0, " can not open the file [" + snpInfoFile + "] to read.");
    LOGGER << "Reading LD Block LDM SNP info from file [" + snpInfoFile + "]." << endl;
    string header;
    string id, allele1, allele2;
    int chr, physPos,ld_n;
    double genPos;
    double allele1Freq;
    int idx = 0;
    SnpInfo *snp;
    snpInfoMap.clear();
    snpInfoVec.clear();
    chromosomes.clear();
    getline(in, header);
    int imputeAf = 0;
    while (in >> chr >> id >> genPos >> physPos >> allele1 >> allele2>>allele1Freq>>ld_n) {
        snp = new SnpInfo(idx++, id, allele1, allele2, chr, genPos, physPos);
        snp->ld_n = ld_n;
        snp->af = allele1Freq;
        snpInfoVec.push_back(snp);
        chromosomes.insert(snp->chrom);
        if (snpInfoMap.insert(pair<string, SnpInfo*>(id, snp)).second == false) {
            LOGGER.e(0, " Duplicate SNP ID found: \"" + id + "\".");
        }
    }
    in.close();
    numSnps = snpInfoVec.size();
    LOGGER << numSnps << " SNPs to be included from [" + snpInfoFile + "]." << endl;
};

///////////// Step 2.3 read LD matrix of gene regions info (3 functions)
void Data::readEigenMatGeneInfoFile( const string &infoFile){
    ifstream in(infoFile.c_str());
    if (!in) LOGGER.e(0, " can not open the file [" + infoFile + "] to read.");
    LOGGER << "Reading xQTL LDM gene info from file [" + infoFile + "]." << endl;
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
                gene->start = start;
                gene->end   = end;
                gene->cisSnpNameVec.clear();
                gene->gene2CisSnpVec.clear();
                gene->cisSnpIDSetInGene.clear();
                gene->cisSnpID2IdxMapInGene.clear();
            }
            // Step 1.2 find wheter given gene eigen snp exists in GWAS ld reference or not
            if(snpInfoMap.find(snp) == snpInfoMap.end()){ 
                // if not, we need to remove this gene to make sure the consistence between GWAS LD and gene LD.
                LOGGER.w(0,"gene " + gene->ensemblID + " was removed since eQTL " + snp + " from gene can not be found from GWAS LD reference.");
                gene->kept = false;
            } else {
                // Found, but will be removed since snp->included is false due to 1) --keep-block
                // here I use eqtl class since in the readEigenMatGeneSnpInfoFile,
                // we already pass snp->included value to eqtl class
                iterEqtl = eqtlInfoMap.find(snp);
                if(!iterEqtl->second->included) gene->kept = false;
            }
            // Step 1.3 continue to add eQTLs into gene
            if(snpCount < snpNum){
                gene->cisSnpNameVec.push_back(snp);
                gene->gene2CisSnpVec.push_back(snpCount -1 );
                gene->cisSnpIDSetInGene.insert(snp);
                gene->cisSnpID2IdxMapInGene.insert(pair<string, int>(snp,snpCount -1 ));
                
            }
            // Step 1.5 add the last eQTL into gene 
            if(snpCount == snpNum){
                gene->cisSnpNameVec.push_back(snp);
                gene->gene2CisSnpVec.push_back(snpCount -1 );
                gene->cisSnpIDSetInGene.insert(snp);
                gene->cisSnpID2IdxMapInGene.insert(pair<string, int>(snp,snpCount -1 ));

                gene->numSnpInGene = snpNum;

                geneInfoVec.push_back(gene);
                if (geneInfoMap.insert(pair<string, GeneInfo*>(id, gene)).second == false) {
                    LOGGER.e(0, " Duplicate LD block ID found: \"" + id + "\".");}
                snpCount = 1;
                continue;  // snpCount ++ not works
            }
            snpCount++;
        } // end of non-dup
    } //  end of while loop
    in.close();
    numGenes = geneInfoVec.size();
    LOGGER << numGenes << " genes to be included from [" + infoFile + "]." << endl;
}
void Data::readEigenMatGeneSnpInfoFile(const string &snpInfoFile,const unsigned includeChr){
    ifstream in(snpInfoFile.c_str());
    if (!in) LOGGER.e(0, " can not open the file [" + snpInfoFile + "] to read.");
    LOGGER << "Reading xQTL LDM SNP info from file [" + snpInfoFile + "]." << endl;
    EqtlInfo *eqtl;
    map<string, EqtlInfo*>::iterator iterEqtl;
    map<string, SnpInfo*>::iterator iterSnp;
    SnpInfo * snp;
    string header;
    string id, allele1, allele2;
    int chr, physPos,ld_n;
    double genPos;
    double allele1Freq;
    int idx = 0;
    getline(in, header);
    int imputeAf = 0;
    while (in >> chr >> id >> genPos >> physPos >> allele1 >> allele2>>allele1Freq>>ld_n) {
        if(chr != includeChr && includeChr != 0) continue;
        eqtl = new EqtlInfo(idx++, id, allele1, allele2, chr, genPos, physPos,allele1Freq);
        eqtl->ld_n  = ld_n;  // here only LD sample size
        eqtlInfoVec.push_back(eqtl);
        eqtl->af = allele1Freq;
        if (eqtlInfoMap.insert(pair<string, EqtlInfo*>(id, eqtl)).second == false) {
            LOGGER.e(0, " Duplicate SNP ID found: \"" + id + "\".");
        }
        // check snp->included
        
        // if eqtl doesn't exist compared with gwas LD, this gene will be removed in  function.
        iterSnp = snpInfoMap.find(id);
        if(iterSnp != snpInfoMap.end()){ //found
            // check allele flip. and report an error is allele flip occurs.
            snp = iterSnp->second;
            if(!snp->included) eqtl->included = false; // this will be worked when using --flag;
            if (allele1 != iterSnp->second->a1 && allele2 != iterSnp->second->a2) {
                LOGGER.e(0,"The allele1 (allele2) strand of [ " + id +  " ] is different between GWAS LD and xQTL LD \n " +
                           "Please regenerate GWAS LD and xQTL LD using same LD reference genotypes.");
            }
        } else { // not found
            LOGGER.e(0,"The eqtl [" + id +  " ] cannot be found between GWAS LD and xQTL LD \n " +
                "Please re-generate GWAS LD and xQTL LD using same LD reference genotypes.");
        }
    }
    in.close();
    numeQTLs = eqtlInfoVec.size();
    LOGGER << numeQTLs << " SNPs to be included from [" + snpInfoFile + "]." << endl;
}
void Data::readEigenMatBinFile(const string &geneEigenMatrixFile, const double eigenCutoff,const string LDMatType,const double diag_mod){
    string filename = geneEigenMatrixFile + ".eigen." + LDMatType + ".bin" ; 
    FILE *fp = fopen(filename.c_str(), "rb");
    if(!fp){LOGGER.e(0, " can not open the file [" + filename + "] to read.");}
    
    GeneInfo * gene;  // gene from eqtl summary statistics
    LDBlockInfo * ldblock;
    vector<int>numSnpInRegion;
    int numRegion = 0;
    vector<bool> regionKeptList;

    if(LDMatType == "ldblock"){
        LOGGER << "Reading LD Block LDM binary info from file [" + filename + "]." << endl;
        numSnpInRegion.resize(numLDBlocks);
        regionKeptList.resize(numLDBlocks);
        for(int i = 0; i < numLDBlocks;i++){
            ldblock = ldBlockInfoVec[i];
            numSnpInRegion[i] = ldblock->numSnpInBlock;
            regionKeptList[i] = ldblock->kept;
        }
        eigenValLdBlock.resize(numLDBlocks);
        eigenVecLdBlock.resize(numLDBlocks);
    } else if(LDMatType == "gene"){
        LOGGER << "Reading xQTL LDM binary info from file [" + filename + "]." << endl;
        numSnpInRegion.resize(numGenes);
        regionKeptList.resize(numGenes);
        for(int i = 0; i < numGenes; i++){
            gene= geneInfoVec[i];
            numSnpInRegion[i] = gene->numSnpInGene;
            regionKeptList[i] = gene->kept;
        }
        eigenValGene.resize(numGenes); 
        eigenVecGene.resize(numGenes); 
    }
    int eigenIdx = 0;
    numRegion = numSnpInRegion.size();
    auto startTime = std::chrono::steady_clock::now();
    for(int i = 0; i < numRegion; i++){
        Gadget::showProgressBar(i, numRegion, startTime,"Read GWAS LD blocks");
        int32_t cur_mknum = 0;
        int32_t cur_k = 0;
        float sumPosEigVal = 0;
        float oldEigenCutoff =0;
        // 1. marker number
        if(fread(&cur_mknum, sizeof(int32_t), 1, fp) != 1){
            LOGGER.e(0,"Read " + filename + " error (m)");
        }
        if(cur_mknum != numSnpInRegion[i]){
            LOGGER << "Current number: " << cur_mknum << " numSnpInRegion: " << numSnpInRegion[i] << endl;
            LOGGER.e(0,"In region  " + to_string(i) + ", inconsistent marker number to marker information in " + filename);
        }
        // 2. ncol of eigenVec (number of eigenvalues)
        if(fread(&cur_k, sizeof(int32_t), 1, fp) != 1){
            LOGGER.e(0,"In region  " + to_string(i) + ", error about number of eigenvalues in  " + filename);
        }
        // 3. sum of all positive eigenvalues
        if(fread(&sumPosEigVal, sizeof(float), 1, fp) != 1){
            LOGGER.e(0,"In region  " + to_string(i) + ", error about sumPosEigVal in " + filename);
        }
        // 4. eigenCutoff
        if(fread(&oldEigenCutoff, sizeof(float), 1, fp) != 1){
            LOGGER.e(0,"In region  " + to_string(i) + ", error about oldEigenCutoff used in " + filename);
        }
        // 5. eigenvalues
        VectorXf lambda(cur_k);
        if(fread(lambda.data(), sizeof(float), cur_k, fp) != cur_k){
            LOGGER.e(0,"In region  " + to_string(i) + ",size error about eigenvalues in " + filename);
        }
        //6. read eigenvector
        MatrixXf U(cur_mknum, cur_k);
        uint64_t nElements = (uint64_t)cur_mknum * (uint64_t)cur_k;
        if(fread(U.data(), sizeof(float), nElements, fp) != nElements){
            LOGGER << "fread(U.data(), sizeof(double), nElements, fp): " << fread(U.data(), sizeof(float), nElements, fp) << endl;
            LOGGER << "nEle: " << nElements << " U.size: " << U.size() <<  " U.col: " << U.cols() << " row: " << U.rows() << endl;
            LOGGER.e(0,"In region  " + to_string(i) + ",size error about eigenvectors in " + filename);
        }
        // 8. select eigenvectors based on proportion
        MatrixXd eigenVector;
        VectorXd eigenValue;
        bool haveValue = false;
        int revIdx = 0;
        if(oldEigenCutoff != eigenCutoff & i == 0){
            LOGGER.w(0,"Current proportion of variance in region is set as " + to_string(eigenCutoff)+ ". But the proportion of variance is set as " + to_string(oldEigenCutoff) + " in "  + filename + ".\n");
        }
        if(eigenCutoff < oldEigenCutoff){
            truncateEigenMatrix(sumPosEigVal, eigenCutoff, lambda.cast<double>(), U.cast<double>(), eigenValue, eigenVector);
        } else {
            eigenValue = lambda.cast<double>();
            eigenVector = U.cast<double>();
        }

        if(LDMatType == "ldblock"){
            if(regionKeptList[i] ){
                eigenValLdBlock[eigenIdx] = eigenValue;
                eigenVecLdBlock[eigenIdx] = eigenVector;
                eigenIdx++;
                if(eigenIdx == numKeptLDBlocks) break;
            }
        } else if(LDMatType == "gene"){
            if(regionKeptList[i]){
                eigenValGene[eigenIdx] = eigenValue ; 
                eigenVecGene[eigenIdx] = eigenVector; 
                eigenIdx++;
                if(eigenIdx == numKeptGenes) break;
            }
        } 
    }
    if(LDMatType == "ldblock"){
        LOGGER << numLDBlocks << " LD Block LDM bin to be included from [" + filename + "]." << endl;
    } else if(LDMatType == "gene"){
        LOGGER << numKeptGenes << " gene LDM bin to be included from [" + filename + "]." << endl;
    } 
}


void Data::readSBayesRCBlockLdmInfoFile(const string &infoFile){
    // Read bim file: recombination rate is defined between SNP i and SNP i-1
    ifstream in(infoFile.c_str());
    if (!in) LOGGER.e(0,"can not open the file [" + infoFile + "] to read.");
    LOGGER << "Reading GWAS LDM info from file [" + infoFile + "]." << endl;
    ldBlockInfoVec.clear();
    ldBlockInfoMap.clear();
        
    string header;
    string id;
    int  chr, blockStart, blockEnd, snpNum;
    int idx = 0;
    int snpCount =  1;
    string snpName;
    string startSnpID, endSnpID;
    LDBlockInfo *ldblock;
    getline(in, header);
    while (in >> id >> chr >> blockStart >> startSnpID >> blockEnd >> endSnpID >> snpNum) {
        
        ldblock = new LDBlockInfo(idx++, id, chr);
        ldblock->startSnpIdx = blockStart;
        ldblock->endSnpIdx   = blockEnd;
        ldblock->numSnpInBlock = snpNum;
        
        ldBlockInfoVec.push_back(ldblock);
        ldBlockInfoMap[id] = ldblock;

    }
    in.close();
    numLDBlocks = (unsigned) ldBlockInfoVec.size();
    LOGGER << numLDBlocks << " GWAS LD Blocks to be included from [" + infoFile + "]." << endl;
}

void Data::readSBayesRCBlockLdmSnpInfoFile(const string &snpInfoFile){
    ifstream in(snpInfoFile.c_str());
    if (!in) LOGGER.e(0,"can not open the file [" + snpInfoFile + "] to read.");
    LOGGER << "Reading GWAS LDM SNP info from file [" + snpInfoFile + "]." << endl;
    snpInfoVec.clear();
    snpInfoMap.clear();
    map<string,int> ld2snpMap;
    string header;
    string id, allele1, allele2;
    int chr, physPos,ld_n;
    double genPos;
    double allele1Freq;
    int idx = 0;
    int index;
    string blockID;
    LDBlockInfo *ldblock;
    getline(in, header);
    while (in >> chr >> id >> index >> genPos >> physPos >> allele1 >> allele2 >> allele1Freq >> ld_n >> blockID) {
        SnpInfo *snp = new SnpInfo(idx++, id, allele1, allele2, chr, genPos, physPos);
        
//        LOGGER << idx << " " << id << " " << blockID << endl;
        
        snp->af = allele1Freq;
        snp->ld_n = ld_n;
        snp->block = blockID;
        snpInfoVec.push_back(snp);
        chromosomes.insert(snp->chrom);
        
        ldblock = ldBlockInfoMap[blockID];
        
        if (ld2snpMap.insert(pair<string, int>(blockID + "_" + id, idx)).second == false) {
            LOGGER.e(0,"Duplicate LDBlock-SNP pair found: \"" + blockID + "_" + id + "\".");
        } else{
            ldblock->snpNameVec.push_back(id);
            ldblock->snpInfoVec.push_back(snp);
        }

        if (snpInfoMap.insert(pair<string, SnpInfo*>(id, snp)).second == false) {
            LOGGER.e(0,"Duplicate SNP ID found: \"" + id + "\".");
        }
    }
    in.close();
    numSnps = (unsigned) snpInfoVec.size();
    LOGGER << numSnps << " SNPs to be included from [" + snpInfoFile + "]." << endl;
}
