#include "Omics.hpp"

void Omics::covertUkbPppPGWASToBesdFormat(Data &data, const string rsidMapFiles, const string proteinMapFile, const string proteinPath,
    const bool &makeBesdBool,const bool &makeBesdSmrBool,const bool &makeQueryBool,const bool &makeSumstatsFormatBool,const bool &makeQueryCojoMaBool,
    const double pValueThreshold, const string grchType, const string title){
    data.readSNPRSIDMaps(rsidMapFiles,grchType,title);
    data.readProMaps(proteinMapFile,pValueThreshold,title);
    // data.saveBesdFormat(title,false);
    if(makeBesdBool) data.saveBesdFormat(title,false); // add averaged gene sample size columns to epi file
    if(makeBesdSmrBool) data.saveBesdFormat(title,makeBesdSmrBool); // formal smr format, no additional N column in the esi file
    if(makeQueryBool) data.saveQueryMoQTLInfo(title); // output all gene and eqtl info in a gz plain format.
    // data.saveQueryMoQTLInfo(title);
}

void Omics::convertFlistFormat2BESDFormat(Data &data,const string title, const string geneInfoFile, const unsigned includeChr,
        const bool makeBesdBool,const bool makeBesdSmrBool, const bool makeQueryBool,const string subGenePath){
    if(geneInfoFile.empty()){
        LOGGER.e(0,"Can not find [ " + geneInfoFile + " ].");
    }

    if(!geneInfoFile.empty()) data.readFlistFile(geneInfoFile,subGenePath);
    if(makeBesdBool) data.saveBesdFormat(title,false); // add averaged gene sample size columns to epi file
    if(makeBesdSmrBool) data.saveBesdFormat(title,makeBesdSmrBool); // formal smr format, no additional N column in the esi file
    if(makeQueryBool) data.saveQueryMoQTLInfo(title); // output all gene and eqtl info in a gz plain format.
}

void Omics::convertBesdQuery2OtherFormat(Data &data,const string &eigenMatrixFile,
        const string &besdFile,const string &eqtlSummaryQueryFile,const string &geneSamSizeFile, const string &title,const string &includeGeneFile, const string excludeGeneFile,
        const string &specificGeneID, const string &includeSnpFile,const unsigned includeChr,const double &pValueThreshold, 
        const double &cisRegionWind, const bool &useCisWinBool, const int &slurmArrayLimit,const bool &makeBesdBool,
        const bool &makeBesdSmrBool,const bool &makeQueryBool,const bool &makeSumstatsFormatBool,
        const bool &makeQueryCojoMaBool, const bool &makeAnnoBool,const bool &isAnnoBinary, 
        const double &eigenCutoff){
    /////////////////////////////////////////////////////////////
    //// Step 1. Read necessary information
    /////////////////////////////////////////////////////////////
    /// Step 1.1 GWAS LD matrix info, in the integrative analysis, we need to make sure
    /// SNPs in xQTLs match that in GWAS LD matrix
    bool hasSnpInfo = false; // 
    if(eigenMatrixFile.size()!= 0){
        bool qcBool = true;
        data.readEigenMatrix(eigenMatrixFile,eigenCutoff,qcBool);
        data.includeMatchedSnp();
        hasSnpInfo = true;
    }
    /// Step 1.2  read snp and gene info from BESD format and then
    /// we need to read binary besd file to extract beta and se after 
    /// QC step
    bool smrBesdBool = false;
    bool makeGeneLDEigenBool  = false; // used in besd file
    bool hasGeneLdmInfo = false;
    if(!besdFile.empty()){
        data.readEsiFile(besdFile + ".esi",smrBesdBool,false);
        data.readEpiFile(besdFile + ".epi");
    }

    /// Step 1.3 read query format, now all information will be added,
    /// but we need to further QC
    if(! eqtlSummaryQueryFile.empty()) {
        data.geneInfoVec.clear(); data.geneInfoMap.clear();
        data.eqtlInfoVec.clear(); data.eqtlInfoMap.clear();
        data.readQeuryGZFormat(eqtlSummaryQueryFile + ".query.gz",data.geneInfoVec,data.geneInfoMap,data.eqtlInfoVec,data.eqtlInfoMap,smrBesdBool,hasSnpInfo,hasGeneLdmInfo);
        // summary data first becasue above function don't have summary data
        data.numeQTLs = data.eqtlInfoVec.size();
        data.numGenes = data.geneInfoVec.size();
        data.keptGeneInfoVec  = data.makeKeptGeneInfoVec(data.geneInfoVec);
        data.numKeptGenes = (unsigned) data.keptGeneInfoVec.size();
        data.incdEqtlInfoVec = data.makeIncdEqtlInfoVec(data.eqtlInfoVec);
        data.numIncdEqtls = (unsigned) data.incdEqtlInfoVec.size();
    }
    /////////////////////////////////////////////////////////////
    //// Step 2. QC: Select xQTL subsets based on various options
    /////////////////////////////////////////////////////////////
    // Specific chromsome
    if (includeChr) data.includeChreQTL(includeChr);
    // Specific genes, extract or exclude genes
    if(!includeGeneFile.empty()) data.includeGene(includeGeneFile);
    if(!excludeGeneFile.empty()) data.excludeGene(excludeGeneFile);
    // here we only use one gene.
    if(specificGeneID.size()) data.includeSpecificGene(specificGeneID);
    // Specific snps
    if(!includeSnpFile.empty()) data.includexQTLSnp(includeSnpFile);
    // data.includeMatchedEqtl(); // construct chrInfoMap;
    /// Now we select SNPs based on BESD format or query format
    // read besd info and select eqtl based on p-value,
    // we don't need to do this in the query format
    if(!besdFile.empty()) data.readBesdFile(besdFile + ".besd",hasSnpInfo,makeGeneLDEigenBool,smrBesdBool); // 
    // match select snplist based on cis-windows;
    if(useCisWinBool)data.filterMoQTLBasedOnCisWindow(cisRegionWind);
    if(pValueThreshold < 1) data.filterMoQTLBasedOnPvalue(pValueThreshold);
    /////////////////////////////////////////////////////////////
    /// Step 3. additional information
    // if a snp has different per-SNP sample size for different genes, we need to do here
    if(!geneSamSizeFile.empty()){
        data.readGeneSampleSizeFile(geneSamSizeFile);
        // QC use per-SNP-sample size
        // data.
    }
    /////////////////////////////
    /// summary
    data.keptGeneInfoVec  = data.makeKeptGeneInfoVec(data.geneInfoVec);
    data.numKeptGenes = (unsigned) data.keptGeneInfoVec.size();
    data.incdEqtlInfoVec = data.makeIncdEqtlInfoVec(data.eqtlInfoVec);
    data.numIncdEqtls = (unsigned) data.incdEqtlInfoVec.size();
    //////////////////////////////
    // output various formats
    /////////////////////////////
    // please note that when saving xqtl data into specific format, three parameters are
    // quite important, gene->cisSnpNameVec, gene->cisSnpIDSetInGene,
    // gene->cisSnpID2IdxMapInGene
    if(makeAnnoBool) data.saveAnnoPlainMatFormat(title, isAnnoBinary,hasSnpInfo); // 
    if(makeBesdBool) data.saveBesdFormat(title,false); // add sample size columns to esi file
    if(makeBesdSmrBool) data.saveBesdFormat(title,makeBesdSmrBool); // formal smr format, no additional N column in the esi file
    if(makeQueryBool) data.saveQueryMoQTLInfo(title); // output all gene and eqtl info in a gz plain format.
    if(makeQueryCojoMaBool) data.saveQueryCojoMaMoQTLInfo(title,slurmArrayLimit);
    if(makeSumstatsFormatBool) {data.saveSumstatsFormat(title);}
}


void Omics::resamplexQTLSampleSize(Data &data,const string eigenMatrixFile,const string geneEigenMatrixFile,
    const string besdFile,const string title, const string &includeGeneFile,const string specificGeneID,
    const int reSamType, const bool makeQueryCojoMaBool,
    const bool makeXqtlNBool, const int reEffectSampleSize,const int slurmArrayLimit, const double eigenCutoff){
    bool hasSnpInfo = false; // 
    // cout << "eigen: " << eigenMatrixFile.size() << endl;
    bool qcBool = true;
    data.readEigenMatrix(eigenMatrixFile,eigenCutoff,qcBool);
    // Select genes from files
    if(!includeGeneFile.empty()) data.includeGene(includeGeneFile);
    // here we only use one gene.
    if(specificGeneID.size()) data.includeSpecificGene(specificGeneID);
    data.includeMatchedSnp();
    hasSnpInfo = true;
    // here we read xQTL LD matrix and then read besd format, please note, we need to regenerate 
    // xQTL LD and besd format with same eQTL snplist.
    data.reSampleXqtlEffect(reSamType,besdFile,geneEigenMatrixFile,reEffectSampleSize,eigenCutoff);
    //////////////////////////////
    // output various formats
    /////////////////////////////
    if(makeQueryCojoMaBool) data.saveQueryCojoMaMoQTLInfo(title,slurmArrayLimit);
}


void Omics::inputIndInfo(Data &data, const string &bedFile, const string &phenotypeFile, const string &keepIndFile, 
        const unsigned keepIndMax,const unsigned mphen, const string &covariateFile, const string &randomCovariateFile,
        const string &residualDiagFile,const string &geneInfoFile, const string &keepIndGeneFile,const string subGenePath,bool ldBool){
    data.readFamFile(bedFile + ".fam");
    data.readPhenotypeFile(phenotypeFile, mphen);
    data.readCovariateFile(covariateFile);
    data.readRandomCovariateFile(randomCovariateFile);
    data.readResidualDiagFile(residualDiagFile);
    data.keepMatchedInd(keepIndFile, keepIndMax,ldBool);
    if(! geneInfoFile.empty()) {
        data.readGeneInfoFile(geneInfoFile,subGenePath); // here we need to read phenotypes from molecular traits
        data.keepMatchedIndGene(keepIndGeneFile, keepIndMax);
    }
}

// individual genotype
void Omics::inputSnpInfo(Data &data, const string &bedFile, const string &geneInfoFile,const double cisRegionWind,
                        const string &includeSnpFile, const string &excludeSnpFile, const string &excludeRegionFile, 
                        const unsigned includeChr, const bool excludeAmbiguousSNP, const string &skeletonSnpFile, 
                        const string &geneticMapFile,  const string &ldBlockInfoFile, const unsigned includeBlock, 
                        const unsigned flank, const double mafmin, const double mafmax, const bool noscale, 
                        const bool readGenotypes){
    data.readBimFile(bedFile + ".bim");
    if (!includeSnpFile.empty()) data.includeSnp(includeSnpFile);
    if (!excludeSnpFile.empty()) data.excludeSnp(excludeSnpFile);
    if (includeChr) data.includeChr(includeChr);
    if (excludeAmbiguousSNP) data.excludeAmbiguousSNP();
    if (!excludeRegionFile.empty()) data.excludeRegion(excludeRegionFile);
    if (!skeletonSnpFile.empty()) data.includeSkeletonSnp(skeletonSnpFile);
    if (!geneticMapFile.empty()) data.readGeneticMapFile(geneticMapFile);
    if (!ldBlockInfoFile.empty()) data.readLDBlockInfoFile(ldBlockInfoFile);
    data.includeMatchedSnp();
    if (includeBlock) data.includeBlock(includeBlock);
    if (readGenotypes) data.readBedFile(noscale, bedFile + ".bed");
    // data.includeMatchedSnp need to be used first!
    if(!geneInfoFile.empty()){
        data.mapGwasSnpToGeneCisRegion(bedFile + ".bed",noscale,cisRegionWind);
        data.ConstructGenePheAndwAcorr();
    } 
}

void Omics::inputSnpInfo(Data &data, const string &bedFile, const string &gwasSummaryFile, const double afDiff, const double mafmin, const double mafmax, const double pValueThreshold, const bool sampleOverlap, const bool imputeN, const bool noscale){
    data.readFamFile(bedFile + ".fam");
    data.readBimFile(bedFile + ".bim");

    data.keptIndInfoVec = data.makeKeptIndInfoVec(data.indInfoVec);
    data.numKeptInds =  (unsigned) data.keptIndInfoVec.size();
    
    data.readGwasSummaryFile(gwasSummaryFile, afDiff, mafmin, mafmax, pValueThreshold, imputeN, true);
    data.includeMatchedSnp();
    data.readBedFile(noscale, bedFile + ".bed");
    data.buildSparseMME(sampleOverlap, noscale);
}

// this function use original LD matrix as LD reference
void Omics::inputSnpInfo(Data &data, const string &includeSnpFile, const string &excludeSnpFile, const string &excludeRegionFile, const string &gwasSummaryFile, const string &ldmatrixFile, const unsigned includeChr, const bool excludeAmbiguousSNP, const string &skeletonSnpFile, const string &geneticMapFile, const double genMapN, const unsigned flank, const string &eQTLFile, const string &ldscoreFile, const string &windowFile, const bool multiLDmat, const bool excludeMHC, const double afDiff, const double mafmin, const double mafmax, const double pValueThreshold, const double rsqThreshold, const bool sampleOverlap, const bool imputeN, const bool noscale, const bool binSnp, const bool readLDMfromTxtFile){
    if (multiLDmat)
        data.readMultiLDmatInfoFile(ldmatrixFile);
    else
        data.readLDmatrixInfoFile(ldmatrixFile + ".info");
    if (!includeSnpFile.empty()) data.includeSnp(includeSnpFile);
    if (!excludeSnpFile.empty()) data.excludeSnp(excludeSnpFile);
    if (includeChr) data.includeChr(includeChr);
    if (excludeAmbiguousSNP) data.excludeAmbiguousSNP();
    if (!excludeRegionFile.empty()) data.excludeRegion(excludeRegionFile);
    if (excludeMHC) data.excludeMHC();
    if (!skeletonSnpFile.empty()) data.includeSkeletonSnp(skeletonSnpFile);
    if (!geneticMapFile.empty()) data.readGeneticMapFile(geneticMapFile);
    if (!ldscoreFile.empty()) data.readLDscoreFile(ldscoreFile);
    if (!windowFile.empty()) data.readWindowFile(windowFile);
    if (!gwasSummaryFile.empty()) data.readGwasSummaryFile(gwasSummaryFile, afDiff, mafmin, mafmax, pValueThreshold, imputeN, true);
    data.includeMatchedSnp();
    if (readLDMfromTxtFile) {
        data.readLDmatrixTxtFile(ldmatrixFile + ".txt");
    } else {
//        if (geneticMapFile.empty()) {
            if (multiLDmat)
                data.readMultiLDmatBinFile(ldmatrixFile);
            else
                data.readLDmatrixBinFile(ldmatrixFile + ".bin");
//        } else {
//            if (multiLDmat)
//                data.readMultiLDmatBinFileAndShrink(ldmatrixFile, genMapN);
//            else
//                data.readLDmatrixBinFileAndShrink(ldmatrixFile + ".bin");
//        }
    }
    
    if (rsqThreshold < 1.0 && !binSnp) {
        data.filterSnpByLDrsq(rsqThreshold);
        data.includeMatchedSnp();
        if (geneticMapFile.empty()) {  // need to read LD data again after LD filtering
            if (multiLDmat)
                data.readMultiLDmatBinFile(ldmatrixFile);
            else
                data.readLDmatrixBinFile(ldmatrixFile + ".bin");
        } else {
            if (multiLDmat)
                data.readMultiLDmatBinFileAndShrink(ldmatrixFile, genMapN);
            else
                data.readLDmatrixBinFileAndShrink(ldmatrixFile + ".bin");
        }
    }
    if (!gwasSummaryFile.empty()) data.buildSparseMME(sampleOverlap, noscale);
    if (!windowFile.empty()) data.binSnpByWindowID();
}


void Omics::clearGenotypes(Data &data){
    data.X.resize(0,0);
}

