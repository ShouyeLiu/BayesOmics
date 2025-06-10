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
    } else {
        data.ConstructGwasEqtlGeneMaps();
    }
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

// this function read eigen matrices 
void Omics::inputSnpInfo(Data &data, const string &includeSnpFile, const string &excludeSnpFile, const string &excludeRegionFile,const string &gwasSummaryFile, 
                        const string &eqtlSummaryFile, const string &eqtlSummaryQueryFile, const string &includeGeneFile, const string &geneSamSizeFile, 
                        const string &eigenMatrixFile, const string &geneEigenMatrixFile, const string &ldBlockInfoFile,
                        const unsigned includeChr, const bool excludeAmbiguousSNP, const unsigned flank, const string &eQTLFile, const string &ldscoreFile,
                        const double eigenCutoff, const double geneEigenCutoff, const bool excludeMHC,
                        const double afDiff, const double mafmin, const double mafmax, const double pValueThreshold, const double rsqThreshold,
                        const bool sampleOverlap, const bool imputeN, const bool noscale, const bool readLDMfromTxtFile, const bool imputeSummary, 
                        const unsigned includeBlock, const string includeBlockID){
    // Step 1. read LD Block eigen matrix
    data.readEigenMatrix(eigenMatrixFile, eigenCutoff);
    // Step 2. do QC process for GWAS parts
    if (!includeSnpFile.empty()) data.includeSnp(includeSnpFile); // this may not work in eigen matrix part
    if (!excludeSnpFile.empty()) data.excludeSnp(excludeSnpFile); // this may not work in eigen matrix part
    if (includeChr) data.includeChr(includeChr);
    if (includeBlock) data.includeBlock(includeBlock);
    if (!includeBlockID.empty()) data.includeBlockID(includeBlockID);
    if (excludeAmbiguousSNP) data.excludeAmbiguousSNP();   // this may not work in eigen matrix part
    if (!excludeRegionFile.empty()) data.excludeRegion(excludeRegionFile); // this may not work in eigen matrix part
    if (excludeMHC) data.excludeMHC();  // this may not work in eigen matrix part
    if (!ldscoreFile.empty()) data.readLDscoreFile(ldscoreFile); 
    // Step 3. read gwas summary statistics
    if (!gwasSummaryFile.empty()) {
        bool removeOutlierN = imputeSummary;
        data.readGwasSummaryFile(gwasSummaryFile, afDiff, mafmin, mafmax, pValueThreshold, imputeN, removeOutlierN);
        if (imputeSummary) {
            data.readEigenMatrixBinaryFile(eigenMatrixFile, eigenCutoff);
            data.impG(includeBlock);
            return;
        }
        data.includeMatchedSnp();  // very necessary
        if ( !eqtlSummaryFile.empty() || !eqtlSummaryQueryFile.empty()){
            // read xQTL LD info: baseline model
            data.readEigenMatrixGene(geneEigenMatrixFile,eigenCutoff,"gene",includeChr);
            // Here is the place to do QC. For example, remove genes or xQTLs 
            /// inlcude genes
            if(!includeGeneFile.empty()) data.includeGene(includeGeneFile);
            /// end of QC
            data.readEqtlSummaryFileFromLD2BESD(eqtlSummaryFile, eqtlSummaryQueryFile,includeChr, afDiff, mafmin, mafmax, pValueThreshold, imputeN, geneEigenCutoff);
            // if(!geneSamSizeFile.empty()){data.readGeneSampleSizeFile(geneSamSizeFile);}
        }
    }
    // Step 4. build model parameter.
    if(!gwasSummaryFile.empty()){
        bool gwasInfoOnly = false;
        if(eqtlSummaryFile.empty() && eqtlSummaryQueryFile.empty() ){ 
            gwasInfoOnly = true;
        }
        // here read GWAS LD low-rank matrix 
        data.buildMMEigen(eigenMatrixFile, sampleOverlap, eigenCutoff, noscale,gwasInfoOnly);  
        if(!eqtlSummaryFile.empty() || !eqtlSummaryQueryFile.empty()){
            bool haveEqtlInfo = true;
            data.buildMMGeneEigen(sampleOverlap,geneEigenCutoff, noscale,haveEqtlInfo);
        } else {
        }
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

double Omics::tuneEigenCutoff(Data &data, const Options &opt){
    cout << "\nFinding the best eigen cutoff from [" << opt.eigenCutoff.transpose() << "] based on pseudo summary data validation." << endl;
    
    Gadget::Timer timer;
    timer.setTime();
    
    unsigned numKeptInds = data.numKeptInds;    
    data.numKeptInds = data.pseudoGwasNtrn;
    
    unsigned size = opt.eigenCutoff.size();
    VectorXd cor(size);
    VectorXd rel(size);
    
    cout << boost::format("%10s %25s %20s\n") % "Cutoff" % "Prediction accuracy (r)" % "Relative accuracy";
    
    for (unsigned i=0; i<size; ++i) {
        double cutoff = opt.eigenCutoff[i];
        cout << boost::format("%10s") % cutoff;

        data.readEigenMatrixBinaryFileAndMakeWandQ(opt.eigenMatrixFile, cutoff, data.pseudoGwasEffectTrn, data.pseudoGwasNtrn, false, false);
        // data.readEigenMatrixBinaryFile(opt.eigenMatrixFile, cutoff);
        // data.constructWandQ(data.pseudoGwasEffectTrn, data.pseudoGwasNtrn);

        data.initVariances(opt.heritability, opt.propVarRandom);
        bool print = false;
        Model *modeli = new ApproxBayesR(data, data.lowRankModel, data.varGenotypic, data.varResidual, opt.pis, opt.piPar, opt.gamma, opt.estimatePi, opt.estimateSigmaSq, opt.noscale, opt.originalModel, opt.overdispersion, opt.estimatePS, opt.spouseCorrelation, opt.diagnosticMode, opt.robustMode, opt.algorithm, print);
        
        vector<McmcSamples*> mcmcSampleVeci;
        MCMC mcmc;

        unsigned chainLength = 150;
        unsigned burnin = 100;
        unsigned thin = 1;
        mcmcSampleVeci = mcmc.run(*modeli, chainLength, burnin, thin, print, opt.outputFreq, opt.title, print, print);

        VectorXd betaMean;
        
        for (unsigned i=0; i<mcmcSampleVeci.size(); ++i) {
            McmcSamples *mcmcSamples = mcmcSampleVeci[i];
            if (mcmcSamples->label == "SnpEffects") {
                betaMean = mcmcSamples->posteriorMean;
            }
        }
        
        // compute prediction accuracy
        cor[i] = betaMean.dot(data.b_val) / sqrt(betaMean.squaredNorm() * data.varPhenotypic);
        rel[i] = cor[i]/cor[0];
        
        cout << boost::format("%25s %20s\n") % cor[i] % rel[i];

    }
    
    data.numKeptInds = numKeptInds;
    
    int bestCutoff_index;
    cor.maxCoeff(&bestCutoff_index);
    double bestCutoff = opt.eigenCutoff[bestCutoff_index];

    if (cor[0] < 0) {
        if (rel.maxCoeff() > -1.25)
            bestCutoff = opt.eigenCutoff[0];
    } else {
        if (rel.maxCoeff() < 1.25)
            bestCutoff = opt.eigenCutoff[0];
    }
    
    timer.getTime();

    if (bestCutoff == opt.eigenCutoff.minCoeff()) {
        cout << "==============================================" << endl;
        cout << "Warning: the best eigen cutoff is the minimum value in the tuning set. We suggest expand the tuning set by including lower candidate values, e.g. --ldm-eigen-cutoff 0.995,0.9,0.8,0.7,0.6  (time used: " << timer.format(timer.getElapse()) << ")." << endl;
        cout << "==============================================" << endl;
    } else {
        cout << bestCutoff << " is selected to be the eigen cutoff to continue the analysis (time used: " << timer.format(timer.getElapse()) << ")."  << endl;
    }
    
    return bestCutoff;
}

double Omics::tuneGeneEigenCutoff(Data &data, const Options &opt){
    cout << "\nFinding the best eigen cutoff from [" << opt.eigenCutoff.transpose() << "] based on pseudo summary data validation." << endl;
    
    Gadget::Timer timer;
    timer.setTime();
    
    unsigned numKeptInds = data.numKeptInds;    
    data.numKeptInds = data.pseudoGwasNtrn;
    
    unsigned size = opt.eigenCutoff.size();
    VectorXd cor(size);
    VectorXd rel(size);
    
    cout << boost::format("%10s %25s %20s\n") % "Cutoff" % "Prediction accuracy (r)" % "Relative accuracy";
    
    for (unsigned i=0; i<size; ++i) {
        double cutoff = opt.geneEigenCutoff[i];
        cout << boost::format("%10s") % cutoff;

        data.UseGeneEigenMakeWAndQgene(cutoff, data.pseudoGwasEffectTrn, data.pseudoGwasNtrn, false, false);
        //data.readEigenMatrixBinaryFile(opt.eigenMatrixFile, cutoff);
        //data.constructWandQ(data.pseudoGwasEffectTrn, data.pseudoGwasNtrn);

        data.initeQTLVariances(opt.cisHeritability,opt.propVarRandom);
        bool print = false;
        // Model *modeli = new ApproxBayesR(data, data.lowRankModel, data.varGenotypic, data.varResidual, opt.pis, opt.piPar, opt.gamma, opt.estimatePi, opt.estimateSigmaSq, opt.noscale, opt.originalModel, opt.overdispersion, opt.estimatePS, opt.spouseCorrelation, opt.diagnosticMode, opt.robustMode, opt.algorithm, print);
        
        vector<McmcSamples*> mcmcSampleVeci;
        MCMC mcmc;

        unsigned chainLength = 150;
        unsigned burnin = 100;
        unsigned thin = 1;
        // mcmcSampleVeci = mcmc.run(*modeli, chainLength, burnin, thin, print, opt.outputFreq, opt.title, print, print);

        VectorXd betaMean;
        
        // for (unsigned i=0; i<mcmcSampleVeci.size(); ++i) {
        //     McmcSamples *mcmcSamples = mcmcSampleVeci[i];
        //     if (mcmcSamples->label == "SnpEffects") {
        //         betaMean = mcmcSamples->posteriorMean;
        //     }
        // }
        
        // compute prediction accuracy
        cor[i] = betaMean.dot(data.b_val) / sqrt(betaMean.squaredNorm() * data.varPhenotypic);
        rel[i] = cor[i]/cor[0];
        
        cout << boost::format("%25s %20s\n") % cor[i] % rel[i];

    }
    
    data.numKeptInds = numKeptInds;
    
    int bestCutoff_index;
    cor.maxCoeff(&bestCutoff_index);
    double bestCutoff = opt.eigenCutoff[bestCutoff_index];

    if (cor[0] < 0) {
        if (rel.maxCoeff() > -1.25)
            bestCutoff = opt.eigenCutoff[0];
    } else {
        if (rel.maxCoeff() < 1.25)
            bestCutoff = opt.eigenCutoff[0];
    }
    
    timer.getTime();

    if (bestCutoff == opt.eigenCutoff.minCoeff()) {
        cout << "==============================================" << endl;
        cout << "Warning: the best eigen cutoff is the minimum value in the tuning set. We suggest expand the tuning set by including lower candidate values, e.g. --ldm-eigen-cutoff 0.995,0.9,0.8,0.7,0.6  (time used: " << timer.format(timer.getElapse()) << ")." << endl;
        cout << "==============================================" << endl;
    } else {
        cout << bestCutoff << " is selected to be the eigen cutoff to continue the analysis (time used: " << timer.format(timer.getElapse()) << ")."  << endl;
    }
    
    return bestCutoff;
}

Model* Omics::buildModel(Data &data, const string &bedFile, const string &gwasFile, const bool &haveXqtlDataBool, const string &bayesType,const string mcmcType, 
                const bool eieoLatent, const bool sampleVareBool, const bool sampleVarEpsBool,
                const unsigned windowWidth, const double heritability, const double cisHeritability, const double propVarRandom, 
                const double pi, const double piEffEqtl,const double piGenicGwas,const double piGenicEqtl,
                const double piEffNonEqtl,const VectorXd &piEffEqtlVec,
                const VectorXd &piEffNonEqtlVec, const double piTheta, const double piAlpha,const double piBeta, 
                const bool estimatePi, const bool noscale,const VectorXd &pis, const VectorXd &piPar, 
                const VectorXd &gamma,const bool estimateSigmaSq,const double phi, const string &algorithm, 
                const double overdispersion,const bool estimatePS,const double icrsq, const double spouseCorrelation, 
                const bool diagnosticMode, 
                const bool originalModel, const bool perSnpGV, const bool robustMode){
    data.initVariances(heritability, propVarRandom);
    if (bayesType == "CO" || bayesType == "RO"){
        if(haveXqtlDataBool){
        // if (!geneInfoFile.empty() || !eqtlSummaryFile.empty()){
        // if ( !eqtlSummaryFile.empty() || !geneInfoFile.empty()){
            // need to initial eQTL variances. 
            data.initeQTLVariances(cisHeritability,propVarRandom);
        }
    }
    // Summary level algorithms
    if (!gwasFile.empty()) {
      if (bayesType == "C"){
        return new ApproxBayesC(data, data.lowRankModel, data.varGenotypic, data.varResidual, data.varRandom,
                                pi, piAlpha, piBeta, estimatePi, noscale, phi, overdispersion, estimatePS, icrsq,
                                spouseCorrelation, diagnosticMode, robustMode);
      } else if(bayesType == "CO"){
        return new ApproxBayesCO(data,mcmcType,eieoLatent,sampleVareBool,sampleVarEpsBool, data.varGenotypic,
                                 data.varResidual,data.varRandom,heritability,cisHeritability ,pi, piEffEqtl,
                                 piGenicGwas, piGenicEqtl,
                                 piTheta,piEffNonEqtl, piAlpha, piBeta, estimatePi, noscale, phi, overdispersion, 
                                 estimatePS, icrsq, spouseCorrelation, diagnosticMode, robustMode);
      } else if (bayesType == "R"){
        return new ApproxBayesR(data, data.lowRankModel, data.varGenotypic, data.varResidual, pis, piPar, gamma,
                                estimatePi, estimateSigmaSq, noscale, originalModel, overdispersion, estimatePS, 
                                spouseCorrelation, diagnosticMode, robustMode, algorithm);
      } else if(bayesType == "RO"){
        // return new 
        return new ApproxBayesRO(data,mcmcType,eieoLatent, sampleVareBool,sampleVarEpsBool, data.varGenotypic, 
                                data.varResidual,data.varRandom,heritability,cisHeritability ,pi,piEffEqtl,
                                piGenicGwas, piGenicEqtl,piEffNonEqtl,
                                piEffEqtlVec,piEffNonEqtlVec,piTheta, piAlpha, piBeta, estimatePi, noscale, phi,
                                overdispersion, estimatePS, icrsq, spouseCorrelation,diagnosticMode, robustMode,pis, 
                                piPar, gamma, originalModel,algorithm);
      } else{
        LOGGER.e(0," Wrong bayes type: " + bayesType + " in the summary-data-based Bayesian analysis.");
      }
    } // end of summary level methods
    // Individual level algorithms
    if (bayesType == "C") {
        data.readBedFile(noscale, bedFile + ".bed");
        return new BayesC(data, data.varGenotypic, data.varResidual, data.varRandom, pi, piAlpha, piBeta, estimatePi, noscale, algorithm);
    } else if (bayesType == "R") {
        data.readBedFile(noscale, bedFile + ".bed");
        return new BayesR(data, data.varGenotypic, data.varResidual, data.varRandom, pis, piPar, gamma, estimatePi, noscale, originalModel, algorithm);
    } else if (bayesType == "CO") {
        data.readBedFile(noscale, bedFile + ".bed");
        // if ( !geneInfoFile.empty()){
        if(haveXqtlDataBool){
            data.buildBayesOmicsMME(bedFile + ".bed",noscale,true);
        } else {
            data.buildBayesOmicsMME(bedFile + ".bed", noscale,false);
        }
        return new BayesCO(data,mcmcType,eieoLatent, data.varGenotypic, data.varResidual,data.varRandom,heritability,cisHeritability,
                           pi,piEffEqtl,piTheta,piEffNonEqtl, piAlpha, piBeta, estimatePi, noscale, phi, overdispersion, estimatePS, 
                           icrsq, spouseCorrelation, diagnosticMode, robustMode);
    } else if (bayesType == "RO") {
        data.readBedFile(noscale, bedFile + ".bed");
        // if ( !geneInfoFile.empty()){
        if(haveXqtlDataBool){
            data.buildBayesOmicsMME(bedFile + ".bed",noscale,true);
        } else {
            data.buildBayesOmicsMME(bedFile + ".bed", noscale,false);
        }
        return new BayesRO(data,mcmcType,eieoLatent, data.varGenotypic, data.varResidual,data.varRandom,heritability,cisHeritability,
                           pi,piEffEqtl,piEffNonEqtl,piEffEqtlVec,piEffNonEqtlVec,piTheta,piAlpha, piBeta,estimatePi, noscale, phi, 
                           overdispersion, estimatePS, icrsq, spouseCorrelation, diagnosticMode, robustMode,pis, piPar, gamma, 
                           originalModel,algorithm);
    } else {
        LOGGER.e(0," Wrong bayes type: " + bayesType);
    }

    return 0;
}

vector<McmcSamples*> Omics::runMcmc(Model &model, const unsigned chainLength, const unsigned burnin, const unsigned thin, const unsigned outputFreq, const string &title, const bool writeBinPosterior, const bool writeTxtPosterior){
    MCMC mcmc;
    return mcmc.run(model, chainLength, burnin, thin, true, outputFreq, title, writeBinPosterior, writeTxtPosterior);
}

void Omics::saveMcmcSamples(const vector<McmcSamples*> &mcmcSampleVec, const string &filename){
    for (unsigned i=0; i<mcmcSampleVec.size(); ++i) {
        McmcSamples *mcmcSamples = mcmcSampleVec[i];
        if (mcmcSamples->label == "SnpEffects" )  continue;
        if (mcmcSamples->label == "WindowDelta") continue;
        // mcmcSamples->writeDataTxt(filename);
    }
}
void Omics::outputResults(const Data &data, const vector<McmcSamples*> &mcmcSampleVec, const string &bayesType,const string &mcmcType, const bool noscale, const string &filename){
    vector<McmcSamples*> mcmcSamplesPar;
    bool keepCisHsq = false;
    bool keepGwasCisHsq = false;
    for (unsigned i=0; i<mcmcSampleVec.size(); ++i) {
        McmcSamples *mcmcSamples = mcmcSampleVec[i];
        if (mcmcSamples->label == "SnpEffects") {
            McmcSamples *pip = NULL;
            for (unsigned i=0; i<mcmcSampleVec.size(); ++i) {
                if (mcmcSampleVec[i]->label == "PIP") {
                    pip = mcmcSampleVec[i];
                    break;
                }
            }
            if (pip == NULL) {
                data.outputSnpResults(mcmcSamples->posteriorMean, mcmcSamples->posteriorSqrMean, mcmcSamples->lastSample, mcmcSamples->pip, noscale, filename + ".snpRes");
            } else
                data.outputSnpResults(mcmcSamples->posteriorMean, mcmcSamples->posteriorSqrMean, mcmcSamples->lastSample, pip->posteriorMean, noscale, filename + ".snpRes");
        }
        else if (mcmcSamples->label == "CovEffects") {
            if (mcmcSamples->datMat.size()) data.outputFixedEffects(mcmcSamples->datMat, filename + ".covRes");
        }
        else if (mcmcSamples->label == "RandCovEffects") {
            data.outputRandomEffects(mcmcSamples->datMat, filename + ".randCovRes");
        }
        else if (mcmcSamples->label == "WindowDelta") {
            data.outputWindowResults(mcmcSamples->posteriorMean, filename + ".window");
        }
        else {
            mcmcSamplesPar.push_back(mcmcSamples);
        }
    }
    if(bayesType == "CO" || bayesType == "RO"){
        outputGeneParameter(data,mcmcSampleVec,mcmcType,filename + ".genePar"); 
    }
    outputPartitionedSnpResults(data,mcmcSampleVec,mcmcType,noscale,filename + ".gene.snpRes.gz");

}


void Omics::outputGeneParameter(const Data &data,const vector<McmcSamples*> &mcmcSampleVec,const string &mcmcType,const string &filename) {
    /// select dataset
    // vector<McmcSamples*> mcmcSamplesPar;
    VectorXd cisMean, cisSqrMean;
    VectorXd genicCisEnrichMean, genicCisEnrichSqrMean;
    VectorXd gwasCisMean,gwasCisSqrMean;
    VectorXd genicGwasEnrichMean, genicGwasEnrichSqrMean;
    VectorXd geneEffMean, geneEffSqrMean;

    for (unsigned i=0; i<mcmcSampleVec.size(); ++i) {
        McmcSamples *mcmcSamples = mcmcSampleVec[i];
        if(mcmcSamples->label == "GeneEffects"){
            // data.outputGeneEffectResults(mcmcSamples->posteriorMean, mcmcSamples->posteriorSqrMean, mcmcSamples->lastSample, pip->posteriorMean,mcmcType,filename + ".geneRes");
            geneEffMean = mcmcSamples->posteriorMean;
            geneEffSqrMean = mcmcSamples->posteriorSqrMean;
        }
        else if (mcmcSamples->label == "cisHsq"){ 
            cisMean = mcmcSamples->posteriorMean;
            cisSqrMean = mcmcSamples->posteriorSqrMean;
        }
        else if (mcmcSamples->label == "gwasCisHsq"){
            gwasCisMean = mcmcSamples->posteriorMean;
            gwasCisSqrMean = mcmcSamples->posteriorSqrMean;
        }  
        else if (mcmcSamples->label == "genicGwasEnrich"){
            genicGwasEnrichMean = mcmcSamples->posteriorMean;
            genicGwasEnrichSqrMean = mcmcSamples->posteriorSqrMean;
        } 
        else if (mcmcSamples->label == "genicEqtlEnrich"){
            genicCisEnrichMean = mcmcSamples->posteriorMean;
            genicCisEnrichSqrMean = mcmcSamples->posteriorSqrMean;
        } 
        else {
            // mcmcSamplesPar.push_back(mcmcSamples);
        }
    }
    

    /// output dataset
    ofstream out(filename.c_str());
    out << boost::format("%6s %20s %6s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s")
    % "Id"
    % "GeneID"
    % "Chrom"
    % "GeneEff"
    % "GeneEffSE"
    % "CTHsq"
    % "CTHsqSE"
    % "CTHsqEn"
    % "CTHsqEnSE"
    % "CGHsq"
    % "CGHsqSE"
    % "CGHsqEn"
    % "CGHsqEnSE";
    out << endl;
    for (unsigned i=0, idx=0; i< data.numKeptGenes; ++i) {
        GeneInfo *gene = data.keptGeneInfoVec[i];
        out << boost::format("%6s %20s %6s %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e")
        % (idx+1)
        % gene->ensemblID
        % gene->chrom
        // % gene->midPhyPos
        % geneEffMean[i]
        % sqrt(geneEffSqrMean[i] - geneEffMean[i]*geneEffMean[i])
        % gwasCisMean[i]
        % sqrt(gwasCisSqrMean[i] - gwasCisMean[i]*gwasCisMean[i])
        % genicGwasEnrichMean[i]
        % sqrt(genicGwasEnrichSqrMean[i] - genicGwasEnrichMean[i]*genicGwasEnrichMean[i])
        % cisMean[i]
        % sqrt(cisSqrMean[i] - cisMean[i]*cisMean[i])
        % genicCisEnrichMean[i]
        % sqrt(genicCisEnrichSqrMean[i] - genicCisEnrichMean[i]*genicCisEnrichMean[i]);
        out << endl;
        ++idx;
    }
    out.close();
}


void Omics::outputPartitionedSnpResults(const Data &data,const vector<McmcSamples*> &mcmcSampleVec,const string &mcmcType,const bool noscale, const string &filename){
    // output gene-snp pair
    map<string,McmcSamples*> eQTLJointVecMap,snpEffectVecMap,eQTLJointVecPipMap,snpEffectVecPipMap;
    boost::regex patternEqtl("EQTLJointVec_(.*)");
    boost::regex patternSNP("SnpJointVec_(.*)");
    boost::regex patternEqtlPip("deltaEQTL_(.*)");
    boost::regex patternGWASPip("deltaGWAS_(.*)");
    smatch matches;
    string geneID;
    for(unsigned i = 0; i < mcmcSampleVec.size(); i++){
        // snp
        if (regex_match(mcmcSampleVec[i]->label, matches, patternSNP)) {
            if (matches.size() > 1 && mcmcSampleVec[i]->label != "SnpJointVec_nonEqtl" ) snpEffectVecMap.insert(pair<string,McmcSamples*>(matches[1].str(),mcmcSampleVec[i]));
        } 
        // pip for snp 
        if (regex_match(mcmcSampleVec[i]->label, matches, patternGWASPip)) {
            if (matches.size() > 1  && mcmcSampleVec[i]->label != "deltaGWAS_nonEqtl") snpEffectVecPipMap.insert(pair<string,McmcSamples*>(matches[1].str(),mcmcSampleVec[i]));
        } 
        if (regex_match(mcmcSampleVec[i]->label, matches, patternEqtl)) {
            if (matches.size() > 1) eQTLJointVecMap.insert(pair<string,McmcSamples*>(matches[1].str(),mcmcSampleVec[i]));
        } 
        // pip for eqtl
        if(mcmcType == "AIAO"){
            eQTLJointVecPipMap = snpEffectVecPipMap;
        } else {
            if (regex_match(mcmcSampleVec[i]->label, matches, patternEqtlPip)) {
                if (matches.size() > 1) eQTLJointVecPipMap.insert(pair<string,McmcSamples*>(matches[1].str(),mcmcSampleVec[i]));
            } 
        }
    }
    if(snpEffectVecMap.size() ==0 || eQTLJointVecMap.size() == 0){
        LOGGER.w(0,"Partitioned gwas snp and xqtl effect is zero when generating SNP results in the .gene.snpRes.gz file ");
        return;
    }

    // save results
    std::ofstream file(filename, std::ios_base::out | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::output> outGZ;
    outGZ.push(boost::iostreams::gzip_compressor());
    outGZ.push(file);
    std::ostream out(&outGZ);
    // ofstream out(filename.c_str());

    out << boost::format("%-15s %-12s %6s %12s %6s %6s %12s %-15s %-15s %-15s %-15s %-15s %-15s")
    % "GeneID"
    % "SNPID"
    % "Chrom"
    % "Position"
    % "A1"
    % "A2"
    % "A1Frq"
    % "beta"
    % "beta_se"
    % "beta_pip"
    % "alpha"
    % "alpha_se"
    % "alpha_pip";
    out << endl;
    // use geneMat

    for(unsigned j = 0; j < data.numKeptGenes;j++){
        vector<string> snpIDInGene = data.gene2cisSnpIDMap.at(j);
        string geneID = data.geneEffectNames[j];
        for(unsigned i = 0; i < snpIDInGene.size();i++){
            EqtlInfo * eqtl = data.eqtlInfoMap.find(snpIDInGene[i])->second;
            SnpInfo * snp = data.snpInfoMap.find(snpIDInGene[i])->second;
            // if(abs(snpEffectVecMap.at(geneID)->posteriorMean[i]) - 1e-6 < 0 &&
            //    abs(eQTLJointVecMap.at(geneID)->posteriorMean[i]) - 1e-6 < 0){
            //     continue; // here we only keep nonzero gene-snp pair
            //    }
            // here we not weight 2pq to avoid allele frequency discrepancy issue
            // double sqrt2pq = sqrt(2.0 * eqtl->af * (1- eqtl->af));
            double sqrtScaleFactoreQTL = (data.scalingeQTLFactorVecVec[j][i]);
            double sqrtScaleFactorGWAS = (snp->scaleFactor);
            //// beta effect
            double betaEffect = (eqtl->flipped ? - snpEffectVecMap.at(geneID)->posteriorMean[i] : snpEffectVecMap.at(geneID)->posteriorMean[i]);
            // double lastBeta = (snp->flipped ? -lastSample[idx] : lastSample[idx]);
            double betaSE = sqrt(snpEffectVecMap.at(geneID)->posteriorSqrMean[i]-snpEffectVecMap.at(geneID)->posteriorMean[i] * snpEffectVecMap.at(geneID)->posteriorMean[i]);
            /// alpha effect 
            double alphaEffect = (eqtl->flipped ? -eQTLJointVecMap.at(geneID)->posteriorMean[i] : eQTLJointVecMap.at(geneID)->posteriorMean[i]);
            // double lastBeta = (snp->flipped ? -lastSample[idx] : lastSample[idx]);
            double alphaSE = sqrt(eQTLJointVecMap.at(geneID)->posteriorSqrMean[i]-eQTLJointVecMap.at(geneID)->posteriorMean[i] * eQTLJointVecMap.at(geneID)->posteriorMean[i]);

            // here we need to outline 
            out << boost::format("%-15s %-12s %6s %12s %6s %6s %12.6f %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e %-15.6e")
            % geneID
            % eqtl->rsID
            % eqtl->chrom
            % eqtl->physPos
            % (eqtl->flipped ? eqtl->a2 : eqtl->a1)
            % (eqtl->flipped ? eqtl->a1 : eqtl->a2)
            % (eqtl->flipped ? 1.0-eqtl->af : eqtl->af)
            % (noscale ? betaEffect : betaEffect/sqrtScaleFactorGWAS)
            % (noscale ? betaSE : betaSE/sqrtScaleFactorGWAS)
            // % snpEffectVecPipMap.at(geneID)->posteriorMean[i]
            % -9
            % (noscale ? alphaEffect : alphaEffect/sqrtScaleFactoreQTL)
            % (noscale ? alphaSE : alphaSE/sqrtScaleFactoreQTL)
            // % eQTLJointVecPipMap.at(geneID)->posteriorMean[i];
            % -9;
            out << endl;

            // cout << " pip " << snpEffectVecPipMap.at(geneID)->pip << " ";
        }
    }
    boost::iostreams::close(outGZ);
    // out.close();
}

McmcSamples* Omics::inputMcmcSamples(const string &mcmcSampleFile, const string &label, const string &fileformat){
    LOGGER << "reading MCMC samples for " << label << endl;
    McmcSamples *mcmcSamples = new McmcSamples(label);
    if (fileformat == "bin") mcmcSamples->readDataBin(mcmcSampleFile + "." + label);
    if (fileformat == "txt") mcmcSamples->readDataTxt(mcmcSampleFile + ".Par", label);
    return mcmcSamples;
}

void Omics::clearGenotypes(Data &data){
    data.X.resize(0,0);
}

void Omics::pip2p(const Data &data, const VectorXd &pip, const double propNull, VectorXd &pval){
    VectorXd pipSrt = pip;
    std::sort(pipSrt.data(), pipSrt.data() + pipSrt.size(), greater <double> ());
    double cumsum = 0.0;
    double pvali = 0.0;
    double numNull = data.numIncdSnps*propNull;
    map<double, double> pip2pMap;
    for (unsigned i=0; i<data.numIncdSnps; ++i){
        cumsum += pipSrt[i];
        pvali = (i+1 - cumsum) / numNull;
        pip2pMap[pipSrt[i]] = pvali;
    }
    for (unsigned i=0; i<data.numIncdSnps; ++i){
        pval[i] = pip2pMap[pip[i]];
    }
}


