// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BayesOmics, a statistical genetics software package
// developed by Shouye Liu.
//
// Copyright (C) 2025 Shouye Liu <syliu.xue@foxmail.com>
//
// This file is licensed under the GNU General Public License v3.0 or (at your option)
// any later version. You may redistribute and/or modify it under the terms of the GPL.
//
// BayesOmics is distributed WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the LICENSE file for details.


#ifndef PROGRAMOPTIONS_HPP
#define PROGRAMOPTIONS_HPP

#include <vector>
#include <string>
#include <utility>
#include<iostream>
#include<Logger.hpp>
#include <omp.h>

#include <Eigen/Core>
#include <Eigen/Eigen>
#include <Eigen/Dense>

#include <boost/program_options.hpp>
#include <boost/format.hpp>

//#include "DataMatrix.hpp"
#include "Gadgets.hpp"

using namespace std;
using namespace boost;
using namespace Eigen;
const double Megabase = 1000000.0;
class Options
{
    public:
    virtual ~Options()= default;;
    void readProgramOptions(int argc, const char * const argv[]);
    void printCommandLine(int argc, const char * const argv[]);

    private:
    bool boostProgramOptionsRoutine(int argc, const char* const argv[]);
    void customCheck();
    void setThread(void);

    public:
    // Here we define various option variables
    boost::program_options::variables_map vm;
    // unsigned seed;
    //std::string famFile;
    static const double MIX_PARAM_ESTIMATE_FLAG; // flag for estimating mixture params using CV

    // debug mode
    bool diagnosticMode; // for sbayes debug
    bool dataMode; 
    // Basic settings
    std::string title;
    std::string analysisType;
    int numThreads;
    int seed;
    string algorithm;
    string optionFile;
    int slurmArrayLimit;

    /// Data management 
    string rsidMapFiles;
    string proteinMapFile;
    string proteinPath;
    string grchType;
    bool makeAnnoBool;
    bool isAnnoBinary;
    bool makeEigenGeneBasedOnGeneSnpPair;
    // besd format
    bool haveXqtlDataBool;
    bool makeBesdBool;
    bool makeBesdSmrBool;
    bool makeQueryBool;
    bool makeSumstatsFormatBool;
    bool makeQueryCojoMaBool;
    bool makeXqtlNBool;
    int reEffectSampleSize;
    int reSamType;
    bool useCisWinBool;
    bool mergeBesdBool;
    bool mergeEigenGeneBool;

    /// Genotype-related settings
    unsigned windowWidth; // in mega-base unit
    double cisRegionWind; // used to define cis-region window 
    unsigned includeChr;  // chromosome to include
    unsigned flank;
    unsigned includeBlock;  // block to include
    string includeBlockID;
    bool multiLDmat;
    bool writeLdmTxt;      // write ldm to txt file
    bool readLdmTxt;      // read ldm from a txt file
    bool excludeMHC;  // exclude SNPs in the MHC region
    bool directPrune; // direct prune ldm
    bool jackknife;   // jackknife estimate for LD sampling variance
    bool excludeAmbiguousSNP;  // exlcude ambiguous SNPs with A/T or G/C alleles
    bool sampleOverlap;  // whether LD ref is the same as GWAS sample
    bool imputeN;  // impute per-SNP sample size
    string bedFile;
    string alleleFreqFile;
    string includeSnpFile;
    string excludeSnpFile;
    string excludeRegionFile;
    string snpResFile;
    string gwasSummaryFile;
    string ldmatrixFile;
    string skeletonSnpFile;
    string ldscoreFile;
    string snpRange;
    string partParam;
    string outLDmatType;
    string windowFile;
    string residualDiagFile;
    string eigenMatrixFile;
    string geneEigenMatrixFile;
    string ldBlockInfoFile;
    string geneListFile; // read gene infomation for eigen-decomposition in SBayesOmics
    bool outInfoOnly;
    string geneInfoFile; // read plist file in BayesOmics 
    string keepIndGeneFile;

    /// Phenotype-related settings
    unsigned mphen; // triat order id in the phenotype file for analysis
    unsigned keepIndMax;  // the maximum number of individuals kept for analysis
    string phenotypeFile;
    string covariateFile;
    string randomCovariateFile;
    string keepIndFile;

    /// gene-related settings
    string subGenePath;
    string eQTLFile;
    string bedGeneFile;
    string eqtlSummaryFile;
    string eqtlSummaryQueryFile;
    string geneticMapFile;
    string geneSamSizeFile;
    VectorXd eigenCutoff;
    VectorXd geneEigenCutoff;
    string includeGeneFile;
    string excludeGeneFile;
    string includeSpecificGeneID; // 
    /// MCMC settings
    string mcmcType;
    
    bool eieoLatentBool;
    bool sampleVareBool;
    bool sampleVarEpsBool;
    unsigned numChains;
    int chainLength;
    unsigned burnin;
    unsigned outputFreq;
    unsigned thin;  // save every this th sampled value in MCMC
    bool writeBinPosterior;
    bool writeTxtPosterior;
    bool outputResults;
    string mcmcSampleFile;

    /// Bayesian-model settings
    string bayesType;
    double pi;
    double piAlpha;
    double piBeta;
    double piEffEqtl;
    double piTheta;
    double piEffNonEqtl;
    double piGenicGwas; // piEffEieo1
    double piGenicEqtl; // piEffEieo2;
    double heritability;
    double cisHeritability;
    double propVarRandom;  // proportion of variance explained by random covariate effects
    double LDthreshold;  // used to define the two ends of per-SNP LD window in the banded LD matrix
    
    double effpopNE;  // for shrunk LDM
    double genMapN;   // for shrunk LDM
    double cutOff;    // for shrunk LDM
    double chisqThreshold;  // significance threshold for nonzero LD chi-square test
    double phi;   // a shrinkage parameter for the heritability estimate in sbayes
    double overdispersion;
    double icrsq;  // average inter-chromosome r^2 across SNPs
    double spouseCorrelation;
    double afDiff; // filtering SNPs by the allele frequency difference in LD and GWAS samples
    double mafmin;  // lower bound of maf
    double mafmax;  // upper bound of maf
    double lambda;  // for conjugate gradient
    double rsqThreshold;
    double pValueThreshold;
    bool estimatePS;  // estimate population stratification in sbayes
    // Bayes R defauls
    unsigned ndists; // Number of distributions for base Bayes R
    VectorXd gamma;  // Default scaling parameters for Bayes R
    VectorXd pis;    // Default pis for Bayes R

    // BayesOmics R defaults
    VectorXd piEffEqtlVec;    // Default pis for Bayes R
    VectorXd piEffNonEqtlVec;    // Default pis for Bayes R
    
    // hyperparameters for the prior distributions
    VectorXd piPar;
    Vector2f piNDCpar;

    bool estimatePi;
    bool estimateSigmaSq; // variance of SNP effects
    bool estimateScale;
    bool noscale;
    bool simuMode; // simulation mode
    bool originalModel; // original BayesR model
    bool twoStageModel;  // two-step approach for estimating X-chr dosage model and G by sex
    bool binSnp;  // bin SNPs
    bool robustMode;  // use the robust parameterisation in SBayes models
    bool perSnpGV;
    bool mergeLdm;
    bool imputeSummary;
    

    Options(){
        // debug
        diagnosticMode = false; // for sbayes
        dataMode = false;
        // Basic settings
        title = "BayesOmics";
        analysisType = "bayes";
        numThreads = 1;
        seed = 2023;
        algorithm = "";
        optionFile = "";
        slurmArrayLimit = 1000;

        /// Data management 
        rsidMapFiles = "";
        proteinMapFile = "";
        proteinPath = "";
        grchType = "hg38";
        makeAnnoBool = false;
        isAnnoBinary = false;
        makeEigenGeneBasedOnGeneSnpPair = false;
        makeBesdBool = false;
        haveXqtlDataBool = false;
        makeBesdSmrBool = false;
        makeQueryBool = false;
        makeSumstatsFormatBool = false;
        makeQueryCojoMaBool = false;
        makeXqtlNBool = false;
        reEffectSampleSize = 0;
        reSamType = 0; // 0: resample by imputation; 1: resample by cross-trait ldsc
        useCisWinBool = false;
        mergeBesdBool = false;
        mergeEigenGeneBool = false;
        


        /// Genotype-related settings
        windowWidth = 0 * Megabase; // in mega-base unit
        cisRegionWind = 0.1 * Megabase; // used to define cis-region window
        includeChr = 0;  // chromosome to include
        flank = 0;
        includeBlock = 0;  // block to include
        includeBlockID= "";
        multiLDmat = false;
        writeLdmTxt = false;          // write ldm to txt file
        readLdmTxt = false;          // read ldm from a txt file
        excludeMHC = false;  // exclude SNPs in the MHC region
        jackknife = false;   // jackknife estimate for LD sampling variance
        directPrune = false;
        excludeAmbiguousSNP = false;  // exlcude ambiguous SNPs with A/T or G/C alleles
        sampleOverlap = false;  // whether LD ref is the same as GWAS sample
        imputeN = false;  // impute per-SNP sample size
        bedFile = "";
        alleleFreqFile = "";
        includeSnpFile = "";
        excludeSnpFile = "";
        excludeRegionFile = "";
        snpResFile = "";
        gwasSummaryFile = "";
        ldmatrixFile = "";
        skeletonSnpFile = "";
        ldscoreFile = "";
        snpRange = "";
        partParam = "";
        outLDmatType = "";
        windowFile = "";
        residualDiagFile = "";
        eigenMatrixFile = "";
        geneEigenMatrixFile = "";
        ldBlockInfoFile = "";
        geneListFile = "";
        outInfoOnly = false;
        geneInfoFile = ""; // BayesOmics 
        keepIndGeneFile = "";

        /// Phenotype-related settings
        mphen = 1; // triat order id in the phenotype file for analysis
        keepIndMax = UINT_MAX;  // the maximum number of individuals kept for analysis
        phenotypeFile = "";
        covariateFile = "";
        randomCovariateFile = "";
        keepIndFile = "";

        /// gene-related settings
        subGenePath = "";
        eQTLFile = "";
        bedGeneFile = "";
        eqtlSummaryFile = "";
        eqtlSummaryQueryFile = "";
        geneticMapFile = "";
        geneSamSizeFile = "";
        eigenCutoff.resize(5);
        eigenCutoff << 0.9999, 0.995, 0.99, 0.95, 0.9;
        geneEigenCutoff.resize(5);
        geneEigenCutoff << 0.9999, 0.995, 0.99, 0.95, 0.9;
        includeGeneFile = "";
        excludeGeneFile = "";
        includeSpecificGeneID = "";

        propVarRandom = 0.05;  // proportion of variance explained by random covariate effects
        LDthreshold = 0.0;  // used to define the two ends of per-SNP LD window in the banded LD matrix

        // Shrunk matrix defaults
        effpopNE  = 11490.672741;
        cutOff    = 1e-5;
        genMapN = 183; // Sample size of CEU population
        
        chisqThreshold = 10;  // significance threshold for nonzero LD chi-square test
        phi = 0;   // a shrinkage parameter for the heritability estimate in sbayes
        overdispersion = 0;
        icrsq = 0;  // average inter-chromosome r^2 across SNPs
        spouseCorrelation = 0;
        afDiff = 0.25; // filtering SNPs by the allele frequency difference in LD and GWAS samples
        mafmin = 0;  // lower bound of maf
        mafmax = 0;  // upper bound of maf
        lambda = 1e6;  // for conjugate gradient
        rsqThreshold = 1.0;
        pValueThreshold = 1.0;
        estimatePS = false;  // estimate population stratification in sbayes

  
        estimatePi = true;
        estimateSigmaSq = true; // variance of SNP effects
        estimateScale = false;
        noscale = false;
        simuMode = false; // simulation mode
        originalModel = false; // original BayesR model
        twoStageModel = false;  // two-step approach for estimating X-chr dosage model and G by sex
        binSnp = false;  // bin SNPs
        robustMode = false;  // use the robust parameterisation in SBayes models
        perSnpGV = false;
        mergeLdm = false;
        imputeSummary = false;
    }
    }; // end of class Options


#endif
