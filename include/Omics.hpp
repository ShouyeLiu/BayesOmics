// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BayesOmics, a statistical genetics software package
// developed by Shouye Liu.
//
// Portions of this file are adapted from GCTB,
// originally licensed under the MIT License. See below for the original license.
//
// Original source: https://cnsgenomics.com/software/gctb/#Download (accessed 20 June 2024)
//
// Modifications and additional code:
// Copyright (C) 2025 Shouye Liu <syliu.xue@foxmail.com>
//
// This file is licensed under the GNU General Public License v3.0 or (at your option)
// any later version. You may redistribute and/or modify it under the terms of the GPL.
//
// BayesOmics is distributed WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the LICENSE file for details.


#ifndef omics_hpp
#define omics_hpp

#include "Data.hpp"
#include "Options.hpp"
#include "Model.hpp"
#include "Mcmc.hpp"
#include "Logger.hpp"
#include "Gadgets.hpp"
#include "ModelBayesCO.hpp"
#include "ModelSBayesCO.hpp"
#include "ModelBayesRO.hpp"
#include "ModelSBayesRO.hpp"


class Omics {
public:
    Options &opt;

    Omics(Options &options): opt(options){};
    
    void inputIndInfo(Data &data, const string &bedFile, const string &phenotypeFile, const string &keepIndFile,
                      const unsigned keepIndMax, const unsigned mphen, const string &covariateFile, 
                      const string &randomCovariateFile, const string &residualDiagFile, const string &geneInfoFile, 
                      const string &keepIndGeneFile,const string subGenePath = "", bool ldBool = false);

    void inputSnpInfo(Data &data, const string &bedFile, const string &geneInfoFile,const double cisRegionWind, const string &includeSnpFile, const string &excludeSnpFile, const string &excludeRegionFile,
                      const unsigned includeChr, const bool excludeAmbiguousSNP, const string &skeletonSnpFile, const string &geneticMapFile, const string &ldBlockInfoFile, const unsigned includeBlock,
                      const unsigned flank, const double mafmin, const double mafmax, const bool noscale, const bool readGenotypes);
   
    void inputSnpInfo(Data &data, const string &includeSnpFile, const string &excludeSnpFile, const string &excludeRegionFile,
                      const string &gwasSummaryFile, const string &ldmatrixFile, const unsigned includeChr, const bool excludeAmbiguousSNP,
                      const string &skeletonSnpFile, const string &geneticMapFile, const double genMapN, const unsigned flank, const string &eQTLFile, const string &ldscoreFile, const string &windowFile,
                      const bool multiLDmatrix, const bool excludeMHC, const double afDiff, const double mafmin, const double mafmax, const double pValueThreshold, const double rsqThreshold, const bool sampleOverlap, const bool imputeN, const bool noscale, const bool binSnp, const bool readLDMfromTxtFile);
    // this function read eigen matrices
    void inputSnpInfo(Data &data, const string &includeSnpFile, const string &excludeSnpFile, const string &excludeRegionFile,
                            const string &gwasSummaryFile,const string &eqtlSummaryFile, const string &eqtlSummaryQueryFile,
                            const string &includeGeneFile,const string &geneSamSizeFile, const string &eigenMatrixFile, 
                            const string &geneEigenMatrixFile, const string &ldBlockInfoFile,
                            const unsigned includeChr, const bool excludeAmbiguousSNP, const unsigned flank, const string &eQTLFile, const string &ldscoreFile,
                            const double eigenCutoff, const double geneEigenCutoff, const bool excludeMHC,
                            const double afDiff, const double mafmin, const double mafmax, const double pValueThreshold, const double rsqThreshold,
                            const bool sampleOverlap, const bool imputeN, const bool noscale, const bool readLDMfromTxtFile, const bool imputeSummary, 
                            const unsigned includeBlock,const string includeBlockID);
    
    void truncBlockEigen(double cutShresh=1e-6);

    void covertUkbPppPGWASToBesdFormat(Data &data, const string rsidMapFiles, const string proteinMapFile, const string proteinPath,
        const bool &makeBesdBool,const bool &makeBesdSmrBool,const bool &makeQueryBool,const bool &makeSumstatsFormatBool,const bool &makeQueryCojoMaBool,
        const double pValueThreshold, const string grchType, const string title);
    void convertFlistFormat2BESDFormat(Data &data,const string title, const string geneInfoFile,const unsigned includeChr,
            const bool makeBesdBool,const bool makeBesdSmrBool, const bool makeQueryBool,const string subGenePath);

    void convertBesdQuery2OtherFormat(Data &data,const string &eigenMatrixFile,
        const string &besdFile,const string &eqtlSummaryQueryFile,const string &geneSamSizeFile, const string &title,const string &includeGeneFile, const string excludeGeneFile,
        const string &specificGeneID, const string &includeSnpFile,const unsigned includeChr,const double &pValueThreshold, 
        const double &cisRegionWind, const bool &useCisWinBool, const int &slurmArrayLimit,const bool &makeBesdBool,
        const bool &makeBesdSmrBool,const bool &makeQueryBool,const bool &makeSumstatsFormatBool,
        const bool &makeQueryCojoMaBool, const bool &makeAnnoBool,const bool &isAnnoBinary, 
        const double &eigenCutoff = 0.9995);

    // resample xQTL sample size
    void resamplexQTLSampleSize(Data &data,const string eigenMatrixFile,const string geneEigenMatrixFile,
    const string besdFile,const string title,const string &includeGeneFile,const string specificGeneID,
    const int reSamType, const bool makeQueryCojoMaBool, 
    const bool makeXqtlNBool, const int reEffectSampleSize, const int slurmArrayLimit, const double eigenCutoff);
    
    // tune for eigen matrix cutoff in LD Block
    double tuneEigenCutoff(Data &data, const Options &opt);
    // tune for eigen matrix cutoff in gene region
    double tuneGeneEigenCutoff(Data &data, const Options &opt);

    void inputSnpInfo(Data &data, const string &bedFile, const string &gwasSummaryFile, const double afDiff, const double mafmin, const double mafmax, const double pValueThreshold, const bool sampleOverlap, const bool imputeN, const bool noscale);

    Model* buildModel(Data &data, const string &bedFile, const string &gwasFile, const bool &haveXqtlDataBool, const string &bayesType,const string mcmcType, 
                const bool eieoLatent, const bool sampleVareBool, const bool sampleVarEpsBool,
                const unsigned windowWidth, const double heritability, const double cisHeritability, const double propVarRandom, 
                const double pi, const double piEffEqtl,const double piGenicGwas,const double piGenicEqtl,
                const double piEffNonEqtl,const VectorXd &piEffEqtlVec,const VectorXd &piEffNonEqtlVec, 
                const double piTheta,const double piAlpha,const double piBeta, const bool estimatePi, const bool noscale,
                const VectorXd &pis, const VectorXd &piPar, const VectorXd &gamma,const bool estimateSigmaSq,const double phi, 
                const string &algorithm, const double overdispersion,
                const bool estimatePS,const double icrsq, const double spouseCorrelation, const bool diagnosticMode, 
                const bool originalModel, const bool perSnpGV, const bool robustMode);

    vector<McmcSamples*> runMcmc(Model &model, const unsigned chainLength, const unsigned burnin, const unsigned thin, const unsigned outputFreq, const string &title, const bool writeBinPosterior, const bool writeTxtPosterior);
    void saveMcmcSamples(const vector<McmcSamples*> &mcmcSampleVec, const string &filename);
    void outputResults(const Data &data, const vector<McmcSamples*> &mcmcSampleVec, const string &bayesType,const string &mcmcType, const bool noscale, const string &filename);
    void outputPartitionedSnpResults(const Data &data,const vector<McmcSamples*> &mcmcSampleVec,const string &mcmcType,const bool noscale, const string &filename);
    void outputGeneParameter(const Data &data,const vector<McmcSamples*> &mcmcSampleVec,const string &mcmcType,const string &filename);
    McmcSamples* inputMcmcSamples(const string &mcmcSampleFile, const string &label, const string &fileformat);
    void clearGenotypes(Data &data);
    void pip2p(const Data &data, const VectorXd &pip, const double propNull, VectorXd &pval);

};

#endif /* amber_hpp */
