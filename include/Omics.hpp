#ifndef amber_hpp
#define amber_hpp

#include "Data.hpp"
#include "Options.hpp"
#include "Logger.hpp"
#include "Gadgets.hpp"


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
    
    // void inputSnpInfo(Data &data, const string &bedFile, const string &gwasSummaryFile, const double afDiff, const double mafmin, const double mafmax, const double pValueThreshold, const bool sampleOverlap, const bool imputeN, const bool noscale);
   
    void inputSnpInfo(Data &data, const string &includeSnpFile, const string &excludeSnpFile, const string &excludeRegionFile,
                        const string &gwasSummaryFile, const string &ldmatrixFile, const unsigned includeChr, const bool excludeAmbiguousSNP,
                        const string &skeletonSnpFile, const string &geneticMapFile, const double genMapN, const unsigned flank, const string &eQTLFile, const string &ldscoreFile, const string &windowFile,
                        const bool multiLDmatrix, const bool excludeMHC, const double afDiff, const double mafmin, const double mafmax, const double pValueThreshold, const double rsqThreshold, const bool sampleOverlap, const bool imputeN, const bool noscale, const bool binSnp, const bool readLDMfromTxtFile);

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

    void inputSnpInfo(Data &data, const string &bedFile, const string &gwasSummaryFile, const double afDiff, const double mafmin, const double mafmax, const double pValueThreshold, const bool sampleOverlap, const bool imputeN, const bool noscale);

    void clearGenotypes(Data &data);

};

#endif /* amber_hpp */
