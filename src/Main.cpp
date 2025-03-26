#include <iostream>
#include <algorithm>
// #include <cstdlib>
#include <functional>  // std::bind
#include "Logger.hpp"
#include "Options.hpp"
#include "Gadgets.hpp"
#include "Data.hpp"
#include "Omics.hpp"

// #include <boost/format.hpp>
// #include <boost/algorithm/string.hpp>
// #include <boost/iostreams/filtering_stream.hpp>
// #include <boost/iostreams/filter/gzip.hpp>

using std::bind;
using std::map;
using std::to_string;
using std::function;
using std::placeholders::_1;
using std::placeholders::_2;
using std::placeholders::_3;

void logRedirect(bool flag_outFile, string curTime){
    std::function<void (int, const string&, const string&)> log;
    if(flag_outFile){
        log = std::bind(&Logger::l, LOGGER_P, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    }else{
        log = std::bind(&Logger::m, LOGGER_P, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
    }

    log(0, "*******************************************************************", "");
    log(0, "* Bayesian Inference on Omics Data", "");
    log(0, "* Compiled on "+std::string(__DATE__)+" at "+std::string(__TIME__) + " " +TIME_ZONE +" (version " + std::string(BayesOmics_VERSION)+" " + Gadget::getOSName()+")", "");
    // log(0, "* Git Hash: "+ Gadget::getGitCommitHash() ,"");
    log(0, "* Git Hash: "+ std::string(GIT_COMMIT_HASH),"");
    log(0, "* (C) 2023-present, The University of Queensland","");
    log(0, "* Authors: Shouye Liu and Jian Zeng", "");
    log(0, "* Please report bugs to Shouye Liu <shouye.liu@uq.edu.au>", "");
    log(0, "*******************************************************************", "");
    // log(0, "Today quote:\n "+ Gadget::getCurrentQuotes(),"");
    log(0, string("at ") + curTime + "", "Analysis started");
    log(0, "Hostname: " + Gadget::getHostName(), "");
    log(0, "", "");
}

int main(int argc, char *argv[]){
    // precision
    std::cout.precision(10);

    Gadget::Timer timer;
    timer.setTime();
    string curTime = timer.getDate();
    logRedirect(false,curTime);
    LOGGER.ts("main");
    if(argc == 1){
        LOGGER.m(0, "Error: no analysis has been launched by the option(s)");
        LOGGER.m(0, "Run command [BayesOmics64 --help] for specific options.");
        LOGGER.m(0, "and see online documentation(" + string(DOC_ONLINE) + ") for details.");
        return 1;
    }
    Options opt;
    try {
        std::cout.precision(20);
        opt.readProgramOptions(argc,argv);
        Data data;
        data.title = opt.title;
        data.burnIn = opt.burnin;
        LOGGER.open(opt.title + ".log");
        logRedirect(true,curTime);
        opt.printCommandLine(argc, argv); 
        Omics omics(opt);

        if (opt.seed){
            Stat::seedEngine(opt.seed);
            arma::arma_rng::set_seed(opt.seed);
        } 
        else  {
            Stat::seedEngine(011415);  // fix the random seed if not given due to the use of MPI
            arma::arma_rng::set_seed(011415);
        }     

        /***** SET UP SNPDATA *****/
        bool readGenotypes;
        if (opt.analysisType == "DataManagement"){
            // Here we transformat plain file into query or besd format
            if(!opt.geneInfoFile.empty()) omics.convertFlistFormat2BESDFormat(data,opt.title, opt.geneInfoFile,opt.includeChr,opt.makeBesdBool,opt.makeBesdSmrBool, opt.makeQueryBool,opt.subGenePath);
            if(!opt.rsidMapFiles.empty()) omics.covertUkbPppPGWASToBesdFormat(data,opt.rsidMapFiles,opt.proteinMapFile,opt.proteinPath,
            opt.makeBesdBool,opt.makeBesdSmrBool,opt.makeQueryBool,opt.makeSumstatsFormatBool,opt.makeQueryCojoMaBool,
            opt.pValueThreshold, opt.grchType,opt.title);
            // when managing xQTL data, only two formats, including besd and query format, are supported in the final analysis.
            if(!opt.eqtlSummaryFile.empty() || !opt.eqtlSummaryQueryFile.empty()) {
                if((opt.makeAnnoBool || opt.makeBesdBool || opt.makeQueryBool || opt.makeSumstatsFormatBool ||
                    opt.makeQueryCojoMaBool || opt.makeBesdSmrBool) 
                    & !opt.makeXqtlNBool){
                        omics.convertBesdQuery2OtherFormat(data,opt.eigenMatrixFile,opt.eqtlSummaryFile,opt.eqtlSummaryQueryFile, 
                            opt.geneSamSizeFile,
                            opt.title,opt.includeGeneFile,opt.excludeGeneFile,opt.includeSpecificGeneID,opt.includeSnpFile,opt.includeChr,
                            opt.pValueThreshold,opt.cisRegionWind,opt.useCisWinBool, opt.slurmArrayLimit,
                            opt.makeBesdBool,opt.makeBesdSmrBool,opt.makeQueryBool, opt.makeSumstatsFormatBool,
                            opt.makeQueryCojoMaBool,opt.makeAnnoBool,
                            opt.isAnnoBinary,opt.eigenCutoff.maxCoeff());
                    }
                if(opt.mergeBesdBool && !opt.eqtlSummaryFile.empty()) data.mergeMultiBesdData(opt.eqtlSummaryFile,opt.title);
            }
            if(opt.mergeEigenGeneBool) data.mergeMultiEigenMat(opt.title,opt.geneListFile,"gene",opt.eigenCutoff.maxCoeff());
            
        } // end of DataManagement
        else if (opt.analysisType == "LDmatrix") {
            readGenotypes = false;
            if (opt.ldmatrixFile.empty()) { // make LD matrix from genotypes

                /***** Step 1. read individual and genotype infomation *****/
                bool generateLDBool = true;
                omics.inputIndInfo(data, opt.bedFile, "", opt.keepIndFile, opt.keepIndMax,
                                  opt.mphen, opt.covariateFile, opt.randomCovariateFile,
                                opt.residualDiagFile,opt.geneInfoFile,opt.keepIndGeneFile,opt.subGenePath,generateLDBool);
                omics.inputSnpInfo(data, opt.bedFile,opt.geneInfoFile,opt.cisRegionWind, 
                        opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile, opt.includeChr, 
                        opt.excludeAmbiguousSNP, opt.skeletonSnpFile, opt.geneticMapFile, opt.ldBlockInfoFile, 
                        opt.includeBlock, opt.flank, opt.mafmin, opt.mafmax, opt.noscale, readGenotypes);
                
                /***** Step 2. read individual and genotype infomation *****/
                if (opt.outLDmatType == "shrunk") {
                    data.makeshrunkLDmatrix(opt.bedFile + ".bed", opt.outLDmatType, opt.snpRange, opt.title, opt.writeLdmTxt, opt.effpopNE, opt.cutOff, opt.genMapN);
                } else if (opt.outLDmatType == "block") {
                    data.makeBlockLDmatrix(opt.bedFile + ".bed", opt.outLDmatType, opt.includeBlock, opt.title, opt.writeLdmTxt);
                }
                else {
                    string snpRange = opt.snpRange;
                    if(!opt.partParam.empty()){
                        snpRange = data.partLDMatrix(opt.partParam, opt.title, opt.outLDmatType);
                    }
                    data.makeLDmatrix(opt.bedFile + ".bed", opt.outLDmatType, opt.chisqThreshold, opt.LDthreshold, opt.windowWidth, snpRange, opt.title, opt.writeLdmTxt);
                }
            }
            else { // manipulate an existing LD matrix or merge existing LD matrices
                if (opt.mergeLdm) {
                    data.mergeLdmInfo(opt.outLDmatType, opt.ldmatrixFile);
                }
                else if (opt.directPrune) {
                    data.directPruneLDmatrix(opt.ldmatrixFile, opt.outLDmatType, opt.chisqThreshold, opt.title, opt.writeLdmTxt);
                }
                else if (opt.jackknife) {
                    readGenotypes = true;
                    omics.inputIndInfo(data, opt.bedFile, opt.bedFile + ".fam", opt.keepIndFile, opt.keepIndMax, opt.mphen, 
                    opt.covariateFile, opt.randomCovariateFile, opt.residualDiagFile, opt.geneInfoFile,opt.keepIndGeneFile,opt.subGenePath);
                    omics.inputSnpInfo(data, opt.bedFile,opt.geneInfoFile,opt.cisRegionWind, 
                            opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile, opt.includeChr, 
                            opt.excludeAmbiguousSNP, opt.skeletonSnpFile, opt.geneticMapFile, opt.ldBlockInfoFile, 
                            opt.includeBlock, opt.flank, opt.mafmin, opt.mafmax, opt.noscale, readGenotypes);
                    data.jackknifeLDmatrix(opt.ldmatrixFile, opt.outLDmatType, opt.title, opt.writeLdmTxt);
                }
                else if (opt.binSnp) {
                    omics.inputSnpInfo(data, opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile, "", opt.ldmatrixFile, opt.includeChr, opt.excludeAmbiguousSNP, opt.skeletonSnpFile, opt.geneticMapFile, opt.genMapN, opt.flank, opt.eQTLFile, opt.ldscoreFile, opt.windowFile, opt.multiLDmat, opt.excludeMHC, opt.afDiff, opt.mafmin, opt.mafmax, opt.pValueThreshold, opt.rsqThreshold, opt.sampleOverlap, opt.imputeN, opt.noscale, opt.binSnp, opt.readLdmTxt);
                    data.binSnpByLDrsq(opt.rsqThreshold, opt.title);
                }
                else {
                    omics.inputSnpInfo(data, opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile, "", opt.ldmatrixFile, opt.includeChr, opt.excludeAmbiguousSNP, opt.skeletonSnpFile, opt.geneticMapFile, opt.genMapN, opt.flank, opt.eQTLFile, opt.ldscoreFile, opt.windowFile, opt.multiLDmat, opt.excludeMHC, opt.afDiff, opt.mafmin, opt.mafmax, opt.pValueThreshold, opt.rsqThreshold, opt.sampleOverlap, opt.imputeN, opt.noscale, opt.binSnp, opt.readLdmTxt);
                    data.resizeLDmatrix(opt.outLDmatType, opt.chisqThreshold, opt.windowWidth, opt.LDthreshold, opt.effpopNE, opt.cutOff, opt.genMapN);
                    data.outputLDmatrix(opt.outLDmatType, opt.title, opt.writeLdmTxt);
                }
            }
        }  // end of analysisType: LDmatrix
        else if (opt.analysisType == "LDmatrixEigen") {
            readGenotypes = false;
            /***** Step 2. perform eigen decomposition for gene cis-region LD matrices *****/
            bool calcGeneLDBool = !opt.eqtlSummaryFile.empty() || !opt.eqtlSummaryQueryFile.empty() || !opt.geneListFile.empty();
            if(!opt.ldBlockInfoFile.empty() || calcGeneLDBool){
                /***** Step 1. read individual and genotype infomation *****/
                bool generateLDBool = true;
                omics.inputIndInfo(data, opt.bedFile, "", opt.keepIndFile, opt.keepIndMax,
                    opt.mphen, opt.covariateFile, opt.randomCovariateFile, opt.residualDiagFile,opt.geneInfoFile,opt.keepIndGeneFile,opt.subGenePath,generateLDBool);
                omics.inputSnpInfo(data, opt.bedFile,opt.geneInfoFile,opt.cisRegionWind, 
                    opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile, opt.includeChr, 
                    opt.excludeAmbiguousSNP, opt.skeletonSnpFile, opt.geneticMapFile, opt.ldBlockInfoFile, 
                    opt.includeBlock, opt.flank, opt.mafmin, opt.mafmax, opt.noscale, readGenotypes);
                /***** Step 2. read LD matrix at once *****/
                data.readBedFile(opt.noscale,opt.bedFile + ".bed");
                /***** Step 3. Do eigen-decomposition for LD Blocks or gene cis-region *****/
                if(!opt.ldBlockInfoFile.empty()) data.calcBlockLDEigenDecom(opt.bedFile + ".bed",opt.ldBlockInfoFile,0,opt.title, opt.eigenCutoff.maxCoeff(), opt.outInfoOnly);
                if(calcGeneLDBool ) {
                    data.calcGeneLDEigenDecomBasedBesd(opt.eqtlSummaryFile,opt.eqtlSummaryQueryFile,opt.geneListFile,opt.includeSpecificGeneID, opt.title,opt.cisRegionWind, opt.eigenCutoff.maxCoeff(),opt.diagnosticMode); 
                }             
            }else if(!opt.eigenMatrixFile.empty()) { // merge existing eigen matrices
                //omics.inputSnpInfo(data, opt.includeSnpFile, opt.excludeSnpFile, opt.excludeRegionFile, "", opt.eigenMatrixFile, opt.ldBlockInfoFile, opt.includeChr, opt.excludeAmbiguousSNP, opt.flank, opt.eQTLFile, opt.ldscoreFile, opt.eigenCutoff, opt.excludeMHC, opt.afDiff, opt.mafmin, opt.mafmax, opt.pValueThreshold, opt.rsqThreshold, opt.sampleOverlap, opt.imputeN, opt.noscale, opt.readLdmTxt);
                /***** Step 1. perform eigen decomposition for the blocked LD matrices *****/
                   
            } else{
                LOGGER.e(0,"Something is wrong dur generating or merging low-rank LD process. \nPlease see (" + string(DOC_ONLINE) + ") for details.");
            }
        } // end of analysisType: LDmatrixEigen
        else if (opt.analysisType == "MergeGwasSummary") {
            if (opt.outLDmatType == "block") {
                data.mergeBlockGwasSummary(opt.gwasSummaryFile, opt.title);
            }
        } // end of analysisType: MergeGwasSummary
        else {
            LOGGER.e(0,"Wrong analysis type: " + opt.analysisType, "");
        } // end of analysisType

        /******************************************************/
        LOGGER.i(0, "");
        timer.getTime();
        LOGGER.i(0,  "at " + timer.getDate(), "Analysis finished");
        LOGGER << "Computational time: "  << timer.format(timer.getElapse()) << endl;

    } catch (const string &err_msg) {
        LOGGER.e(0, err_msg);
    } catch (const char *err_msg) {
        LOGGER.e(0, string(err_msg));
    } // end of try function
    return 0;
} // End of main function