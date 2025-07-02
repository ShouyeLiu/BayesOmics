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



#ifndef sbayesE_hpp
#define sbayesE_hpp

#include <stdio.h>
#include "Data.hpp"
#include "Mcmc.hpp"
#include "Model.hpp"

class ApproxBayesCO : public ApproxBayesC {
    public:
    class EQTLJointMat : public ParamMat {
        public:
        vector<string> geneNames;
        vector<string> cisSnpIDVec;
        unsigned numGenes;
        unsigned numSnpInCis;
        EQTLJointMat(const vector<string> &cisSnpIDVec,
        const vector<string> &geneNames,
        const string &lab = "EQTLJointMat"): ParamMat(lab,cisSnpIDVec,geneNames),
        numSnpInCis(int(cisSnpIDVec.size())),
        numGenes(int(geneNames.size())){   
        }
    };
    class SnpEffectMat : public ParamMat {
        public:
        vector<string> geneNames;
        vector<string> snpNames;
        unsigned numGenes;
        unsigned numSnps;
        SnpEffectMat(const vector<string> &snpNames,
        const vector<string> &geneNames,
        const string &lab = "SnpEffectMat"): ParamMat(lab,snpNames,geneNames),
        numSnps(int(snpNames.size())),
        numGenes(int(geneNames.size())){   
        }
    };

    class EQTLJointVec : public vector<ParamSet*> {
        public:
        vector<string> colnames;
        vector<string> geneNames;
        map<int, vector<string>> gene2cisSnpIDMap;
        unsigned numGenes;
        EQTLJointVec(const vector<string> &geneNames,
        const  map<int, vector<string>> &gene2cisSnpIDMap,
        const string &lab = "EQTLJointVec"):
        numGenes(int(gene2cisSnpIDMap.size())),geneNames(geneNames),gene2cisSnpIDMap(gene2cisSnpIDMap) {
            colnames.resize(numGenes);
            for(unsigned i = 0; i< numGenes; i++){
                colnames[i] = lab + geneNames[i];
                this->push_back(new ParamSet(colnames[i], gene2cisSnpIDMap.at(i)));
            }
        } 
    };

    class SnpEffectVec : public vector<ParamSet*> {
        public:
        vector<string> colnames;
        vector<string> gwasAndgeneNames;
        map<int,vector<string>> gwas2SnpIDMap;
        unsigned numGenes;// here including gwas and genenames
        SnpEffectVec(const vector<string> &gwasAndgeneNames,
            const map<int,vector<string>> &gwas2SnpIDMap,
        const string &lab = "SnpEffectVec"):
        numGenes(int(gwas2SnpIDMap.size())),gwasAndgeneNames(gwasAndgeneNames),gwas2SnpIDMap(gwas2SnpIDMap){
            colnames.resize(numGenes);
            for(unsigned i = 0; i< numGenes; i++){
                colnames[i] = lab + gwasAndgeneNames[i];
                this->push_back(new ParamSet(colnames[i], gwas2SnpIDMap.at(i)));
            }
        } 
    };

    class SigmaSqBetaNonEqtl : public Parameter, public Stat::InvChiSq {
        // variance of snp effects has a scaled-inverse chi-square prior
    public:
        double df;  // hyperparameter nub;
        double scale;     // hyperparameter
        double vargPrior;
        int numNonEqtl;
        bool noscale;
        double scaleRunningMean; // for running mean
        double intergenicRatio;

        SigmaSqBetaNonEqtl(const double vg, const unsigned numSnps, const unsigned numNonEqtl, 
        const double piEffNonEqtl, const bool noscale, const string &lab = "varbNEqtl")
        :vargPrior(vg),numNonEqtl(numNonEqtl), Parameter(lab), df(4){
            if(numNonEqtl != 0){
                intergenicRatio = (double) numNonEqtl/ (double) numSnps;
                // value = vg / ( snp2pq.mean() * numNonEqtl * piEffNonEqtl);  // derived from prior knowledge on Vg and pi
                value = vargPrior* intergenicRatio /((double) numNonEqtl * piEffNonEqtl);

                scale = 0.5 * value;
                scaleRunningMean = scale;
                // LOGGER << "vargPrior: " << vargPrior << " numNonEqtl: " << numNonEqtl << " numSnps: " << numSnps << " piEffNonEqtl: " << piEffNonEqtl << endl;
            }
        }
        void sampleFromFC(const double &snpEffSumSq, const unsigned &numSnpEff);
        void sampleFromFC(const double &snpEffSumSq, const unsigned &numSnpEff,const unsigned hsqBetaInt);
        void sampleFromFC(const double &snpEffSumSq, const unsigned &numSnpEff,const double varg,const double piEffEqtl,const unsigned iter, const unsigned burnIn);
    };

    class SigmaSqBetaEqtl : public Parameter, public Stat::InvChiSq {
        // variance of snp effects has a scaled-inverse chi-square prior
    public:
        const double df;  // hyperparameter
        double scale;     // hyperparameter
        bool noscale;  // no scaling on the genotypes
        double vargPrior;
        int numEqtl;
        double scaleRunningMean;
        double genicRatio;

        SigmaSqBetaEqtl(const double vg, const unsigned numSnps, const unsigned numEqtl, const double piEffEqtl, 
        const bool noscale, const string &lab = "varbEqtl")
        :vargPrior(vg),numEqtl(numEqtl), Parameter(lab), df(4), noscale(noscale){
            if(numEqtl != 0){
                genicRatio = (double) numEqtl/ (double) numSnps;
                value = vg * genicRatio /( (double) numEqtl * piEffEqtl );
                scale = 0.5f * value;  
                scaleRunningMean = scale;
            }            
        }
        // void sampleFromFC(const double snpEffSumSq, const unsigned numSnpEff);
        // void sampleFromPrior(const double piEffEqtl);
        void sampleFromFC(const double &snpEffSumSq, const unsigned &numSnpEff,const double varg,const double piEffNonEqtl,const unsigned iter, const unsigned burnIn);
        void sampleFromFC(VectorXd ssqAlphaEqtlPG, const VectorXd &numNonZerosEqtlVecAcrossGenesPostIW);
    };

    class SigmaSqAlphaVec : public ParamSet,public Stat::InvChiSq, public vector<ApproxBayesC::VarEffects*> {
        public:
        VectorXd values;
        vector<string> geneNames;
        int numGenes;
        double df;
        VectorXd scaleVec;
        VectorXd scaleRunMeanVec;
        VectorXd numEqtlPG;

        SigmaSqAlphaVec(const vector<string> geneNames,const VectorXd varPhenotypiceQTL, const double h2eQTL,
                const map<int,vector<int>> &gene2cisSnpMap, 
                const double piEffEqtl,const bool noscale):
            ParamSet("SigmaSqAlphaVec", geneNames)
            , geneNames(geneNames),
            numGenes(geneNames.size()){
            values.setZero(numGenes);
            scaleVec.setZero(numGenes);
            numEqtlPG.setZero(numGenes);
            scaleRunMeanVec.setZero(numGenes);
            for(unsigned i = 0; i < geneNames.size(); ++i){
                vector<int> snpsInGene = gene2cisSnpMap.at(i);
                numEqtlPG[i] = snpsInGene.size();
                values[i] = varPhenotypiceQTL[i] * h2eQTL /( numEqtlPG[i] * piEffEqtl );
                scaleVec[i] = values[i]* 0.5;
                scaleRunMeanVec[i] = scaleVec[i];
            }
            df = 4;
        }
        void sampleFromPrior();
        void sampleFromFC(VectorXd ssqAlphaEqtlPG, const VectorXd &numNonZerosEqtlVecAcrossGenesPostIW,const VectorXd varg,const double piEffEqtl,const unsigned iter, const unsigned burnIn);
        
    };

    class SigmaSqTheta : public Parameter, public Stat::InvChiSq {
        // variance of snp effects has a scaled-inverse chi-square prior
    public:
        const double df;  // hyperparameter
        double scale;     // hyperparameter

        SigmaSqTheta(const double h2snp,const double h2eQTL, const unsigned numGenes, const double piTheta, const string &lab = "SigmaSqTheta")
        :Parameter(lab), df(4){
            value =  0.5f * h2snp/(numGenes * piTheta * h2eQTL);  // derived from prior knowledge
            // value = 1.0;
            scale = 0.5f * value;
            // cout << "SigmaSqTheta: " << value << " scale: " << scale << " piTheta: " << piTheta << " h2snp: " << h2snp;
            // cout << " h2eQTL " << h2eQTL << endl;
        }
        void sampleFromFC(const double snpEffSumSq, const unsigned numSnpEff);
    };

    class SigmaSqMat : public vector<ParamMat*>, public Stat::Gamma, public Stat::InvChiSq {
        public:
        vector<string> geneNames; // sigmaSqMat is symmatic matix 
        unsigned numGenes;
        vector<MatrixXd> sigmaSqMats;
        vector<MatrixXd> sigmaSqInvMats;
        vector<MatrixXd> varcovPriors;
        vector<MatrixXd> iwScaleMat; // average scale matrix as a whole
        // average scale matrix separately
        vector<MatrixXd> iwScaleMatCorr; // for runnign mean
        VectorXd scaleMatBetaSigmaSquareVec; // for running mean
        VectorXd scaleMatAlphaSigmaSquareVec; // for running mean
        /// ///////////
        VectorXd sigmaSqDetLogVec;
        VectorXd sigmaSqAlphaPM; // sigmaSq Alpha per xQTL
        VectorXd scaleAlphaPM; 
        unsigned numNonZeroTheta;
        // VectorXd geneEffects;
        // double nnG;
        double sigmaSqBetaEqtlPM;
        double sigmaSqAlphaAll;
        // priors
        double scaleBetaEqtl, nub,nua, scaleAlphaAll;

        SigmaSqMat (const vector<string> &geneNames,const string &lab = "SigmaSqMats"):
        numGenes(int(geneNames.size()))
        , geneNames(geneNames){
            sigmaSqMats.resize(numGenes);
            sigmaSqInvMats.resize(numGenes);
            varcovPriors.resize(numGenes);
            iwScaleMat.resize(numGenes);
            iwScaleMatCorr.resize(numGenes);
            // geneEffects.resize(numGenes);
            scaleMatBetaSigmaSquareVec.resize(numGenes);
            scaleMatAlphaSigmaSquareVec.resize(numGenes);
            sigmaSqAlphaPM.resize(numGenes);
            scaleAlphaPM.resize(numGenes);
            
            sigmaSqDetLogVec.setZero(numGenes);
            for(unsigned i = 0; i < numGenes; ++i){
                vector<string> sigmaSqName = {"gwas",geneNames[i]};
                this->push_back(new ParamMat(geneNames[i],sigmaSqName,sigmaSqName));
                sigmaSqMats[i].setZero(2,2);
                sigmaSqInvMats[i].setZero(2,2);
                varcovPriors[i].setZero(2,2);
                iwScaleMat[i].setZero(2,2);
                iwScaleMatCorr[i].setZero(2,2);
                // sigmaSqDetLogVec.setZero();
            }
            scaleBetaEqtl = 0;
            scaleAlphaAll = 0;
            sigmaSqBetaEqtlPM =0;
            sigmaSqAlphaAll =0;
            nub = 4.0;
            nua = 4.0;
            scaleMatBetaSigmaSquareVec.setZero();
            scaleMatAlphaSigmaSquareVec.setZero();
            sigmaSqAlphaPM.setZero();
            scaleAlphaPM.setZero();
        }
        void sampleFromFCIWIndSMatPriorAIAO(int iter, int burnIn, VectorXd &geneEffects,const double &ssqBetaEqtl,const double &ssqAlphaEqtl, 
        const double & ssqBetaTotalGenic, const unsigned &numNonNullBetaTotGenic,
         const  VectorXd ssqBetaEqtlPG,const  VectorXd ssqAlphaEqtlPG, vector<Matrix2d> ssqEqtlMat,
         const VectorXd &numNonZerosEqtlVec,const VectorXd &numNonZerosEqtlVecAcrossGenesPostIW,
         const double hsqGenic,const VectorXd cisHsqMean ,const bool messageBool);
        void sampleFromFCIWSMatCorrPriorAIAO(int iter, int burnIn,double sigmaSqBetaEqtlPM,VectorXd sigmaSqAlphaPM, VectorXd &geneEffects,const double &ssqBetaEqtl,const double &ssqAlphaEqtl, 
         const double & ssqBetaTotalGenic, const unsigned &numNonNullBetaTotGenic,
         const  VectorXd ssqBetaEqtlPG,const  VectorXd ssqAlphaEqtlPG, vector<Matrix2d> ssqEqtlMat,
         const VectorXd &numNonZerosEqtlVec , const VectorXd &numNonZerosEqtlVecAcrossGenesPostIW,
         const double hsqGenic,const double cisHsqMean ,const bool messageBool);

        void sampleFromFCIWIndSMatPriorEIEO(int iter, int burnIn, const double &sigmaSqBetaEqtlPM,const VectorXd &sigmaSqAlphaPM,const vector<Matrix2d> & ssqEqtlMat, double &geneticCov, 
         VectorXd &geneEffects, const VectorXd &numEqtlPG,
         const bool messageBool);

        void setPrior(const double &sigmaSqBetaEqtl, const VectorXd &sigmaSqAlphaVec);
       // VectorXd compute
    };
   
     class SigmaSqMatResidual: public vector<ParamMat*>, public Stat::Gamma, public Stat::InvChiSq {
        public:
        vector<string> geneNames; // sigmaSqMat is symmatic matix 
        unsigned numGenes;
        vector<MatrixXd> sigmaSqMats;
        vector<MatrixXd> sigmaSqInvMats;
        vector<MatrixXd> varcovPriors;
        double scaleBetaEqtl, nub,nua, scaleAlphaAll;

        SigmaSqMatResidual (const vector<string> &geneNames,const string &lab = "SigmaSqMatsRes"):
        numGenes(int(geneNames.size()))
        , geneNames(geneNames){
            sigmaSqMats.resize(numGenes);
            sigmaSqInvMats.resize(numGenes);
            varcovPriors.resize(numGenes);
            for(unsigned i = 0; i < numGenes; ++i){
                vector<string> sigmaSqName = {"gwas",geneNames[i]};
                this->push_back(new ParamMat(geneNames[i],sigmaSqName,sigmaSqName));
                sigmaSqMats[i].setZero(2,2);
                sigmaSqInvMats[i].setZero(2,2);
                varcovPriors[i].setZero(2,2);
            }
            nub = 4.0;
            nua = 4.0;
        }
        void sampleFromFC(const vector<VectorXd> &wbcorrGene,const vector<VectorXd> &wAcorrGene,const VectorXd &varEps, const VectorXd &vareVec, 
        const map<int,vector<int> > &ldblock2gwasSnpMap, const map<string, int> &geneID2IdxMap,
        const map<int,string> &gwasSnpIdx2snpIDMap,const map<string,vector<string> > &gwasSnpID2geneIDMap);
        void setPrior(const double &gwasVare, const VectorXd &varEps);
       // VectorXd compute
    };
    class SnpEffects : public ApproxBayesC::SnpEffects, public Stat::MultiNormal {
    public:
        VectorXd betaTotal;     // save sample squres of full conditional normal distribution regardless of delta values
        VectorXd betaTotalLatent;
        VectorXd snpEffectVec;
        Vector2d numNonZerosEqtlVec;    // number of nonzero effects in genic regions
        VectorXd numNonZerosEqtlVecAcrossGenesPostIW; 
        VectorXd numNonZerosGenicGwasVec; // used to cauculate genic gwas enrichment.
        VectorXd numNonZerosGenicEqtlVec; // used to cauculate genic gwas enrichment.
        unsigned numNonZerosNonEqtl;    // number of nonzero effects in between gene regions
        unsigned numNonNullEqtl;  // number of SNPs with at least one nonzero effect in genic regions
        unsigned numNonNullSnpTot;  // total number of SNPs with nonzero effects in the genome
        unsigned numNonNullBetaTotGenic; // number of SNPs with at least one nonzero gwas effect
        double numNonNullSnpPerGene;
        double nnG;                  // number of genes with at least one null zero eQTLs.
        Vector4d numSnpCompVec; // store for nsnp00-11 
        // zhili's pipeline
        VectorXd ssqBlocks;
        vector<VectorXd> whatBlocks;
        double vargInt;
        double vargGenic;
        vector<VectorXd> whatBlocksInt;
        vector<VectorXd> whatBlocksGen;

        // snp effect
        vector<Matrix2d> ssqEqtlMat;
        double ssqNonEqtl;
        double ssqBetaEqtl;
        double ssqAlphaEqtl;
        double ssqBetaTotalGenic; // ssq for betaTotal in the genic region
        VectorXd ssqAlphaEqtlPG; // ssq for alpha per gene
        VectorXd ssqBetaEqtlPG;


        VectorXd betaTotalMean; // debug
        // calculate heritability directly
        VectorXd ghat;  // used to calculate total snp heritability
        map<int, VectorXd> gwhatMap; // used to calculate 
        map<int, VectorXd> gwhatGwasMap; // used to calculate mediated heritability

        SnpEffects(const vector<string> &header): ApproxBayesC::SnpEffects(header){
            betaTotal.setZero(size);
            betaTotalLatent.setZero(size);
            numNonZerosEqtlVec.setZero(2);
            ssqNonEqtl = 0;
            ssqBetaEqtl = 0;
            ssqAlphaEqtl = 0;
            ssqBetaTotalGenic = 0;
            badSnps.setZero(size);

        }                                          
        void sampleFromFCAIAO(Data data, const vector<MatrixXd> &QblocksMat, const int iter,const bool diagnose, const string title,
                vector <VectorXd> &wcorrBlocks, vector<VectorXd> &wAcorr,vector<VectorXd> &wbcorrGene,
                const vector <MatrixDat> &Qblocks, const vector<MatrixDat> &Qgene, 
                const map<int,vector<int> > &ldblock2gwasSnpMap,
                SnpEffectVec &snpEffectVec, EQTLJointVec &eQTLJointVec,
                const map<int,string> &gwasSnpIdx2snpIDMap, const map<string, int> &geneID2IdxMap, const map<string,vector<string> > &gwasSnpID2geneIDMap, 
                const map<string, int> &cisSnpID2IdxMap, const double &sigmaSqBetaNonEqtl,  SigmaSqMat &sigmaSqMats,SigmaSqMatResidual &sigmaSqMatRes, const double &piEffEqtl, const double &piEffNonEqtl,
                const VectorXd &nGWAS, const vector<VectorDat> &neQTL,const VectorXd &varEps, const VectorXd &vareVec);

        void sampleFromFCEIEO(Data data,const vector<MatrixXd> &QblocksMat,const int iter,const bool diagnose, const string title, vector <VectorXd> &wcorrBlocks, 
                vector<VectorXd> &wAcorr, const vector <MatrixDat> &Qblocks, const vector<MatrixDat> &Qgene, const map<int,vector<int> > &ldblock2gwasSnpMap, 
                SnpEffectVec &snpEffectVec, EQTLJointVec &eQTLJointVec, SnpEffectVec &snpEffectVecLatent, EQTLJointVec &eQTLJointVecLatent,
                const map<int,string> &gwasSnpIdx2snpIDMap, const map<string, int> &geneID2IdxMap, const map<string,vector<string> > &gwasSnpID2geneIDMap, 
                const map<string, int> &cisSnpID2IdxMap, const double &sigmaSqBetaNonEqtl,  SigmaSqMat &sigmaSqMats, const double &piEffEieo1, 
                const double &piEffEieo2, const double &piEffNonEqtl, const VectorXd &nGWAS, const vector<VectorDat> &neQTL,const VectorXd &varEps, const VectorXd &vareVec);
    };

    class GeneEffects :public Stat::InvChiSq, public Stat::Normal, public ParamSet{
        public:
            vector<string> geneNames;
            unsigned numGenes;
            double propMed;
            double vareMed;
            unsigned nnGene;
            double  ssqGene;
            VectorXd deltaGene;
            VectorXd valuesAdjust;
            VectorXd geneEffectScaleFactorMean ;
            GeneEffects(const vector<string> geneNames):
            ParamSet("GeneEffects",geneNames),
            numGenes(int(geneNames.size())),
            geneNames(geneNames){
                propMed = 0.0;
                nnGene = 0;
                deltaGene.setZero(numGenes);
                valuesAdjust.setZero(numGenes);
                geneEffectScaleFactorMean.setZero(numGenes);
                vareMed = 0.01;
            }

        void sampleFromeFC(Data data, int iter,const string title ,VectorXd &betaTotal, MatrixXd &eQTLJointMat,
                            const map<int,vector<int> > &gene2gwasSnpMap,const map<int,vector<int>> &gene2cisSnpMap,
                            double &sigmaSqTheta, double &vareMed, const double &piTheta, VectorXd &deltaGene);

        void sampleFromeFC(Data data, int iter ,const string title, VectorXd &betaTotal, EQTLJointVec &eQTLJointVec,
                            const map<int,vector<int> > &gene2gwasSnpMap,
                            const map<int,string> &gwasSnpIdx2snpIDMap,
                            double &sigmaSqTheta, double &vareMed, const double &piTheta, VectorXd &deltaGene );
    };
   //////// heritability related parameter
    class MediatedHeritability : public Parameter {
        // compute heritability based on sampled values of genotypic and residual variances
        // strictly speaking, this is not a model parameter
        public:
        MediatedHeritability(const string &lab = "hsqMed"): Parameter(lab){};
        void compute(const double &geneVar, const double vary){
            value = geneVar/vary;
        }
    };
    class Heritability : public Parameter {
        // compute heritability based on sampled values of genotypic and residual variances
        // strictly speaking, this is not a model parameter
    public:
        double valueGenic;
        Heritability(const string &lab = "hsq"): Parameter(lab){
            valueGenic = 0;
            
        };
        void compute(const double genVar, const double vary){
            value = genVar/vary;
        }
    };

    class HeritabilityCis : public ParamSet {
        // heritability from gene expression to SNPs
    public:
         VectorXd valueGenicVec;
        // double cisHsqMean;
        HeritabilityCis(const vector<string> &header, const string &lab = "cisHsq"): ParamSet(lab,header){
            valueGenicVec.setZero(header.size());
        };
        void compute(const VectorXd geneVar, const VectorXd varyQTl){
            values = geneVar.array() / varyQTl.array();
        }
    };
    class HeritabilityGwasCis : public ParamSet {
        // partition GWAS heritability to different genes
    public:
        // double cisHsqMean;
        HeritabilityGwasCis(const vector<string> &header, const string &lab = "gwasCisHsq"): ParamSet(lab,header){};
        void compute(const VectorXd geneVar, const double vary){
            values = geneVar.array() / vary;
        }
    };

    class HeritabilityEnrich : public ParamSet {
        public:
        HeritabilityEnrich(const vector<string> &header, const string &lab = "HsqEnich"):
        ParamSet(lab,header){
            values.setZero(header.size());
        };
        void compute(const VectorXd hsqVec, const VectorXd numNonZerosVec){
            values = hsqVec.array() / numNonZerosVec.array();
        }
    };

    class GenotypicVar : public ApproxBayesC::GenotypicVar {
        public:
        GenotypicVar(const double varg, const unsigned n): ApproxBayesC::GenotypicVar(varg,n){}
        void compute(VectorXd &betaTotal, const vector <MatrixDat> &Qblocks,const map<int,vector<int> > &ldblock2gwasSnpMap);
        void compute(int niter,VectorXd &betaTotal, const vector<MatrixXd> &Qblocks,const map<int,vector<int> > &ldblock2gwasSnpMap);
    };

    class GenotypicVarGene : public Parameter {
        public:
        unsigned numGenes;
        GenotypicVarGene(const vector<string> &header,const string &lab = "vargGene"):numGenes(header.size()), Parameter(lab){
        };
        void compute(VectorXd &geneEffects, EQTLJointVec &eQTLJointVec,const vector<MatrixDat> &Qgene);
    };

    class GenotypicVarGeneCis : public ParamSet {
        public:
        VectorXd valueBetaVec;
        unsigned numGenes;
        GenotypicVarGeneCis(const vector<string> &header, const string &lab = "vargGeneCis"): numGenes(header.size()),
        ParamSet(lab,header){
            valueBetaVec.setZero(header.size());
        };
        void compute( EQTLJointVec &eQTLJointVec,const vector<MatrixDat> &Qgene);
        void computeHsqBeta(SnpEffectVec &snpEffectVec,const vector<MatrixDat> &Qgene);

    };
    
    class ResidualVar : public Stat::InvChiSq{
        public: 
        double df;  // hyperparameter
        double scale;     // hyperparameter
        VectorXd valueVec; 
        double value;   
        unsigned numBlocks;

        double vary;
        float threshold; // sbayesrc's way;
        float mean;   
        ResidualVar(const map<int,vector<int> > &ldblock2gwasSnpMap,const double &vary):numBlocks(ldblock2gwasSnpMap.size()),vary(vary){
            valueVec.setZero(numBlocks);
            for(int i = 0; i < numBlocks; i++){
                valueVec[i] = vary ;
            }
            value = vary;
            df = 4.0;
            scale = 0.5 * vary;
            threshold = 1.1;
            mean = vary;

        }
        void sampleFromFC(int iter,const double &varPhenotypic, const vector <VectorXd> &wcorrBlocks, const vector<VectorXd> &gwasEffectInBlock, 
                                              const vector <MatrixDat> &Qblocks, const VectorXd &nGWAS, const VectorXd &betaTotal,
                                              const map<int,vector<int> > &ldblock2gwasSnpMap);
        void sampleFromFC(int iter,vector<VectorXd> &wcorrBlocks, const vector<VectorXd> &whatBlocks, VectorXd &ssqBlocks, const VectorXd &nGWASblocks, const VectorXd &numEigenvalBlock);

    };

    class Rounding : public ApproxBayesC::Rounding {
    public:
        void computeRcorr(vector<MatrixXd> &LDBlockUs,vector<VectorXd> &LDBlockLambdas, 
        const vector<ChromInfo*> &chromInfoVec,const VectorXd &snpEffects, VectorXd &rcorr);
    };

    class ResidualVareEQTL : public ParamSet, public vector<ApproxBayesC::ResidualVar*>,public Stat::InvChiSq {
        public:
        // VectorXd values;
        vector<string> geneNames;
        unsigned numGenes;
        double df;  // hyperparameter
        // double scale;     // hyperparameter
        VectorXd scales;
        
        ResidualVareEQTL(const vector<string> geneNames,const VectorXd &varGenotypiceQTL ,const VectorXd vare, const unsigned nobs, const double icrsq):
            ParamSet("ResidualVareEQTL",geneNames),numGenes(int(geneNames.size())), geneNames(geneNames){
            scales.resize(numGenes);
            for(unsigned i = 0; i < numGenes; ++i){
                // this->push_back(new ApproxBayesC::ResidualVar(vare[i],nobs,icrsq));
                // values[i] = (*this)[i]->value;
                values[i] = varGenotypiceQTL(i);
                scales[i] = 0.5 * values[i];
            }   
            df = 4.0;         
        }
        void sampleFromFC(const VectorXd &varPhenotypiceQTL,const vector<VectorXd> &wAcorr, const vector<VectorXd> &eQTLEffAcrossGenes, 
                                             const vector<MatrixDat> &Qgene, const vector<VectorDat> &neQTL, EQTLJointVec &eQTLJointVec);
    };


    public:

    string mcmcType;
    bool eieoLatent;
    bool sampleVareBool;
    bool sampleVarEpsBool;
    bool message;
    Stat::Normal normal;

    vector<VectorXd> wAcorr;
    vector<VectorXd> wbcorrGene;
    vector <VectorXd> wcorrBlocks;
    
    MediatedHeritability medHsq;
    Heritability hsq;
    HeritabilityCis cisHsq;
    HeritabilityGwasCis gwasCisHsq;
    SnpEffects snpEffects;
    // EQTLJointMat eQTLJointMat;
    // EQTLJointMat eQTLJointMatLatent;
    // SnpEffectMat snpEffectMat;
    // SnpEffectMat snpEffectMatLatent;
    // effects
    EQTLJointVec eQTLJointVec;
    EQTLJointVec eQTLJointVecLatent;
    SnpEffectVec snpEffectVec;
    SnpEffectVec snpEffectVecLatent;


    SigmaSqBetaNonEqtl sigmaSqBetaNonEqtl;
    SigmaSqBetaEqtl  sigmaSqBetaEqtl;
    SigmaSqTheta sigmaSqTheta;
    
    // VarSnpEffects sigmaSq;  // single value for gwas
    SigmaSqAlphaVec sigmaSqAlphaVec;   // vector vlaues for various genes
    SigmaSqMat sigmaSqMats;  // matrix values for various genes
    SigmaSqMatResidual sigmaSqMatRes; // sample residuals
    GeneEffects geneEffectVec;
    // EQTLJointVec eQTLJointVec;
    ResidualVareEQTL varEps; // residuals for eQTL of genes
    ResidualVar vare;
    GenotypicVar varg;
    GenotypicVarGene vargGene;
    GenotypicVarGeneCis vargGeneCis;

    HeritabilityEnrich genicGwasEnrich;
    HeritabilityEnrich genicEqtlEnrich;

    Rounding  rounding;
    BayesC::Pi piEffEqtl;   // pi for nonzero effects in genic regions (overlaping eQTL are allowed to have multiple nonzero effects)
    BayesC::Pi piEffEieo1; // used in EIEO model
    BayesC::Pi piEffEieo2; // used in EIEO model
    BayesC::Pi piTheta;
    BayesC::Pi piEffNonEqtl;   // pi for nonzero effects in between gene regions (thus the same of nonnull SNPs in intergenic regions)
    NumNonZeroSnp nnzBtw;  // number of nonzero effects in between gene regions
    NumNonZeroSnp nnzGen;  // number of nonzero effects in genic regions for both gwas and eQTL total 11
    NumNonZeroSnp nnsGen;  // number of SNPs with at least one nonzero effect in genic regions for eQTL
    NumNonZeroSnp nnsTot;  // total number of SNPs with nonzero effects in the genome for GWAS 
    NumNonZeroSnp nnsPG;   // average number of SNPs with nonzero effects per gene
    NumNonZeroSnp nnGene;     // nonzero gene effect;
    // NumNonZeroSnp nnzTheta;
    // eieo model
    NumNonZeroSnp nnEqtl;
    NumNonZeroSnp nnEqtlOverlap; // number of beta effect in the genic region
    NumNonZeroSnp nsnp00; // no gwas and no eqtl
    NumNonZeroSnp nsnp10; // have gwas and no eqtl
    NumNonZeroSnp nsnp01; // no gwas and have eqtl
    NumNonZeroSnp nsnp11; //  have gwas and have eqtl

    Parameter cisHsqMean;  // mean value of cis-heritability for gene expression
    Parameter vareMean;    // mean value of trait residual variance across LD blocks
    Parameter varEpsMean;  // mean value of expression residual variance across genes
    Parameter sigmaSqAlpha;
    Parameter geneticCorr;
    // calculate running mean
    // double genVarNonEqtlPrior; // scale for intergenic region
    // double scaleNonEqtlPrior;



    ApproxBayesCO(const Data &data,const string mcmcType,const bool eieoLatent,  const bool sampleVareBool, const bool sampleVarEpsBool, const double varGenotypic, const double varResidual, const double varRandom ,
                const double h2snp, const double h2eQTL,const double pival,const double piEffEqtlVal, const double piGenicGwas,const double piGenicEqtl, const double piThetaVal,
                const double piEffNonEqtlVal, const double piAlpha, const double piBeta, const bool estimatePi, const bool noscale,
                const double phi, const double overdispersion, const bool estimatePS, const double icrsq, const double spouseCorrelation,
                const bool diagnosticMode, const bool robustMode, const bool randomStart = true, const bool message = true)
    : ApproxBayesC(data,true, varGenotypic, varResidual,varRandom, pival, piAlpha, piBeta, estimatePi, noscale, phi, overdispersion, estimatePS, icrsq, spouseCorrelation, diagnosticMode, robustMode, randomStart, false)
    , wcorrBlocks(data.wcorrBlocks)
    , wAcorr(data.wAcorr)
    , wbcorrGene(data.wbcorrGene)
    , snpEffects(data.snpEffectNames)
    , geneEffectVec(data.geneEffectNames)
    // , snpEffectMat(data.snpEffectNames,data.gwasAndGeneEffectNames)
    // , snpEffectMatLatent(data.snpEffectNames,data.gwasAndGeneEffectNames)
    // , eQTLJointMat(data.cisSnpIDVec,data.geneEffectNames)
    // , eQTLJointMatLatent(data.cisSnpIDVec,data.geneEffectNames)
    , eQTLJointVec(data.geneEffectNames,data.gene2cisSnpIDMap,"EQTLJointVec_")
    , eQTLJointVecLatent(data.geneEffectNames,data.gene2cisSnpIDMap,"EQTLJointVecLatent_")
    , snpEffectVec(data.gwasAndGeneEffectNames,data.gwas2SnpIDMap,"SnpJointVec_")
    , snpEffectVecLatent(data.gwasAndGeneEffectNames,data.gwas2SnpIDMap,"SnpJointVecLatent_")
    , piEffEqtl(piEffEqtlVal, piAlpha, piBeta,"piEffEqtl")
    , piTheta(piThetaVal,piAlpha,piBeta,"piTheta")
    , piEffNonEqtl(piEffNonEqtlVal, piAlpha, piBeta,"piEffNonEqtl")
    , piEffEieo1(piGenicGwas, piAlpha, piBeta,"piEffEieo1")
    , piEffEieo2(piGenicEqtl, piAlpha, piBeta,"piEffEieo2")
    , mcmcType(mcmcType)
    , eieoLatent(eieoLatent)
    , message(message)
    , sampleVareBool(sampleVareBool)
    , sampleVarEpsBool (sampleVarEpsBool)
    , nnzGen("NnzGen")
    , nnzBtw("NnzBtw")
    , nnsGen("NnsEqtl")
    , nnsTot("NnsTot")
    , nnEqtl("NnEqtl")
    , nnEqtlOverlap("NnEqtlOverlap")
    , nnsPG("NnsPG")
    , nnGene("NnGene")
    , nsnp00("Nsnp00")
    , nsnp10("Nsnp10")
    , nsnp01("Nsnp01")
    , nsnp11("Nsnp11")
    , sigmaSqAlpha("varaEqtl")
    , geneticCorr("geneticCorRho")
    , sigmaSqBetaNonEqtl(varGenotypic,data.numIncdSnps,data.numNonEqtl,piEffNonEqtlVal,  noscale, "varbNEqtl")
    , sigmaSqBetaEqtl(varGenotypic,data.numIncdSnps,data.numEqtlOverlap,piEffEqtlVal,noscale, "varbEqtl")
    , sigmaSqTheta(h2snp,h2eQTL,data.numGenes,piThetaVal)
    , sigmaSqAlphaVec(data.geneEffectNames, data.varPhenotypiceQTL,h2eQTL,data.gene2cisSnpMap, piEffEqtlVal, noscale)
    , sigmaSqMats(data.geneEffectNames)
    , sigmaSqMatRes(data.snpEffectNames)
    // , genVarNonEqtlPrior(varGenotypic)
    // , noscale(noscale)
    // , scaleNonEqtlPrior(sigmaSqBetaNonEqtl.scale)
    , varg(varGenotypic, data.numKeptInds)
    , vargGene(data.geneEffectNames)
    , vargGeneCis(data.geneEffectNames)
    , cisHsq(data.geneEffectNames)
    , gwasCisHsq(data.geneEffectNames)
    , genicGwasEnrich(data.geneEffectNames,"genicGwasEnrich")
    , genicEqtlEnrich(data.geneEffectNames,"genicEqtlEnrich")
    , vare(data.ldblock2gwasSnpMap,data.varPhenotypic)
    , varEps(data.geneEffectNames,data.varPhenotypiceQTL,data.varResidualeQTL, data.numKeptIndseQTL,icrsq)
    , cisHsqMean("hsqCisMean")
    , vareMean("vareBlkMean")
    , varEpsMean("vareGenMean") {
        // single parameter
        paramVec = {&nnsGen, &nnzBtw, &nnsTot, &nnzGen,&nnEqtlOverlap, &nnsPG, &nnGene, &piEffEqtl, // aiao
            &nnEqtl,&nsnp00, &nsnp01, &nsnp10, &nsnp11, &piEffEieo1, &piEffEieo2, 
            &piEffNonEqtl, &sigmaSqBetaEqtl, &sigmaSqBetaNonEqtl, &sigmaSqAlpha, &geneticCorr,
            &hsq, &medHsq, &cisHsqMean, &vareMean, &varEpsMean, &varg};
        
        // vector parameter
        paramSetVec = {&snpEffects, &geneEffectVec, &cisHsq,&genicEqtlEnrich, &gwasCisHsq, &genicGwasEnrich};
        for(unsigned i = 0; i < eQTLJointVec.numGenes;i++){
            paramSetVec.push_back(eQTLJointVec[i]);
        }
        for(unsigned i = 0; i < snpEffectVec.numGenes;i++){
            paramSetVec.push_back(snpEffectVec[i]);
        }
        // matrix parameter
        for(unsigned i = 0; i < sigmaSqMats.numGenes;++i){
            paramMatVec.push_back(sigmaSqMats[i]);
            paramMatVec.push_back(sigmaSqMatRes[i]);
        }

        if(mcmcType == "AIAO"){
            LOGGER << "SBayesCO-AIAO model is used." << endl;
            paramToPrint = {&nnsTot, &nnzGen,&piEffNonEqtl,&piEffEqtl, &nnGene, &nnsPG, &sigmaSqBetaNonEqtl, &sigmaSqBetaEqtl, &sigmaSqAlpha, &hsq, &medHsq, &cisHsqMean, &vareMean, &varEpsMean};
            //  paramToPrint = {&piEffNonEqtl,&nnsTot, &sigmaSqBetaNonEqtl, &hsq,  &vareMean, &varg};
        }
        if (mcmcType == "EIEO"){
            LOGGER << "SBayesCO-EIEO model is used." << endl;
            for(unsigned i = 0; i < eQTLJointVec.numGenes;i++){
                paramSetVec.push_back(eQTLJointVecLatent[i]);
                paramSetVec.push_back(snpEffectVecLatent[i]);
            }
            paramToPrint = {&nnsTot,&nnzGen, &nnEqtl,&piEffEieo1, &piEffNonEqtl, &nnGene, &nnsPG, &nsnp00, &nsnp10, &nsnp01, &nsnp11, &sigmaSqBetaNonEqtl, &sigmaSqBetaEqtl, &sigmaSqAlpha, &hsq, &medHsq, &cisHsqMean, &vareMean, &varEpsMean};
        }
        if (modelPS) {
            paramVec.push_back(&ps);
            paramToPrint.push_back(&ps);
        }
        if (diagnose) {
            // nro.out.open((data.label+".diagnostics").c_str());
            // paramVec.push_back(&nro);
            // paramToPrint.push_back(&nro);
        }
        if (spouseCorrelation) {
            paramVec.push_back(&covg);
            paramToPrint.push_back(&covg);
        }
        if (message) {
            LOGGER << "\nApproximate BayesCO model fitted." << endl;
            
            if (robustMode) LOGGER << "Using a more robust parameterisation " << endl;
        }
       // if (randomStart) sampleStartVal();
        if (true) setStartVal();

        if(true){
            string outPath = data.label;
            // LOGGER << " vary: " << data.varPhenotypic << " vare: " << vare.valueVec << endl;; 
            if(data.numNonEqtl) {
                LOGGER << "sigmaSqBetaNonEqtlVec " << sigmaSqBetaNonEqtl.value  << endl;
                LOGGER << "sigmaSqBetaNEQTL: " << sigmaSqBetaNonEqtl.value << endl;
                LOGGER << " pi for intergenic region: " << piEffNonEqtl.value << endl;
            }
            if(data.numKeptGenes){
                // LOGGER << " varPhenotypiceQTL: " << data.varPhenotypiceQTL << endl;
                LOGGER  << " sigmaSqBetaEqtl " << sigmaSqBetaEqtl.value << endl;
                LOGGER << "pi for genic regeion: " << piEffEqtl.value << endl;
                // LOGGER << "SigmaSqAlphaVec: " << sigmaSqAlphaVec.values << endl;
                // LOGGER << " varEps: " << varEps.values<< endl;
                string outFile;
                // gwas snplist
                std::ofstream file1, file2;
                outFile = "-prior-cpp.txt";
                file1.open((outPath + outFile).c_str());
                for(unsigned i = 0; i < data.numKeptGenes; i++){
                    GeneInfo * gene = data.keptGeneInfoVec[i];
                    file1 << gene->ensemblID << "\t" << sigmaSqAlphaVec.values[i] << "\t" << varEps.values[i] << "\t" << "SBCO-CPP" << endl;
                }
                file1.close();
            }
        }

    }
    void sampleUnknowns(void);
    void setStartVal(void);
};

#endif /* stratify_hpp */