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

#ifndef MODEL_BAYESCO_HPP
#define MODEL_BAYESCO_HPP

#include "Model.hpp"

class BayesCO : public BayesC {
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
    class DeltaMat : public vector<ParamMat*> {
        public:
        vector<string> geneNames;
        vector<string> snpNames;
        unsigned numGenes;
        unsigned numSnps;
        unsigned numTraits;
        DeltaMat(const vector<string> &snpNames,
        const vector<string> &geneNames,
        const unsigned &numTraits,
        const string &lab = "DeltaMat"):
        numTraits(numTraits), 
        numSnps(int(snpNames.size())),
        numGenes(int(geneNames.size())){   
            for(unsigned i = 0; i < numTraits; ++i){
                this->push_back(new ParamMat(lab,snpNames,geneNames));
                // values[i] = (*this)[i]->values;
            }  
        }
    };
    
    class Intercept : public Parameter, public Stat::Flat {
        // all fixed effects has flat prior
    public:
        Intercept(const VectorXd &y,const string &lab = "Intercept")
        : Parameter(lab){
            value = Gadget::calcMean(y);
        }
        void sampleFromFC(VectorXd &wbcorr, const double &vare, const VectorXd &nGWAS);
    };
    class InterceptEQTL : public ParamSet, public Stat::Flat {
        // all fixed effects has flat prior
    public:
        unsigned numGenes;
        InterceptEQTL(const vector <VectorXd > &genePheVec,const vector<string> &header, const string &lab = "InterceptEQTL")
        :numGenes(header.size()) ,ParamSet(lab, header){
            values = Gadget::calColMeans(genePheVec);
        }
        
        void sampleFromFC(vector<VectorXd> &wAcorr, const map<string, vector<int> > &genePheIdxMap,
        const VectorXd &varEps, const vector<VectorDat> &neQTL);
    };

    class DeltaVec : public vector<ParamSet*> {
        public:
        vector<string> colnames;
        vector<string> geneNames;
        unsigned numGenes;
        map<int, vector<string>> gene2cisSnpIDMap;
        DeltaVec(const vector<string> &geneNames,
        const  map<int, vector<string>> &gene2cisSnpIDMap,
        const string &lab = "DeltaVec"):
        numGenes(int(geneNames.size())),geneNames(geneNames),gene2cisSnpIDMap(gene2cisSnpIDMap){   
            colnames.resize(numGenes);
            for(unsigned i = 0; i< numGenes; i++){
                colnames[i] = lab + geneNames[i];
                 this->push_back(new ParamSet(colnames[i], gene2cisSnpIDMap.at(i)));
            }
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

        SigmaSqBetaNonEqtl(const double vg, const unsigned numSnps, const unsigned numNonEqtl, 
        const double piEffNonEqtl,const VectorXd snp2pq, const bool noscale, const string &lab = "varbNEqtl")
        :vargPrior(vg),numNonEqtl(numNonEqtl), Parameter(lab), df(4){
            if(numNonEqtl != 0){
                // value = vg / ( snp2pq.mean() * numNonEqtl * piEffNonEqtl);  // derived from prior knowledge on Vg and pi
                value = vargPrior*( (double) numNonEqtl/ (double) numSnps)/(snp2pq.mean()* (double) numNonEqtl * piEffNonEqtl);
                // LOGGER << "sigmaSqNE: " << value << endl;
                scale = 0.5f * value;
                LOGGER << "snp2pq sum: " << snp2pq.sum() << " value: " << value << " scale: "  << scale << endl;
            }
        }
        void sampleFromFC(const double snpEffSumSq, const unsigned numSnpEff);
        void sampleFromPrior(const double piEffEqtl);
    };

    class SigmaSqBetaEqtl : public Parameter, public Stat::InvChiSq {
        // variance of snp effects has a scaled-inverse chi-square prior
    public:
        const double df;  // hyperparameter
        double scale;     // hyperparameter
        bool noscale;  // no scaling on the genotypes
        double vargPrior;
        int numEqtl;

        SigmaSqBetaEqtl(const double vg, const unsigned numSnps, const unsigned numEqtl, const double piEffEqtl, 
        const VectorXd snp2pqZ, const bool noscale, const string &lab = "varbEqtl")
        :vargPrior(vg),numEqtl(numEqtl), Parameter(lab), df(4), noscale(noscale){
            // if (noscale) {
            //     value = vg / (numEqtl * piEffEqtl);  // derived from prior knowledge on Vg and pi
            // } else {
            //     value = vg / (numEqtl * piEffEqtl);  // derived from prior knowledge on Vg and pi
            // }
            if(numEqtl != 0){
                value = vg * ((double) numEqtl/ (double) numSnps)/(snp2pqZ.mean() * (double) numEqtl * piEffEqtl);
                // cout << "snp2pqZ mean: " << snp2pqZ.mean() << endl;
                scale = 0.5f * value;  
            }            
        }
        // void sampleFromFC(const double snpEffSumSq, const unsigned numSnpEff);
        // void sampleFromPrior(const double piEffEqtl);
    };

    class SigmaSqAlphaVec : public ParamSet,public Stat::InvChiSq, public vector<BayesC::VarEffects*> {
        public:
        VectorXd values;
        vector<string> geneNames;

        SigmaSqAlphaVec(const vector<string> geneNames,const VectorXd varGenotypiceQTL,
                const map<int,vector<int>> &gene2cisSnpMap, const VectorXd &snp2pqZ, 
                const double piEffEqtl,const bool noscale):
            ParamSet("SigmaSqAlphaVec", geneNames)
            , geneNames(geneNames){
            values.setZero(geneNames.size());
            for(unsigned i = 0; i < geneNames.size(); ++i){
                vector<int> snpsInGene = gene2cisSnpMap.at(i);
                values[i] = varGenotypiceQTL[i] /(snp2pqZ.mean() * (double) snpsInGene.size() * piEffEqtl);
            }
        }
        void sampleFromPrior();
        void sampleFromFC(const MatrixXd & deltaMat, const MatrixXd &eQTLJointMat, const map<int,vector<int>> & gene2cisSnpMap);
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
            // LOGGER << "SigmaSqTheta: " << value << " scale: " << scale << " piTheta: " << piTheta << " h2snp: " << h2snp;
            // LOGGER << " h2eQTL " << h2eQTL << endl;
        }
        void sampleFromFC(const double snpEffSumSq, const unsigned numSnpEff);
    };

    class SigmaSqMat : public vector<ParamMat*>, public Stat::InvChiSq {
        public:
        vector<string> geneNames; // sigmaSqMat is symmatic matix 
        unsigned numGenes;
        vector<MatrixXd> sigmaSqMats;
        vector<MatrixXd> sigmaSqInvMats;
        vector<MatrixXd> varcovPriors;
        VectorXd sigmaSqDetLogVec;
        unsigned numNonZeroTheta;
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
            sigmaSqDetLogVec.setZero(numGenes);
            for(unsigned i = 0; i < numGenes; ++i){
                vector<string> sigmaSqName = {"nonEqtl",geneNames[i]};
                this->push_back(new ParamMat(geneNames[i],sigmaSqName,sigmaSqName));
                sigmaSqMats[i].setZero(2,2);
                sigmaSqInvMats[i].setZero(2,2);
                varcovPriors[i].setZero(2,2);
                // sigmaSqDetLogVec.setZero();
            }
            scaleBetaEqtl = 0.0;
            scaleAlphaAll = 0.0;
            nub = 4.0;
            nua = 4.0;
            sigmaSqBetaEqtlPM = 0.0;
            sigmaSqAlphaAll = 0.0;
        }
        void sampleFromFCInvWishartGeneral(const double &ssqBetaEqtl, const double &ssqAlphaEqtl, vector<Matrix2d> ssqEqtlMat,
         const VectorXd &numNonZerosEqtlVec, const vector<unsigned> numNonZerosEqtlVecAcrossGenesPostIW,const bool messageBool);

        void setPrior(const double &sigmaSqBetaEqtl, const VectorXd &sigmaSqAlphaVec);
       // VectorXd compute
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
            GeneEffects(const vector<string> geneNames):
            ParamSet("GeneEffects",geneNames),
            numGenes(int(geneNames.size())),
            geneNames(geneNames){
                propMed = 0.0;
                nnGene = 0;
                deltaGene.setZero(numGenes);
                vareMed = 0.01;
            }

        void sampleFromeFC(Data data, int iter ,VectorXd &betaTotal, MatrixXd &eQTLJointMat,
                            const map<int,vector<int> > &gene2gwasSnpMap,const map<int,vector<int>> &gene2cisSnpMap,
                            double &sigmaSqTheta, double &vareMed, const double piTheta, VectorXd deltaGene);
    };

   //////// heritability related parameter
    class MediatedHeritability : public Parameter {
        // compute heritability based on sampled values of genotypic and residual variances
        // strictly speaking, this is not a model parameter
        public:
        MediatedHeritability(const string &lab = "hsqMed"): Parameter(lab){};
        void compute(const double &geneVar, const double varg, const double vare){
            value = geneVar/ (varg + vare);
            // LOGGER << "medHsq: " << value << endl;

        }
    };
    class Heritability : public Parameter {
        // compute heritability based on sampled values of genotypic and residual variances
        // strictly speaking, this is not a model parameter
    public:
        Heritability(const string &lab = "hsq"): Parameter(lab){};
        void compute(const double genVar, const double vare){
            value = genVar/(genVar + vare);
        }
    };

    class HeritabilityCis : public ParamSet {
    public:
        // double cisHsqMean;
        HeritabilityCis(const vector<string> &header, const string &lab = "cisHsq"): ParamSet(lab,header){};
        void compute(const VectorXd geneVar, const VectorXd varyQTl){
            values = geneVar.array() / varyQTl.array();
        }
    };

    class GenotypicVar : public ApproxBayesC::GenotypicVar {
        public:
        GenotypicVar(const double varg, const unsigned n): ApproxBayesC::GenotypicVar(varg,n){}
        void compute(VectorXd &betaTotal,const MatrixXd &Z);
        void compute(VectorXd &ghat){
            value = Gadget::calcVariance(ghat);
        }
    };

    class GenotypicVarGene : public Parameter {
        public:
        GenotypicVarGene(const string &lab = "vargGene"): Parameter(lab){};
        void compute(VectorXd &geneEffects,  MatrixXd &eQTLJointMat,const MatrixXd &ZGene, const map<string, vector<int> > &genePheIdxMap, const vector<string> &geneEffectNames, const map<int,vector<int>> &gene2cisSnpMap);
        void compute(const int &nGWAS,const map<int,VectorXd> &gwhatGwasMap, const VectorXd geneEffectVec){
            int numGenes  = gwhatGwasMap.size();
            VectorXd gexp = VectorXd::Zero(nGWAS);
            for(unsigned i = 0; i < numGenes; i++ ){
                gexp += gwhatGwasMap.at(i) * geneEffectVec(i);
            }

            value = Gadget::calcVariance(gexp);
        }

    };

    class GenotypicVarGeneCis : public ParamSet {
        public:
        GenotypicVarGeneCis(const vector<string> &header, const string &lab = "vargGeneCis"): ParamSet(lab,header){};
        void compute( MatrixXd &eQTLJointMat, const MatrixXd &ZGene, const map<string, vector<int> > &genePheIdxMap,
                    const vector<string> &geneEffectNames, const map<int,vector<int>> &gene2cisSnpMap);
        void compute(map<int,VectorXd> &gwhatMap){
            int numGenes  = gwhatMap.size();
            for(unsigned i = 0; i < numGenes; i++ ){
                values[i] = Gadget::calcVariance(gwhatMap.at(i));
            }
        }
    };
    
    class ResidualVar : public BayesC::ResidualVar {
        public:     
        ResidualVar(const double &vare, const unsigned &nobs, const double &icrsq): BayesC::ResidualVar(vare,nobs) {
        }
        void sampleFromFC( const VectorXd &ycorr, const VectorXd &nGWAS);
        void sampleFromFC( const VectorXd &y,const MatrixXd X, const VectorXd &snpEffect, const VectorXd &nGWAS);
    };


    class ResidualVareEQTL : public ParamSet, public vector<BayesC::ResidualVar*>,public Stat::InvChiSq {
        public:
        // VectorXd values;
        vector<string> geneNames;
        
        ResidualVareEQTL(const vector<string> geneNames, const VectorXd vare, const unsigned nobs, const double icrsq):
            ParamSet("ResidualVareEQTL",geneNames), geneNames(geneNames){
            for(unsigned i = 0; i < geneNames.size(); ++i){
                this->push_back(new BayesC::ResidualVar(vare[i],nobs));
                values[i] = (*this)[i]->value;
                // values[i] = 2.207168e-02;
            }            
        }
        void sampleFromFC(const vector<VectorXd> &wAcorr, const vector<VectorDat> &neQTL);
        void sampleFromFC( const vector<VectorXd> &w, const MatrixXd ZGene, EQTLJointVec &eQTLJointVec, const vector<VectorDat> &neQTL,
                                            const map<string, vector<int> > &genePheIdxMap, const vector<string> &geneEffectNames,
                                            const map<int,vector<int>> &gene2cisSnpMap);
    };

    class ResidualMat : public vector<ParamMat*>,public Stat::Gamma, public Stat::InvChiSq {
        public:
        vector<string> geneNames; // sigmaSqMat is symmatic matix 
        unsigned numGenes;
        vector<MatrixXd> ResMats;
        vector<MatrixXd> ResInvMats;
        vector<MatrixXd> varcovPriors;

        double sigmaSqBetaEqtlPM;
        double sigmaSqAlphaAll;

        ResidualMat(const vector<string> &geneNames,const string &lab = "SigmaSqMats"):
        numGenes(int(geneNames.size()))
        , geneNames(geneNames){
            ResMats.resize(numGenes);
            ResInvMats.resize(numGenes);
            varcovPriors.resize(numGenes);
            for(unsigned i = 0; i < numGenes; ++i){
                vector<string> sigmaSqName = {"nonEqtl",geneNames[i]};
                this->push_back(new ParamMat(geneNames[i],sigmaSqName,sigmaSqName));
                ResMats[i].setZero(2,2);
                ResInvMats[i].setZero(2,2);
                varcovPriors[i].setZero(2,2);
                // sigmaSqDetLogVec.setZero();
            }
            sigmaSqBetaEqtlPM = 0;
            sigmaSqAlphaAll = 0;
        }

        void setPrior(const double &sigmaSqBetaEqtl, const VectorXd &sigmaSqAlphaVec);
       // VectorXd compute
    };

    class SnpEffects : public BayesC::SnpEffects, public Stat::MultiNormal {
    // for the ease of sampling, we model the SNP effect to be alpha_j = beta_j * delta_j where beta_j has a univariate normal prior.
    public:
        VectorXd betaTotal;     // save sample squres of full conditional normal distribution regardless of delta values
        VectorXd snpEffectVec;
        vector<int> numSnpComp; // eieo 
        Vector2d numNonZerosEqtlVec;    // number of nonzero effects in genic regions for beta and alpha Pair
        unsigned numNonZerosNonEqtl;    // number of nonzero effects in between gene regions
        unsigned numNonNullEqtl;  // number of SNPs with at least one nonzero effect in genic regions
        unsigned numNonNullSnpTot;  // total number of SNPs with nonzero effects in the genome
        double numNonNullSnpPerGene;
        double nnG;                  // number of genes with at least one null zero eQTLs.
        VectorXd betaTotalMean; // debug
        // EIEO
        Vector4d numSnpCompVec; // store for nsnp00-11 
        // calculate heritability directly
        VectorXd ghat;  // used to calculate total snp heritability
        map<int, VectorXd> gwhatMap; // used to calculate 
        map<int, VectorXd> gwhatGwasMap; // used to calculate mediated heritability
        // snp effect
        VectorXd betaTotalLatent;
        vector<unsigned> numNonZerosEqtlVecAcrossGenesPostIW;
        vector<Matrix2d> ssqEqtlMat;
        double ssqNonEqtl;
        double ssqBetaEqtl;
        double ssqAlphaEqtl;
        double ssqBetaTotalGenic; // ssq for betaTotal in the genic region
        VectorXd ssqAlphaEqtlPG; // ssq for alpha per gene
        VectorXd ssqBetaEqtlPG;
        VectorXd numNonZerosGenicGwasVec; // used to cauculate genic gwas enrichment.
        VectorXd numNonZerosGenicEqtlVec; // used to cauculate genic gwas enrichment.
        unsigned numNonNullBetaTotGenic; // number of SNPs with at least one nonzero gwas effect



        SnpEffects(const vector<string> &header): BayesC::SnpEffects(header, "Gibbs"){
            betaTotal.setZero(size);
            numNonZerosEqtlVecAcrossGenesPostIW.reserve(2);
            numSnpComp.reserve(4);
            betaTotalMean.setZero(size);
            numNonZerosEqtlVec.setZero(2);
        }
        // AIAO model
        void sampleFromFCAIAO(VectorXd &wbcorr, vector<VectorXd> &wAcorr,const vector<MatrixDat> &Z,const vector<MatrixDat> ZGene, const map<string, vector<int> > &genePheIdxMap, 
            SnpEffectVec &snpEffectVec, EQTLJointVec &eQTLJointVec, MatrixXd &deltaMat, const map<int,string> &gwasSnpIdx2snpIDMap, map<string, int> geneID2IdxMap,
            const map<string,vector<string> > &gwasSnpID2geneIDMap, const map<string, int> cisSnpID2IdxMap, const double sigmaSqBetaNonEqtl,
            SigmaSqMat &sigmaSqMats, const double &piEffEqtl, const double &piEffNonEqtl,const VectorXd varEps, const double &vare);

        // EIEO model
        void sampleFromFCEIEO(VectorXd &wbcorr, vector<VectorXd> &wAcorr,const vector<MatrixDat> &Z,const vector<MatrixDat> ZGene, const map<string, vector<int> > &genePheIdxMap, 
            SnpEffectVec &snpEffectVec, EQTLJointVec &eQTLJointVec, SnpEffectVec &snpEffectVecLatent, EQTLJointVec &eQTLJointVecLatent,
            MatrixXd &deltaMat, const map<int,string> &gwasSnpIdx2snpIDMap, map<string, int> geneID2IdxMap,
            DeltaVec deltaVecGWAS,DeltaVec deltaVecEQTL,
            const map<string,vector<string> > &gwasSnpID2geneIDMap, const map<string, int> cisSnpID2IdxMap, const double sigmaSqBetaNonEqtl,
            SigmaSqMat &sigmaSqMats, const double &piEffEieo1, const double &piEffEieo2, const double &piEffNonEqtl,const VectorXd varEps, const double &vare);
    };

public:
    string mcmcType;
    bool eieoLatent;
    vector<VectorXd> wAcorr;
    VectorXd wbcorr;
    
    MediatedHeritability medHsq;
    Heritability hsq;
    HeritabilityCis cisHsq;

    Intercept intercept; // intercept for gwas;
    InterceptEQTL interceptEqtl;

    // EQTLJointMat eQTLJointMat;
    // SnpEffectMat snpEffectMat;
    // EQTLJointMat eQTLJointMatLatent;
    // SnpEffectMat snpEffectMatLatent;
    
    SnpEffects snpEffects;
    EQTLJointVec eQTLJointVec;
    EQTLJointVec eQTLJointVecLatent;
    SnpEffectVec snpEffectVec;
    SnpEffectVec snpEffectVecLatent;

    // BayesC::VarEffects sigmaSqBetaNonEqtl;
    // BayesC::VarEffects sigmaSqBetaEqtl;
    SigmaSqBetaNonEqtl sigmaSqBetaNonEqtl;
    SigmaSqBetaEqtl  sigmaSqBetaEqtl;
    SigmaSqTheta sigmaSqTheta;
    
    SigmaSqAlphaVec sigmaSqAlphaVec;   // vector vlaues for various genes
    SigmaSqMat sigmaSqMats;  // matrix values for various genes
    ResidualMat residualMats; // residual matrix values for various genes;
    GeneEffects geneEffectVec;
    // EQTLJointVec eQTLJointVec;
    DeltaMat deltaMat;
    DeltaMat deltaTrait;
    ResidualVareEQTL varEps; // residuals for eQTL of genes
    ResidualVar vare;
    GenotypicVar varg;
    GenotypicVarGene vargGene;
    GenotypicVarGeneCis vargGeneCis;
    Rounding  rounding;
    BayesC::Pi piEffEqtl;   // pi for nonzero effects in genic regions (overlaping eQTL are allowed to have multiple nonzero effects)
    BayesC::Pi piEffEieo1; // used in EIEO model
    BayesC::Pi piEffEieo2; // used in EIEO model
    BayesC::Pi piTheta;
    BayesC::Pi piEffNonEqtl;   // pi for nonzero effects in between gene regions (thus the same of nonnull SNPs in intergenic regions)
    NumNonZeroSnp nnzBtw;  // number of nonzero effects in between gene regions
    NumNonZeroSnp nnzGen;  // number of nonzero effects in genic regions
    NumNonZeroSnp nnsGen;  // number of SNPs with at least one nonzero effect in genic regions
    NumNonZeroSnp nnsTot;  // total number of SNPs with nonzero effects in the genome
    NumNonZeroSnp nnsPG;   // average number of SNPs with nonzero effects per gene
    NumNonZeroSnp nnGene;     // nonzero gene effect;
    // NumNonZeroSnp nnzTheta;

    Parameter cisHsqMean;  // mean value of cis-heritability for gene expression
    Parameter vareMean;    // mean value of trait residual variance across LD blocks
    Parameter varEpsMean;  // mean value of expression residual variance across genes
    Parameter sigmaSqAlpha;
    // Parameter sigmaSqBetaEqtl;
    // Parameter sigmaSqBetaNonEqtl;

    BayesCO(const Data &data,const string mcmcType,const bool eieoLatent, const double varGenotypic, const double varResidual, const double varRandom ,
                const double h2snp, const double h2eQTL,const double pival,const double piEffEqtlVal, const double piThetaVal,
                const double piEffNonEqtlVal, const double piAlpha, const double piBeta, const bool estimatePi, const bool noscale,
                const double phi, const double overdispersion, const bool estimatePS, const double icrsq, const double spouseCorrelation,
                const bool diagnosticMode, const bool robustMode, const bool randomStart = true, const bool message = true)
    : BayesC(data, varGenotypic, varResidual, varRandom, pival, piAlpha, piBeta, estimatePi, noscale, "Gibbs", false)
    , wbcorr(data.wbcorr)
    , wAcorr(data.wAcorr)
    , intercept(data.y)
    , interceptEqtl(data.genePheVec,data.geneEffectNames)
    , snpEffects(data.snpEffectNames)
    , geneEffectVec(data.geneEffectNames)
    // , eQTLJointVec(data.snpEffectNames,data.snpEffectNames.size()* data.geneEffectNames.size())
    // , snpEffectMat(data.snpEffectNames,data.gwasAndGeneEffectNames)
    // , eQTLJointMat(data.cisSnpIDVec,data.geneEffectNames)
    // , snpEffectMatLatent(data.snpEffectNames,data.gwasAndGeneEffectNames)
    // , eQTLJointMatLatent(data.snpEffectNames,data.geneEffectNames)
    , deltaMat(data.snpEffectNames,data.gwasAndGeneEffectNames,1)
    , deltaTrait(data.snpEffectNames,data.gwasAndGeneEffectNames,2)
    , eQTLJointVec(data.geneEffectNames,data.gene2cisSnpIDMap,"EQTLJointVec_")
    , eQTLJointVecLatent(data.geneEffectNames,data.gene2cisSnpIDMap,"EQTLJointVecLatent_")
    , snpEffectVec(data.gwasAndGeneEffectNames,data.gwas2SnpIDMap,"SnpJointVec_")
    , snpEffectVecLatent(data.gwasAndGeneEffectNames,data.gwas2SnpIDMap,"SnpJointVecLatent_")
    , piEffEqtl(piEffEqtlVal, piAlpha, 1,"piEffEqtl")
    , piTheta(piThetaVal,piAlpha,piBeta,"piTheta")
    , piEffNonEqtl(piEffNonEqtlVal, piAlpha, piBeta,"piEffNonEqtl")
    , piEffEieo1(piEffEqtlVal, piAlpha, piBeta,"piEffEieo1")
    , piEffEieo2(piEffEqtlVal, piAlpha, piBeta,"piEffEieo2")
    , mcmcType(mcmcType)
    , eieoLatent(eieoLatent)
    , nnzGen("nnzGen")
    // , nnzTheta("nnzTheta")
    , nnzBtw("nnzBtw")
    , nnsGen("nnsGen")
    , nnsTot("nnsTot")
    , nnsPG("nnsPG")
    , nnGene("nnGene")
    , sigmaSqAlpha("varaEqtl")
    , sigmaSqBetaNonEqtl(varGenotypic,data.numIncdSnps,data.numNonEqtl,piEffNonEqtlVal, data.snp2pq, noscale, "varbNEqtl")
    , sigmaSqBetaEqtl(varGenotypic,data.numIncdSnps,data.numEqtl,piEffEqtlVal,data.snp2pqEqtl,noscale, "varbEqtl")
    // , sigmaSqBetaEqtl(varGenotypic,data.numEqtl,data.numEqtl,piEffEqtlVal,noscale, "varbEqtl")
    // , sigmaSqBetaNonEqtl(varGenotypic, data.snp2pq, piEffNonEqtlVal,noscale, "varbNEqtl")
    , sigmaSqTheta(h2snp,h2eQTL,data.numGenes,piThetaVal)
    , sigmaSqAlphaVec(data.geneEffectNames, data.varGenotypiceQTL,data.gene2cisSnpMap, data.snp2pqEqtl, piEffEqtlVal, noscale)
    , sigmaSqMats(data.geneEffectNames)
    , residualMats(data.geneEffectNames)
    , varg(varGenotypic, data.numKeptInds)
    , vargGeneCis(data.geneEffectNames)
    , cisHsq(data.geneEffectNames)
    , vare(data.varResidual,data.numKeptInds,icrsq)
    , varEps(data.geneEffectNames,data.varResidualeQTL, data.numKeptIndseQTL,icrsq)
    , cisHsqMean("hsqCisMean")
    , vareMean("vareBlkMean")
    , varEpsMean("vareGenMean") {
        paramVec = {&nnsGen, &nnzBtw, &nnsTot, &nnzGen, &nnsPG, &nnGene, &piEffEqtl,&piEffEieo1, &piEffEieo2, &piEffNonEqtl, &sigmaSqBetaEqtl, &sigmaSqBetaNonEqtl, &sigmaSqAlpha, &hsq, &medHsq, &cisHsqMean, &vareMean, &varEpsMean, &varg};
        paramSetVec = {&snpEffects, &geneEffectVec};
        for(unsigned i = 0; i < sigmaSqMats.numGenes;++i){
            paramMatVec.push_back(sigmaSqMats[i]);
        }
        for(unsigned i = 0; i < residualMats.numGenes;++i){
            paramMatVec.push_back(residualMats[i]);
        }

        for(unsigned i = 0; i < eQTLJointVec.numGenes;i++){
            paramSetVec.push_back(eQTLJointVec[i]);
        }
        for(unsigned i = 0; i < snpEffectVec.numGenes;i++){
            paramSetVec.push_back(snpEffectVec[i]);
        }

        if(mcmcType == "AIAO"){
            for(unsigned i = 0; i < deltaMat.numTraits; ++i){
                paramMatVec.push_back(deltaMat[i]);
            }
            paramToPrint = {&nnsTot, &nnsGen, &nnGene, &nnsPG, &sigmaSqBetaNonEqtl, &sigmaSqBetaEqtl, &sigmaSqAlpha, &hsq, &medHsq, &cisHsqMean, &vareMean, &varEpsMean};
            //  paramToPrint = {&piEffNonEqtl,&nnsTot, &sigmaSqBetaNonEqtl, &hsq,  &vareMean, &varg};
        }
        if (mcmcType == "EIEO"){
            for(unsigned i = 0; i < deltaTrait.numTraits; ++i){
                paramMatVec.push_back(deltaTrait[i]);
            }
            paramToPrint = {&nnsTot, &nnsGen, &nnGene, &nnsPG, &piEffEqtl, &piEffNonEqtl,&varg, &sigmaSqBetaNonEqtl, &sigmaSqBetaEqtl, &sigmaSqAlpha, &hsq, &medHsq, &cisHsqMean, &vareMean, &varEpsMean};
        }
        // if (modelPS) {
        //     paramVec.push_back(&ps);
        //     paramToPrint.push_back(&ps);
        // }
        // if (diagnose) {
        //     nro.out.open((data.label+".diagnostics").c_str());
        //     paramVec.push_back(&nro);
        //     paramToPrint.push_back(&nro);
        // }
        // if (spouseCorrelation) {
        //     paramVec.push_back(&covg);
        //     paramToPrint.push_back(&covg);
        // }
        if (message) {
            // LOGGER << "\nApproximate BayesCO model fitted." << endl;

            if (noscale)
            {
                LOGGER << "Fitting model assuming unscaled genotypes " << endl;
            } else
            {
                LOGGER << "Fitting model assuming scaled genotypes "  << endl;
            }
            if (robustMode) LOGGER << "Using a more robust parameterisation " << endl;
            // LOGGER << "vare: " << vare.value << " vary: " << data.varPhenotypic << " " ; // << varEps.values << endl;
            // // LOGGER << " gwas y: " << endl << data.y.head(2) << endl << " tail " << endl << data.y.tail(2) << endl;
            // // LOGGER << "gwas snp2pq sum " << data.snp2pq.sum() << " eqtl sum " << data.snp2pqEqtl.sum() << endl;
            // LOGGER << "sigmaSqBetaNonEqtl " << sigmaSqBetaNonEqtl.value << " sigmaSqBetaEqtl " << sigmaSqBetaEqtl.value << endl;
            // LOGGER << "piEffEqtl: " << piEffEqtl.value << " piEffNonEqtl: " << piEffNonEqtl.value << endl;
            // LOGGER << "SigmaSqAlphaVec: " << sigmaSqAlphaVec.values << endl;
            // LOGGER << " varEps: " << varEps.values << endl;
            // LOGGER << " varPhenotypiceQTL: " << data.varPhenotypiceQTL << endl;
            // LOGGER << "geneNames: " << endl;

            // for(unsigned i = 0; i < data.geneEffectNames.size(); i++){
            //     LOGGER << data.geneEffectNames[i] << endl;
            // }
        }
       // if (randomStart) sampleStartVal();
        if (true) setStartVal( data);

        //LOGGER << "piEffNonEqtlVal " << piEffNonEqtlVal << " varGenotypic " << varGenotypic << " sigmaSq " << sigmaSq.value << endl;

    }

    void sampleUnknowns(void);
    void setStartVal(Data data);
};


#endif