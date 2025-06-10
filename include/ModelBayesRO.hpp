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

#ifndef MODEL_BAYESRO_HPP
#define MODEL_BAYESRO_HPP

#include "Model.hpp"
#include "ModelBayesCO.hpp"

class BayesRO : public BayesCO {
    // Prior for snp efect pi_1 * N(0, 0) + pi_2 * N(0, sig^2_beta * gamma_2) + pi_3 * N(0, sig^2_beta * gamma_3) + pi_3 * N(0, sig^2_beta * gamma_4)
//     // consider S as unknown to make inference on the relationship between MAF and effect size
public: 

    class SigmaSqBetaNonEqtl : public Parameter, public Stat::InvChiSq {
        // variance of snp effects has a scaled-inverse chi-square prior
    public:
        double df;  // hyperparameter nub;
        double scale;     // hyperparameter
        double vargPrior;
        int numNonEqtl;
        bool noscale;

        VectorXd gamma;
        unsigned ndists;
        VectorXd sigmaSqBetaNonEqtlVec;
        VectorXd sigmaSqNonEqtlInvVec;
        VectorXd sigmaSqBetaNonEqtlLogVec;

        SigmaSqBetaNonEqtl(VectorXd &gamma,const double vg, const unsigned numSnps, const unsigned numNonEqtl, 
        const VectorXd &piEffNonEqtlVec,const VectorXd snp2pq, const bool noscale, const string &lab = "varbNEqtl")
        :vargPrior(vg),numNonEqtl(numNonEqtl), gamma(gamma), ndists(int(gamma.size())) , Parameter(lab), df(4){
            if(numNonEqtl != 0){
                // value = vg / ( snp2pq.mean() * numNonEqtl * piEffNonEqtl);  // derived from prior knowledge on Vg and pi
                value = vargPrior*( (double) numNonEqtl/ (double) numSnps)/(snp2pq.mean()* (double) numNonEqtl * (piEffNonEqtlVec.array() * gamma.array()).sum());
                // LOGGER << "sigmaSqNE: " << value << endl;
                scale = 0.5f * value;
                // LOGGER << "snp2pq sum: " << snp2pq.sum() << " value: " << value << " scale: "  << scale << endl;
                // BayesR
                sigmaSqBetaNonEqtlVec.setZero(ndists);
                sigmaSqNonEqtlInvVec.setZero(ndists);
                sigmaSqBetaNonEqtlLogVec.setZero(ndists);
                for(unsigned k = 1; k < ndists; k++){
                    sigmaSqBetaNonEqtlVec[k] = gamma[k] * value;
                    sigmaSqNonEqtlInvVec[k] = 1.0 / sigmaSqBetaNonEqtlVec[k];
                    sigmaSqBetaNonEqtlLogVec[k] = log(sigmaSqBetaNonEqtlVec[k]);
                }
            }
        }
        void sampleFromFC(const double snpEffSumSq, const unsigned numSnpEff);
    };

    class SigmaSqBetaEqtl : public Parameter, public Stat::InvChiSq {
        // variance of snp effects has a scaled-inverse chi-square prior
    public:
        const double df;  // hyperparameter
        double scale;     // hyperparameter
        bool noscale;  // no scaling on the genotypes
        double vargPrior;
        int numEqtl;

        SigmaSqBetaEqtl(const VectorXd &gamma, const double vg, const unsigned numSnps, const unsigned numEqtl, const VectorXd piEffEqtlVec, 
        const VectorXd snp2pqZ, const bool noscale, const string &lab = "varbEqtl")
        :vargPrior(vg),numEqtl(numEqtl), Parameter(lab), df(4), noscale(noscale){
            if(numEqtl != 0){
                value = vg * ((double) numEqtl/ (double) numSnps)/(snp2pqZ.mean() * (double) numEqtl * (piEffEqtlVec.array() * gamma.array()).sum() );
                scale = 0.5f * value;  
                // cout << "sigmaSqBetaEql: " << value << " scale: " << scale << endl;  //// 1 
            }       
        }
    };

    class SigmaSqAlphaVec : public ParamSet,public Stat::InvChiSq, public vector<BayesC::VarEffects*> {
        public:
        VectorXd values;
        vector<string> geneNames;

        SigmaSqAlphaVec(const VectorXd &gamma, const vector<string> geneNames,const VectorXd varGenotypiceQTL,
                const map<int,vector<int>> &gene2cisSnpMap, const VectorXd &snp2pqZ, 
                const VectorXd piEffEqtlVec,const bool noscale):
            ParamSet("SigmaSqAlphaVec", geneNames)
            , geneNames(geneNames){
            values.setZero(geneNames.size());
            
            // cout << "varGenotypiceQTL: " << varGenotypiceQTL << endl;
           
            for(unsigned i = 0; i < geneNames.size(); ++i){
                vector<int> snpsInGene = gene2cisSnpMap.at(i);
                // cout << "snp2pqZ.mean(): " << snp2pqZ.mean() << endl;
                // cout << "piEffEqtlVec: " << piEffEqtlVec << endl;
                // cout << "gamma: " << gamma << endl;
                // cout << "(piEffEqtlVec.array() * gamma.array() ).sum(): " << (piEffEqtlVec.array() * gamma.array() ).sum() << endl;
                //  cout << "h: " << snp2pqZ.mean() * (double) snpsInGene.size() * (piEffEqtlVec.array() * gamma.array() ).sum() << endl;
                values[i] = varGenotypiceQTL[i] /(snp2pqZ.mean() * (double) snpsInGene.size() * (piEffEqtlVec.array() * gamma.array() ).sum() );
                // cout << " gene: " << geneNames[i] << " SigmaSqAlphaVec: " << values[i] << endl;
            }

        }
        void sampleFromPrior();
        void sampleFromFC(const MatrixXd & deltaMat, const MatrixXd &eQTLJointMat, const map<int,vector<int>> & gene2cisSnpMap);
    };

    class SigmaSqMatVec : public vector<ParamMat*>,public Stat::Gamma, public Stat::InvChiSq {
        // 
        public:
        vector<string> geneNames; // sigmaSqMat is symmatic matix 
        unsigned numGenes;
        VectorXd gamma;
        const unsigned ndists;
        vector< vector<MatrixXd> > sigmaSqMatVec;
        vector< vector<MatrixXd> > sigmaSqInvMatVec;
        vector< vector<MatrixXd> > varcovPriorMatVec;
        vector< VectorXd > sigmaSqDetLogVecVec;
        vector<unsigned> numNonZeroThetaVec;

        double sigmaSqBetaEqtlPM;
        double sigmaSqAlphaAll;
        // priors
        double scaleBetaEqtl;
        double scaleAlphaAll;
        double nub;
        double nua;
        VectorXd sigmaSqAlphaPM; // sigmaSq Alpha per xQTL
        VectorXd scaleAlphaPM; 

        SigmaSqMatVec(const VectorXd &gamma, const vector<string> &geneNames,const string &lab = "SigmaSqMatVec"):
        gamma(gamma),ndists(int(gamma.size())), geneNames(geneNames),numGenes(int(geneNames.size())){
            sigmaSqMatVec.resize(ndists);
            sigmaSqInvMatVec.resize(ndists);
            varcovPriorMatVec.resize(ndists);
            sigmaSqDetLogVecVec.resize(ndists);

            sigmaSqAlphaPM.resize(numGenes);
            scaleAlphaPM.resize(numGenes);

            for (unsigned k = 0; k <ndists; k++) {
                BayesCO::SigmaSqMat *oneGene =  new BayesCO::SigmaSqMat(geneNames);
                sigmaSqMatVec[k] = oneGene->sigmaSqMats;
                sigmaSqInvMatVec[k] = oneGene->sigmaSqInvMats;
                varcovPriorMatVec[k] = oneGene->varcovPriors;
                sigmaSqDetLogVecVec[k] = oneGene->sigmaSqDetLogVec;
                
                delete oneGene;
            } // end of ndists
            
            scaleBetaEqtl = 0.0;
            scaleAlphaAll = 0.0;
            nub = 4.0;
            nua = 4.0;
            sigmaSqBetaEqtlPM = 0;
            sigmaSqAlphaAll = 0;
            sigmaSqAlphaPM.setZero();
            scaleAlphaPM.setZero();
        }  // end of constructor

        void sampleFromFCInvWishartGeneral(VectorXd &gamma, const double &ssqBetaEqtl, const double &ssqAlphaEqtl,vector<Matrix2d> ssqEqtlMat, 
                            const VectorXd &numNonZerosEqtlVec,const vector<unsigned> numNonZerosEqtlVecAcrossGenesPostIW,const bool messageBool);
        void sampleFromFCIWCorr(VectorXd &gamma, const double &ssqBetaEqtl, const double &ssqAlphaEqtl,vector<Matrix2d> ssqEqtlMat, 
                            const VectorXd &numNonZerosEqtlVec,const vector<unsigned> numNonZerosEqtlVecAcrossGenesPostIW,const bool messageBool);

        void setPrior(const VectorXd &gamma, const double &sigmaSqBetaEqtl, const VectorXd &sigmaSqAlphaVec);

    }; // end of class SigmaSqMatVec
    
    class SnpEffects : public BayesCO::SnpEffects {
    public:
        vector<vector<unsigned> > snpset;
        VectorXd gamma;
        unsigned ndists;
        double sum2pq;

        VectorXd nsnpDistNonEqtl;
        VectorXd nsnpDistEqtl;

        double ssqNonEqtl;
        double ssqBetaEqtl;
        double ssqAlphaEqtl;
        vector<Matrix2d> ssqEqtlMat; // need to initialize
        vector<unsigned> numNonZerosEqtlVecAcrossGenesPostIW;

        SnpEffects(VectorXd &gamma, const vector<string> &header): BayesCO::SnpEffects(header),
        gamma(gamma), ndists(int(gamma.size())) {
            sum2pq = 0.0;
            ssqNonEqtl = 0;
            ssqBetaEqtl = 0;
            ssqAlphaEqtl = 0;

            nsnpDistNonEqtl.setZero(ndists);
            nsnpDistEqtl.setZero(ndists);
        }
        void sampleFromFCAIAO(const VectorXd &gamma, VectorXd &wbcorr, vector<VectorXd> &wAcorr, const vector<MatrixDat> &Z, const vector<MatrixDat> ZGene, const map<string, vector<int> > &genePheIdxMap, 
            SnpEffectVec &snpEffectVec, EQTLJointVec &eQTLJointVec, 
                MatrixXd &deltaMat, const map<int,string> &gwasSnpIdx2snpIDMap, map<string, int> geneID2IdxMap,const map<string,vector<string> > &gwasSnpID2geneIDMap, const map<string, int> cisSnpID2IdxMap,
                SigmaSqBetaNonEqtl &sigmaSqBetaNonEqtl, SigmaSqMatVec &sigmaSqMatVec,const VectorXd &piEffEqtlVec, const VectorXd &piEffNonEqtlVec, const VectorXd varEps, const double &vare);
    };

    class ProbMixComps : public vector<Parameter*>, public Stat::Dirichlet {
        // prior probability of a snp being in any of the distributions effect has a dirichlet prior
        // 
    public:
        VectorXd alphaVec;  // hyperparameter
        VectorXd values;
        const unsigned ndist;

        ProbMixComps(const VectorXd &pis, const VectorXd &alphas): ndist(pis.size()){
            for (unsigned i = 0; i<ndist; ++i) {
                 //Parameter * pi = new Parameter("Pi");
                 this->push_back(new Parameter("Pi" + to_string(static_cast<long long>(i + 1))));
            }
            if (alphas.size() != ndist) alphaVec.setOnes(ndist);
            else alphaVec = alphas;
            values = pis;
        }
        
        void sampleFromFC(const VectorXd &snpStore);
    };
    
    class VgMixComps : public vector<Parameter*> {
    public:
        VectorXd values;
        const unsigned ndist;
        unsigned zeroIdx, minIdx;
        
        VgMixComps(const VectorXd &gamma): ndist(gamma.size()){
            double min = 1.0;
            minIdx = 0;
            for (unsigned i = 0; i<ndist; ++i) {
                this->push_back(new Parameter("Vg" + to_string(static_cast<long long>(i + 1))));
                if (gamma[i] == 0) zeroIdx = i;
                else if (gamma[i] < min) {
                    min = gamma[i];
                    minIdx = i;
                }
            }
            values.setZero(ndist);
        }
        
        void compute(const VectorXd &snpEffects, const MatrixXd &Z, const vector<vector<unsigned> > snpset, const double varg);
    };

    class NumSnpMixComps : public vector<Parameter*> {
    public:
        VectorXd values;
        const unsigned ndist;
        
        NumSnpMixComps(const VectorXd &pis): ndist(pis.size()){
            for (unsigned i = 0; i<ndist; ++i) {
                this->push_back(new Parameter("NumSnp" + to_string(static_cast<long long>(i + 1))));
            }
            values.setZero(ndist);
            // cout << " pis: " << pis << endl;
        }
        void getValues(const VectorXd &snpStore);
    };

    class Gammas : public ParamSet {
        // Set of scaling factors for each of the distributions
    public:
        Gammas(const VectorXd &gamma, const vector<string> &header, const string &lab = "gamma"): ParamSet(lab, header){
            values = gamma;
        }
    };

public:
    // VectorXd snpStore;   
    SnpEffects snpEffects;
    // VarEffects sigmaSq;
    // ProbMixComps Pis;
    ProbMixComps piEffEqtlVec;   // extend piEffEqtl to vector
    ProbMixComps piEffNonEqtlVec;
    VgMixComps VgEqtlVec;
    VgMixComps VgNonEqtlVec;
    NumSnpMixComps numSnpsEqtlVec;
    NumSnpMixComps numSnpsNonEqtlVec;
    Gammas gamma;
    SigmaSqMatVec sigmaSqMatVec;
    
    SigmaSqBetaNonEqtl sigmaSqBetaNonEqtl;
    SigmaSqBetaEqtl  sigmaSqBetaEqtl;

    
    SigmaSqAlphaVec sigmaSqAlphaVec;   // vector vlaues for various genes

    
    BayesRO( const Data &data,const string mcmcType,const bool eieoLatent, const double varGenotypic, const double varResidual, const double varRandom ,
                const double h2snp, const double h2eQTL,const double pival,const double piEffEqtlVal, const double piEffNonEqtlVal,const VectorXd &piEffEqtlVecVal,
                const VectorXd &piEffNonEqtlVecVal,const double piThetaVal, const double piAlpha, const double piBeta, const bool estimatePi, const bool noscale,
                const double phi, const double overdispersion, const bool estimatePS, const double icrsq, const double spouseCorrelation,
                const bool diagnosticMode, const bool robustMode, const VectorXd pis, const VectorXd &piPar, 
                VectorXd gamma, const bool originalModel, const string &algorithm, const bool message = true):
        BayesCO(data,mcmcType,eieoLatent, varGenotypic, varResidual, varRandom, h2snp, h2eQTL,  pival, piEffEqtlVal, piThetaVal, piEffNonEqtlVal, 
                piAlpha, piBeta, estimatePi, noscale, phi, overdispersion,  estimatePS, icrsq, spouseCorrelation, diagnosticMode, robustMode,
                true, true)
        //, Pis(pis, piPar),
        , piEffEqtlVec(piEffEqtlVecVal, piPar)
        , piEffNonEqtlVec(piEffNonEqtlVecVal, piPar)
        , numSnpsEqtlVec(piEffEqtlVecVal)
        , numSnpsNonEqtlVec(piEffNonEqtlVecVal)
        , VgEqtlVec(gamma)
        , VgNonEqtlVec(gamma)
        , gamma(gamma, vector<string>(gamma.size()))
        , snpEffects(gamma, data.snpEffectNames)
        , sigmaSqMatVec(gamma,data.geneEffectNames)
        , sigmaSqBetaNonEqtl(gamma,varGenotypic,data.numIncdSnps,data.numNonEqtl,piEffNonEqtlVecVal, data.snp2pq, noscale, "varbNEqtl")
        , sigmaSqBetaEqtl(gamma, varGenotypic,data.numIncdSnps,data.numEqtl,piEffEqtlVecVal,data.snp2pqEqtl,noscale, "varbEqtl")
        , sigmaSqAlphaVec(gamma, data.geneEffectNames, data.varGenotypiceQTL,data.gene2cisSnpMap, data.snp2pqEqtl, piEffEqtlVecVal, noscale)
        {

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
        
            if (message) {

                if (noscale)
                {
                    LOGGER << "Fitting model assuming unscaled genotypes." << endl;
                } else
                {
                    LOGGER << "Fitting model assuming scaled genotypes."  << endl;
                }
                if (robustMode) LOGGER << "Using a more robust parameterisation." << endl;
            }
            if (true) setStartVal( data);

        }  
    void setStartVal(Data data); 
    void sampleUnknowns(void);
};
   

#endif
