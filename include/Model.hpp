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


#ifndef model_hpp   
#define model_hpp

#include <iostream>
// #include <cmath>
#include "Stat.hpp"
#include "Data.hpp"

using namespace std;


class Parameter {
    // base class for a single parameter
public:
    const string label;
    double value;   // sampled value
    
    Parameter(const string &label): label(label){
        value = 0.0;
    }
};

class ParamSet {
    // base class for a set of parameters of same kind, e.g. fixed effects, snp effects ...
public:
    const string label;
    const vector<string> &header;
    std::map<string, int> name2index;
    unsigned size;
    VectorXd values;
        
    ParamSet(const string &label, const vector<string> &header)
    : label(label), header(header), size(int(header.size())){
        values.setZero(size);
        for (unsigned j = 0; j < size; j++) name2index.insert(pair<string, int>(header[j], j));
    }
    double getValue(string nameIdx) const {
        if(name2index.find(nameIdx) == name2index.end()){
            LOGGER.e(0,"Could not find snp "+ nameIdx + " in current paramset using getValue function.");
        }
         return values(name2index.at(nameIdx)); 
         } 
    void setValue(string nameIdx, double svalue){
        if(name2index.find(nameIdx) == name2index.end()){
            LOGGER.e(0,"Could not find snp "+ nameIdx + " in current paramset using setValue function.");
        }
        values(name2index.at(nameIdx)) = svalue;
    }
};

class ParamMat {
    // for matrix
public:
    const string label;
    const vector<string> &rownames;
    const vector<string> &colnames;
    unsigned nrow;
    unsigned ncol;
    MatrixXd values;
        
    ParamMat(const string &label, const vector<string> &rownames, 
    const vector<string> &colnames)
    : label(label),rownames(rownames),colnames(colnames),
    nrow(int(rownames.size())), ncol(int(colnames.size())) {
        //values.setZero(size);
        values.setZero(nrow,ncol);
    }
    void resetZero(){
        values.setZero(nrow,ncol);
    }
};

class Model {
public:
    unsigned numSnps;
        
    vector<ParamSet*> paramSetVec;
    vector<Parameter*> paramVec;
    vector<ParamMat*> paramMatVec;
    vector<Parameter*> paramToPrint;
    vector<ParamSet*> paramSetToPrint;
    
    virtual void sampleUnknowns(void) = 0;
    virtual void sampleStartVal(void) = 0;
};

class BayesC : public Model {
    // model settings and prior specifications in class constructors
public:
    
    class FixedEffects : public ParamSet, public Stat::Flat {
        // all fixed effects has flat prior
    public:
        FixedEffects(const vector<string> &header, const string &lab = "CovEffects")
        : ParamSet(lab, header){}
        
        void sampleFromFC(VectorXd &ycorr, const MatrixXd &X, const VectorXd &XPXdiag, const double vare);
    };
    
    class RandomEffects : public ParamSet, public Stat::Normal {
        // random covariate effects
    public:
        double ssq;  // sum of squares

        RandomEffects(const vector<string> &header, const string &lab = "RandCovEffects")
        : ParamSet(lab, header){}
        
        void sampleFromFC(VectorXd &ycorr, const MatrixXd &W, const VectorXd &WPWdiag, const VectorXd &Rsqrt, const bool weightedRes, const double sigmaSqRand, const double vare, VectorXd &rhat);
    };
    
    class VarRandomEffects : public Parameter, public Stat::InvChiSq {
        // variance of random covariate effects has a scaled-inverse chi-square prior
    public:
        const double df;  // hyperparameter
        double scale;     // hyperparameter

        VarRandomEffects(const double varRandom, const double numRandomEffects, const string &lab = "SigmaSqRand")
        : Parameter(lab), df(4)
        {
            //value = varRandom/numRandomEffects;
            value = varRandom;
            scale = 0.5f*value/numRandomEffects;  // due to df = 4
        }
        
        void sampleFromFC(const double randEffSumSq, const unsigned numRandEff);
    };

    
    class SnpEffects : public ParamSet, public Stat::NormalZeroMixture {
        // all snp effects has a mixture prior of a nomral distribution and a point mass at zero
    public:
        double sumSq;
        unsigned numNonZeros;
        
        VectorXd posteriorMean;
        VectorXd posteriorMeanPIP;
        VectorXd pip;
        
        enum {gibbs, hmc} algorithm;
        
        unsigned cnt;
        double mhr;

        bool shuffle;
        vector<int> snpIndexVec;

        SnpEffects(const vector<string> &header, const string &alg, const string &lab = "SnpEffects")
        : ParamSet(lab, header){
            sumSq = 0.0;
            numNonZeros = 0;
            posteriorMean.setZero(size);
            posteriorMeanPIP.setZero(size);
            pip.setZero(size);
            if (alg=="HMC") algorithm = hmc;
            else algorithm = gibbs;
            cnt = 0;
            mhr = 0.0;

            shuffle = true;  // shuffle SNP order
            snpIndexVec.resize(size);
            for (unsigned i=0; i<size; i++) snpIndexVec[i] = i;
        }
        
        void sampleFromFC(VectorXd &ycorr, const MatrixXd &Z, const VectorXd &ZPZdiag, const VectorXd &Rsqrt, const bool weightedRes,
                          const double sigmaSq, const double pi, const double vare, VectorXd &ghat);
        void gibbsSampler(VectorXd &ycorr, const MatrixXd &Z, const VectorXd &ZPZdiag, const VectorXd &Rsqrt, const bool weightedRes,
                          const double sigmaSq, const double pi, const double vare, VectorXd &ghat);
        void hmcSampler(VectorXd &ycorr, const MatrixXd &Z, const VectorXd &ZPZdiag,
                        const double sigmaSq, const double pi, const double vare, VectorXd &ghat);
        ArrayXd gradientU(const VectorXd &alpha, const MatrixXd &ZPZ, const VectorXd &ypZ,
                        const double sigmaSq, const double vare);
        double computeU(const VectorXd &alpha, const MatrixXd &ZPZ, const VectorXd &ypZ,
                       const double sigmaSq, const double vare);
        
        void sampleFromFC_omp(VectorXd &ycorr, const MatrixXd &Z, const VectorXd &ZPZdiag, const double sigmaSq, const double pi, const double vare, VectorXd &ghat);
        void computePosteriorMean(const unsigned iter);

    };
    
    class SnpPIP : public ParamSet {
    public:
        SnpPIP(const vector<string> &header, const string &lab = "PIP") : ParamSet(lab, header){}
        
        void getValues(const VectorXd &pip){values = pip;}
    };
    
    class VarEffects : public Parameter, public Stat::InvChiSq {
        // variance of snp effects has a scaled-inverse chi-square prior
    public:
        const double df;  // hyperparameter
        double scale;     // hyperparameter
        bool noscale;  // no scaling on the genotypes

        VarEffects(const double vg, const VectorXd &snp2pq, const double pi, const bool noscale, const string &lab = "SigmaSq")
        : Parameter(lab), df(4), noscale(noscale)
        {
            // LOGGER << "To scale or not to scale " << noscale << endl;
            // LOGGER << "Scale value 1 " << value << endl;
            if (noscale) {
                value = vg / (snp2pq.sum() * pi);  // derived from prior knowledge on Vg and pi
            } else {
                value = vg / (snp2pq.size() * pi);  // derived from prior knowledge on Vg and pi
            }
            
            scale = 0.5f*value;  // due to df = 4
            
            //LOGGER << value << " " << vg << " " << snp2pq.sum() << " " << pi << " " << noscale << endl;
        }
        
        void sampleFromFC(const double snpEffSumSq, const unsigned numSnpEff);
        void sampleFromPrior(void);
        void computeScale(const double varg, const VectorXd &snp2pq, const double pi);
        void computeScale(const double varg, const double sum2pq);
        void compute(const double snpEffSumSq, const double numSnpEff);

    };
    
    class ScaleVar : public Parameter, public Stat::Gamma {
        // scale factor of variance variable
    public:
        const double shape;
        const double scale;
        
        ScaleVar(const double val, const string &lab = "Scale"): shape(1.0), scale(1.0), Parameter(lab){
            value = val;  // starting value
        }
        
        void sampleFromFC(const double sigmaSq, const double df, double &scaleVar);
        void getValue(const double val){ value = val; };
    };
    
    class Pi : public Parameter, public Stat::Beta {
        // prior probability of a snp with a non-zero effect has a beta prior
    public:
        const double alpha;  // hyperparameter
        const double beta;   // hyperparameter
        
        Pi(const double pi, const double alpha, const double beta, const string &lab = "Pi"): Parameter(lab), alpha(alpha), beta(beta){  // informative prior
            value = pi;
        }
        
        void sampleFromFC(const unsigned numSnps, const unsigned numSnpEff);
        void sampleFromPrior(void);
        void compute(const double numSnps, const double numSnpEff);
    };
    
    
    class ResidualVar : public Parameter, public Stat::InvChiSq {
        // residual variance has a scaled-inverse chi-square prior
    public:
        const double df;      // hyperparameter
        const double scale;   // hyperparameter
        unsigned nobs;
        
        ResidualVar(const double vare, const unsigned n, const string &lab = "ResVar")
        : Parameter(lab), df(4)
        , scale(0.5f*vare){
            nobs = n;
            value = vare;  // due to df = 4
        }
        
        void sampleFromFC(VectorXd &ycorr);
    };
    
    class GenotypicVar : public Parameter {
        // compute genotypic variance from the sampled SNP effects
        // strictly speaking, this is not a model parameter
    public:
        GenotypicVar(const double varg, const string &lab = "GenVar"): Parameter(lab){
            value = varg;
        };
        void compute(const VectorXd &ghat);
    };
    
    class RandomVar : public Parameter {
        // compute variance explained due to random covariate effects
    public:
        RandomVar(const double varRandom, const string &lab = "RanVar"): Parameter (lab){
            value = varRandom;
        }
        
        void compute(const VectorXd &rhat);
    };
    
    class Heritability : public Parameter {
        // compute heritability based on sampled values of genotypic and residual variances
        // strictly speaking, this is not a model parameter
    public:
        Heritability(const string &lab = "hsq"): Parameter(lab){};
        void compute(const double genVar, const double resVar){
            value = genVar/(genVar+resVar);
        }
    };
    
    class Rounding : public Parameter {
        // re-compute ycorr to eliminate rounding errors
    public:
        unsigned count;
        
        Rounding(const string &lab = "Rounding"): Parameter(lab){
            count = 0;
        }
        void computeYcorr(const VectorXd &y, const MatrixXd &X, const MatrixXd &W, const MatrixXd &Z,
                          const VectorXd &fixedEffects, const VectorXd &randomEffects, const VectorXd &snpEffects,
                          VectorXd &ycorr);
    };
    
    class NumNonZeroSnp : public Parameter {
        // number of non-zero SNP effects
    public:
        NumNonZeroSnp(const string &lab = "NnzSnp"): Parameter(lab){};
        void getValue(const unsigned nnz){ value = nnz; };
    };

    class varEffectScaled : public Parameter {
        // Alternative way to estimate genetic variance: sum 2pq sigmaSq
    public:
        varEffectScaled(const string &lab = "SigmaSqG"): Parameter(lab){};
        void compute(const double sigmaSq, const double sum2pq){value = sigmaSq*sum2pq;};
    };

    
public:
    const Data &data;
    
    VectorXd ycorr;   // corrected y for mcmc sampling
    VectorXd ghat;    // predicted total genotypic values
    VectorXd rhat;    // predicted total random covariate values
    
    bool estimatePi;
    
    FixedEffects fixedEffects;
    RandomEffects randomEffects;
    SnpEffects snpEffects;
    SnpPIP snpPip;
    VarEffects sigmaSq;
    VarRandomEffects sigmaSqRand;
    ScaleVar scale;
    Pi pi;
    ResidualVar vare;
    
    GenotypicVar varg;
    Heritability hsq;
    RandomVar varRand;
    Rounding rounding;
    NumNonZeroSnp nnzSnp;
    
    BayesC(const Data &data, const double varGenotypic, const double varResidual, const double varRandom, const double pival, const double piAlpha, const double piBeta, const bool estimatePi, const bool noscale,
           const string &algorithm = "Gibbs", const bool message = true):
    data(data),
    ycorr(data.y),
    fixedEffects(data.fixedEffectNames),
    randomEffects(data.randomEffectNames),
    sigmaSqRand(varRandom, data.numRandomEffects),
    snpEffects(data.snpEffectNames, algorithm),
    snpPip(data.snpEffectNames),
    sigmaSq(varGenotypic, data.snp2pq, pival, noscale),
    scale(sigmaSq.scale),
    pi(pival, piAlpha, piBeta),
    vare(varResidual, data.numKeptInds),
    varg(varGenotypic),
    varRand(varRandom),
    estimatePi(estimatePi)
    {
        numSnps = data.numIncdSnps;
        paramSetVec = {&snpEffects, &fixedEffects, &snpPip};           // for which collect mcmc samples
        paramVec = {&pi, &nnzSnp, &sigmaSq, &vare, &varg, &hsq};       // for which collect mcmc samples
        paramToPrint = {&pi, &nnzSnp, &sigmaSq, &vare, &varg, &hsq};   // print in order
        if (data.numRandomEffects) {
            paramSetVec.push_back(&randomEffects);
            paramVec.push_back(&sigmaSqRand);
            paramVec.push_back(&varRand);
            paramToPrint.push_back(&varRand);
        }
        paramToPrint.push_back(&rounding);
        if (message) {
            string alg = algorithm;
            if (alg!="HMC") alg = "Gibbs (default)";
            LOGGER << "\nBayesC model fitted. Algorithm: " << alg << "." << endl;
            LOGGER << "scale factor: " << sigmaSq.scale << endl;
        }
    }
    
    void sampleUnknowns(void);
    void sampleStartVal(void);
};


class BayesR : public BayesC {
    // Prior for snp efect pi_1 * N(0, 0) + pi_2 * N(0, sig^2_beta * gamma_2) + pi_3 * N(0, sig^2_beta * gamma_3) + pi_3 * N(0, sig^2_beta * gamma_4)
    // consider S as unknown to make inference on the relationship between MAF and effect size
public:
    
    class SnpEffects : public BayesC::SnpEffects {
    public:
        vector<vector<unsigned> > snpset;
        double sum2pq;
        SnpEffects(const vector<string> &header, const string &alg): BayesC::SnpEffects(header, "Gibbs"){
            sum2pq = 0.0;
        }
        
        void sampleFromFC(VectorXd &ycorr, const MatrixXd &Z, const VectorXd &ZPZdiag, const VectorXd &Rsqrt, const bool weightedRes,
                          const double sigmaSq, const VectorXd &pis,  const VectorXd &gamma,
                          const double vare, VectorXd &ghat, VectorXd &snpStore,
                          const double varg, const bool originalModel);
    };

    class ProbMixComps : public vector<Parameter*>, public Stat::Dirichlet {

        // prior probability of a snp being in any of the distributions effect has a dirichlet prior
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
    
    class VarEffects : public BayesC::VarEffects {
    public:
        VarEffects(const double vg, const VectorXd &snp2pq, const VectorXd &gamma, const VectorXd &pi, const bool noscale, const string &lab = "SigmaSq"):
        BayesC::VarEffects(vg, snp2pq, 1-pi[0], noscale, lab) {
            if (noscale) {
                value = vg / (snp2pq.sum() * gamma.dot(pi));  // derived from prior knowledge on Vg and pi
            } else {
                value = vg / (snp2pq.size() * gamma.dot(pi));  // derived from prior knowledge on Vg and pi
            }
            
            scale = (df-2)/df*value;
        }
        
        void computeScale(const double varg, const VectorXd &snp2pq, const VectorXd &gamma, const VectorXd &pi);
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
    VectorXd snpStore;   
    SnpEffects snpEffects;
    VarEffects sigmaSq;
    ProbMixComps Pis;
    VgMixComps Vgs;
    NumSnpMixComps numSnps;
    Gammas gamma;
    
    bool originalModel;

    BayesR(const Data &data, const double varGenotypic, const double varResidual, const double varRandom, const VectorXd pis, const VectorXd &piPar, const VectorXd gamma, const bool estimatePi, const bool noscale, const bool originalModel,
           const string &algorithm, const bool message = true):
    BayesC(data, varGenotypic, varResidual, varRandom, 1-pis[0], piPar[0], piPar[1], estimatePi, noscale, "Gibbs", false),
    Pis(pis, piPar),
    numSnps(pis),
    Vgs(gamma),
    gamma(gamma, vector<string>(gamma.size())),
    snpEffects(data.snpEffectNames, algorithm),
    sigmaSq(varGenotypic, data.snp2pq, gamma, pis, noscale),
    originalModel(originalModel)
    {
        paramSetVec  = {&snpEffects, &fixedEffects};
        for (unsigned i=0; i<Pis.size(); ++i) { 
           Pis[i]->value=Pis.values[i];
        }
        paramVec     = {&nnzSnp, &sigmaSq, &vare, &varg, &hsq};
        if (originalModel) paramVec.insert(paramVec.begin(), Vgs.begin(), Vgs.end());
        paramVec.insert(paramVec.begin(), numSnps.begin(), numSnps.end());
        paramToPrint = {&sigmaSq, &vare, &varg, &hsq};
        if (originalModel) paramToPrint.insert(paramToPrint.begin(), Vgs.begin(), Vgs.end());
        paramToPrint.insert(paramToPrint.begin(), numSnps.begin(), numSnps.end());
        if (data.numRandomEffects) {
            paramSetVec.push_back(&randomEffects);
            paramVec.push_back(&sigmaSqRand);
            paramVec.push_back(&varRand);
            paramToPrint.push_back(&varRand);
        }
        paramToPrint.push_back(&rounding);
        if (message) {
            string alg = algorithm;
            if (alg!="HMC") alg = "Gibbs (default)";
            LOGGER << "\nBayesR model fitted. Algorithm: " << alg << "." << endl;
            LOGGER << "scale factor: " << sigmaSq.scale << endl;
            // LOGGER << "Gamma: " << gamma.transpose() << endl;
        }
    }   
    void sampleUnknowns(void);
};
    


class ApproxBayesC : public BayesC {
public:
    class FixedEffects : public BayesC::FixedEffects {
    public:
        FixedEffects(const vector<string> &header): BayesC::FixedEffects(header){}
        
        void sampleFromFC(const MatrixXd &XPX, const VectorXd &XPXdiag,
                          const MatrixXd &ZPX, const VectorXd &XPy,
                          const VectorXd &snpEffects, const double vare,
                          VectorXd &rcorr);
    };
    
    class SnpEffects : public BayesC::SnpEffects {
    public:
        double sum2pq;
        
        VectorXd nnzPerChr;
        VectorXd nnzPerBlk;
        VectorXi leaveout;
        VectorXd ssqBlocks;
        VectorXi badSnps;

        SnpEffects(const vector<string> &header): BayesC::SnpEffects(header, "Gibbs"){
            sum2pq = 0.0;
            leaveout.setZero(size);
            badSnps.setZero(size);
        }
        
        void sampleFromFC(VectorXd &rcorr, const vector<SparseVector <double>  > &ZPZsp, const VectorXd &ZPZdiag, const VectorXd &ZPy,
                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                          const VectorXd &se, const VectorXd &tss, VectorXd &varei, const VectorXd &n, const VectorXd &snp2pq, const VectorXd &LDsamplVar,
                          const double sigmaSq, const double pi, const double vare, const double varg, const double ps, const double overdispersion);
        void sampleFromFC(VectorXd &rcorr, const vector<VectorXd> &ZPZ, const VectorXd &ZPZdiag, const VectorXd &ZPy,
                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                          const VectorXd &se, const VectorXd &tss, VectorXd &varei, const VectorXd &n, const VectorXd &snp2pq, const VectorXd &LDsamplVar,
                          const double sigmaSq, const double pi, const double vare, const double varg, const double ps, const double overdispersion);
        void hmcSampler(VectorXd &rcorr, const VectorXd &ZPy, const vector<VectorXd> &ZPZ,
                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                          const double sigmaSq, const double pi, const double vare);
        VectorXd gradientU(const VectorXd &effects, VectorXd &rcorr, const VectorXd &ZPy, const vector<VectorXd> &ZPZ,
                                                     const VectorXi &windStart, const VectorXi &windSize, const unsigned chrStart, const unsigned chrSize,
                                                     const double sigmaSq, const double vare);
        double computeU(const VectorXd &effects, const VectorXd &rcorr, const VectorXd &ZPy, const double sigmaSq, const double vare);
        
        void sampleFromFC(vector<VectorXd> &wcorrBlocks, const vector<MatrixXd> &Qblocks, vector<VectorXd> &whatBlocks,
                          const vector<LDBlockInfo*> keptLdBlockInfoVec, const VectorXd &nGWASblocks, const VectorXd &vareBlocks,
                          const double sigmaSq, const double pi, const double varg, const VectorXd &snp2pq);
    };
    
    class ResidualVar : public BayesC::ResidualVar {
    public:
        const double icrsq;
        
        ResidualVar(const double vare, const unsigned nobs, const double icrsq): BayesC::ResidualVar(vare, nobs), icrsq(icrsq) {}
        
        //void sampleFromFC(VectorXd &rcorr, const SpMat &ZPZinv);
//        void sampleFromFC(const double ypy, const VectorXd &effects, const VectorXd &ZPy, const VectorXd &rcorr, const double varg, const double nnz);
        void sampleFromFC(const double ypy, const VectorXd &effects, const VectorXd &ZPy, const VectorXd &rcorr, const double covg);
//        void sampleFromFCshrink(const double ypy, const VectorXd &effects, const VectorXd &ZPy, const VectorXd &rcorr, const double hsq, const double phi);
        
//        void sampleFromFC2(const double ypy, const VectorXd &effects, const VectorXd &ZPy, const VectorXd &ghat);
        
//        void randomWalkMHsampler(const double ypy, const VectorXd &effects, const VectorXd &ZPy, const VectorXd &rcorr, const VectorXd &ZPZrss, const double sigmaSq, const double pi);
    };
    
    class GenotypicVar : public BayesC::GenotypicVar {
    public:
        const unsigned nobs;
        
        GenotypicVar(const double varg, const unsigned n): BayesC::GenotypicVar(varg), nobs(n){}
//        void compute(const VectorXd &effects, const VectorXd &ZPy, const VectorXd &rcorr);
        void compute(const VectorXd &effects, const VectorXd &ZPy, const VectorXd &rcorr, const double covg);
    };
    
    class BlockGenotypicVar : public ParamSet, public BayesC::GenotypicVar {
    public:
        const unsigned nobs;
        unsigned numBlocks;
        double total;
        
        BlockGenotypicVar(const vector<string> &header, const double varg, const unsigned n, const string &lab = "BlockGenVar"):
        ParamSet(lab, header), BayesC::GenotypicVar(varg), nobs(n){
            numBlocks = header.size();
            total = 0.0;
        }

        void compute(const vector<VectorXd> &whatBlocks);
    };
    
    class BlockResidualVar : public ParamSet, public Stat::InvChiSq {
    public:
        const double df;      // hyperparameter
        const double scale;   // hyperparameter
        
        const double vary;

        unsigned numBlocks;
        double threshold;
        double mean;
        
        BlockResidualVar(const vector<string> &header, const double varPhenotypic, const string &lab = "BlockResVar"):
        ParamSet(lab, header), df(4), scale(0.5f*varPhenotypic), vary(varPhenotypic) {
            values.setConstant(size, varPhenotypic);
            numBlocks = header.size();
            threshold = 1.1;
            mean = varPhenotypic;
        }
        
        void sampleFromFC(vector<VectorXd> &wcorrBlocks, VectorXd &ssqBlocks, const VectorXd &nGWASblocks, const VectorXd &numEigenvalBlock);
    };


    class Rounding : public BayesC::Rounding {
    public:
        void computeRcorr(const VectorXd &ZPy, const vector<SparseVector <double>  > &ZPZsp,
                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                          const VectorXd &snpEffects, VectorXd &rcorr);
        void computeRcorr(const VectorXd &ZPy, const vector<VectorXd> &ZPZ,
                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                          const VectorXd &snpEffects, VectorXd &rcorr);
        void computeGhat(const MatrixXd &Z, const VectorXd &snpEffects, VectorXd &ghat);
    };
    
//    class Overdispersion : public BayesC::ResidualVar {
//    public:
//        Overdispersion(const double vare, const unsigned nobs): BayesC::ResidualVar(vare, nobs, "TauSq"){}
//        
//        void sampleFromFC(const VectorXd &y, const VectorXd &ghat);
//    };
    
    class PopulationStratification : public Parameter, public Stat::InvChiSq {
    public:
        const double df;
        const double scale;
        
        VectorXd chrSpecific;
        
        PopulationStratification(): Parameter("PS"), df(4), scale(0.5){
            chrSpecific.setZero(22);
        }
        
        void compute(const VectorXd &rcorr, const VectorXd &ZPZdiag, const VectorXd &LDsamplVar, const double varg, const double vare, const VectorXd &chisq);
        void compute(const VectorXd &rcorr, const VectorXd &ZPZdiag, const VectorXd &LDsamplVar, const double varg, const double vare, const vector<ChromInfo*> chromInfoVec);
    };
    
    class NumResidualOutlier : public Parameter {
    public:
        ofstream out;
        unsigned iter;
        
        NumResidualOutlier(): Parameter("Nro"){
            iter = 0;
        }
        
        void compute(const VectorXd &rcorr, const VectorXd &ZPZdiag, const VectorXd &LDsamplVar, const double varg, const double vare, const vector<string> &snpName, VectorXi &leaveout, const vector<SparseVector <double>  > &ZPZ, const VectorXd &ZPy, const VectorXd &snpEffects);
    };
    
    class InterChrGenetCov : public Parameter {
    public:
        const double spouseCorrelation;
        const unsigned nobs;
        
        InterChrGenetCov(const double corr, const unsigned nobs): Parameter("GenCov"), spouseCorrelation(corr), nobs(nobs) {}
        
        void compute(const double ypy, const VectorXd &effects, const VectorXd &ZPy, const VectorXd &rcorr);
    };
    
    class NnzGwas : public Parameter {
    public:
        unsigned iter;
        
        NnzGwas(): Parameter("NnzGwas"){
            iter = 0;
        }
        
        void compute(const VectorXd &effects, const vector<SparseVector <double>  > &ZPZ, const VectorXd &ZPZdiag);
    };
    
    class PiGwas : public Parameter {
    public:
        unsigned iter;
        PiGwas(): Parameter("PiGwas"){
            iter = 0;
        }
        
        void compute(const double nnzGwas, const unsigned numSnps);
    };


    class NumBadSnps : public Parameter {
    public:
        double betaThresh;
        vector<string> snpNames;
        vector<string> badSnpName;
        vector<unsigned> badSnpIdx;
        ofstream out;
        
        bool writeTxt;
        
        NumBadSnps(const string &title, const VectorXd &b, const vector<string> &snpNames): Parameter("NumSkeptSnp"), snpNames(snpNames){
            VectorXd abs_b = b.array().abs();
            std::sort(abs_b.data(), abs_b.data() + abs_b.size());
            int index8 = 0.8 * (abs_b.size() - 1);
            betaThresh = abs_b[index8];
            //cout << "Set beta cutoff threshold: " << betaThresh << endl;
            //cout << b.head(10) << endl;
            
            string filename = title + ".skepticalSNPs";
            out.open(filename.c_str());
            writeTxt = true;
        }
        void compute_eigen(VectorXi &badSnps, VectorXd &effects, VectorXd &effectMean, const VectorXd &b, vector<VectorXd> &wcorrBlocks, const vector<MatrixXd> &Qblocks, const vector<LDBlockInfo*> keptLdBlockInfoVec, const int iter);
    };
    

    
public:
    const Data &data;
    const double phi;   // the shrinkage parameter for heritability estimate
    const double overdispersion;
    
    VectorXd rcorr;
    VectorXd varei;   // residual variance specific to each snp
    
    vector <double>  hsqMCMC;
    
    bool sparse;
    bool modelPS;
    bool diagnose;
    bool robustMode;
    
    FixedEffects fixedEffects;
    SnpEffects snpEffects;
    BayesC::VarEffects sigmaSq;
    BayesC::Pi pi;
    ResidualVar vare;
    GenotypicVar varg;
//    BayesC::ResidualVar vare;
    Rounding rounding;
    NumBadSnps nBadSnps;
    varEffectScaled sigmaSqG;
//    Overdispersion tauSq;
    PopulationStratification ps;
    NumResidualOutlier nro;
    InterChrGenetCov covg;
    PiGwas pigwas;
    NnzGwas nnzgwas;
    double genVarPrior;
    double scalePrior;
    bool noscale;
    bool lowRankModel;

    BlockGenotypicVar vargBlk;
    BlockResidualVar vareBlk;
    
    vector<VectorXd> wcorrBlocks;
    vector<VectorXd> whatBlocks;
   
    ApproxBayesC(const Data &data, const bool lowrank, const double varGenotypic, const double varResidual, const double varRandom, const double pival, const double piAlpha, const double piBeta, const bool estimatePi, const bool noscale,
                 const double phi, const double overdispersion, const bool estimatePS, const double icrsq, const double spouseCorrelation,
                 const bool diagnosticMode, const bool robustMode, const bool randomStart = false, const bool message = true)
    : BayesC(data, varGenotypic, varResidual, 0.0, pival, piAlpha, piBeta, estimatePi, noscale, "Gibbs", false)
    , data(data)
    , rcorr(data.ZPy)
    , wcorrBlocks(data.wcorrBlocks)
    , varei(data.tss.array()/data.n.array())
    , fixedEffects(data.fixedEffectNames)
    , snpEffects(data.snpEffectNames)
    , sigmaSq(varGenotypic, data.snp2pq, pival, noscale)
    , pi(pival, piAlpha, piBeta)
    , genVarPrior(varGenotypic)
    , noscale(noscale)
    , scalePrior(sigmaSq.scale)
    , vare(varResidual, data.numKeptInds, icrsq)
    , varg(varGenotypic, data.numKeptInds)
//    , tauSq(varResidual, data.numKeptInds)
    , phi(phi)
    , overdispersion(overdispersion)
    , covg(spouseCorrelation, data.numKeptInds)
    , robustMode(robustMode)
    , vargBlk(data.ldblockNames, varGenotypic, data.numKeptInds)
    , vareBlk(data.ldblockNames, data.varPhenotypic)
    , nBadSnps(data.title, data.b, data.snpEffectNames)
    , lowRankModel(lowrank)
    {
        sparse = data.sparseLDM;
        modelPS = estimatePS;
        diagnose = diagnosticMode;
        paramSetVec = {&snpEffects};
        paramVec = {&pi, &nnzSnp, &sigmaSq, &vare, &varg, &sigmaSqG, &hsq};
        paramToPrint = {&pi, &nnzSnp, &sigmaSq, &vare, &varg, &sigmaSqG, &hsq};
//        if (sparse) {
//            paramVec.push_back(&pigwas);
//            paramVec.push_back(&nnzgwas);
//            paramToPrint.push_back(&pigwas);
//            paramToPrint.push_back(&nnzgwas);
//        }
        if (lowRankModel) {
            paramSetVec.push_back(&vargBlk);
            paramSetVec.push_back(&vareBlk);
        }
        if (modelPS) {
            paramVec.push_back(&ps);
            paramToPrint.push_back(&ps);
        }
        if (diagnose) {
            nro.out.open((data.label+".diagnostics").c_str());
            paramVec.push_back(&nro);
            paramToPrint.push_back(&nro);
        }
        if (spouseCorrelation) {
            paramVec.push_back(&covg);
            paramToPrint.push_back(&covg);
        }
        if (message) {
            LOGGER << "\nSBayesC" << endl;
            if (lowRankModel) {
                LOGGER << "Using the low-rank model" << endl;
            }
            LOGGER << "scale factor: " << sigmaSq.scale << endl;
            if (noscale)
            {
               LOGGER << "Fitting model assuming unscaled genotypes " << endl; 
            } else
            {
               LOGGER << "Fitting model assuming scaled genotypes "  << endl;
            }
            if (robustMode) LOGGER << "Using a more robust parameterisation " << endl;
        }
        if (randomStart) {
            LOGGER << "Sampling starting values for parameters..." << endl;
            sampleStartVal();
        }
    }
    
    void sampleUnknowns(void);
    static void ldScoreReg(const VectorXd &chisq, const VectorXd &LDscore, const VectorXd &LDsamplVar,
                           const double varg, const double vare, double &ps);
    void checkHsq(vector <double>  &hsqMCMC);
};


// -----------------------------------------------------------------------------------------------
// Approximate Bayes R
// -----------------------------------------------------------------------------------------------

class ApproxBayesR : public ApproxBayesC {
    
public:
    
    class SnpEffects : public ApproxBayesC::SnpEffects {
    public:
        vector<vector<unsigned> > snpset;
        VectorXd deltaNZ;
        VectorXd lambdaVec;
        VectorXd uhatVec;
        VectorXd invGammaVec;
        vector<unsigned> deltaNzIdx;
        double sum2pq;

        SnpEffects(const vector<string> &header): ApproxBayesC::SnpEffects(header){
            sum2pq = 0.0;
            deltaNZ.setZero(size);
        }
        
        void sampleFromFC(VectorXd &rcorr, const vector<SparseVector<double> > &ZPZsp, const VectorXd &ZPZdiag, const VectorXd &ZPy,
                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                          const VectorXd &se, const VectorXd &tss, VectorXd &varei, const VectorXd &n, const VectorXd &snp2pq, const VectorXd &LDsamplVar,
                          const double sigmaSq, const VectorXd &pis, const VectorXd &gamma, const double vare, VectorXd &snpStore, 
                          const double varg, const double ps, const double overdispersion,
                          const bool originalModel);
        void sampleFromFC(VectorXd &rcorr, const vector<VectorXd> &ZPZ, const VectorXd &ZPZdiag, const VectorXd &ZPy,
                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                          const VectorXd &se, const VectorXd &tss, VectorXd &varei, const VectorXd &n, const VectorXd &snp2pq, const VectorXd &LDsamplVar,
                          const double sigmaSq, const VectorXd &pis, const VectorXd &gamma, const double vare, VectorXd &snpStore,
                          const double varg, const double ps, const double overdispersion,
                          const bool originalModel);
        
        void sampleFromFC(const VectorXd &ZPy, const SpMat &ZPZsp, const VectorXd &ZPZdiag,
                          VectorXd &rcorr, const VectorXd &LDsamplVar,
                          const double sigmaSq, const VectorXd &pis, const VectorXd &gamma, VectorXd &snpStore,
                          const double varg, const double vare, const double ps, const double overdispersion, const bool originalModel);
        
        void sampleFromFC(const VectorXd &ZPy, const VectorXd &ZPZdiag, const MatrixXd &Z, const double n_ref, const double n_gwas,
                          const double sigmaSq, const VectorXd &pis, const VectorXd &gamma, const double vare,
                          VectorXd &snpStore, VectorXd &ghat, const double varg, const bool originalModel);

        void sampleFromFC(vector<VectorXd> &wcorrBlocks, const vector<MatrixXd> &Qblocks, vector<VectorXd> &whatBlocks,
                          const vector<LDBlockInfo*> keptLdBlockInfoVec, const VectorXd &nGWASblocks, const VectorXd &vareBlocks,
                          const double sigmaSq, const VectorXd &pis, const VectorXd &gamma, VectorXd &snpStore, const double varg,
                          const bool originalModel);

        void adjustByCG(const VectorXd &ZPy, const vector<SparseVector<double> > &ZPZsp, VectorXd &rcorr);
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
        
        void compute(const VectorXd &snpEffects, const VectorXd &ZPy, const VectorXd &rcorr, const vector<vector<unsigned> > snpset, const double varg, const double nobs);
        //void compute(const VectorXd &snpEffects, const vector<SparseVector<double> > &ZPZsp, const vector<vector<unsigned> > snpset, const double varg, const double nobs);
        //void compute(const VectorXd &snpEffects, const vector<VectorXd> &ZPZ, const vector<vector<unsigned> > snpset, const double varg, const double nobs);
    };
    
    
    VectorXd snpStore;   
    SnpEffects snpEffects;
    ApproxBayesC::FixedEffects fixedEffects;
    ApproxBayesC::ResidualVar vare;
    ApproxBayesC::GenotypicVar varg;
    ApproxBayesC::Rounding rounding;
    varEffectScaled sigmaSqG;
    ApproxBayesC::PopulationStratification ps;
    ApproxBayesC::NumResidualOutlier nro;
    ApproxBayesC::InterChrGenetCov covg;
    ApproxBayesC::PiGwas pigwas;
    ApproxBayesC::NnzGwas nnzgwas;
    BayesR::ProbMixComps Pis;
    BayesR::NumSnpMixComps numSnps;
    VgMixComps Vgs;
    ApproxBayesC::BlockGenotypicVar vargBlk;
    ApproxBayesC::BlockResidualVar vareBlk;
    
    BayesR::Gammas gamma;
    double genVarPrior;
    double scalePrior;
    bool noscale;
    bool originalModel;
    bool estimateSigmaSq;
    bool estimateHsq;

    const double overdispersion;
    
    enum {gibbs, cg, mh} algorithm;
        
    
    ApproxBayesR(const Data &data, const bool lowrank, const double varGenotypic, const double varResidual, const VectorXd pis, const VectorXd &piPar, const VectorXd gamma, const bool estimatePi, const bool estimateSigmaSq, const bool noscale, const bool originalModel, const double overdispersion, const bool estimatePS, const double spouseCorrelation, const bool diagnosticMode, const bool robustMode, const string &alg, const bool message = true):
    ApproxBayesC(data, lowrank, varGenotypic, varResidual, 0.0, (1-pis[0]), piPar[0], piPar[1], estimatePi, noscale, 0, overdispersion, estimatePS, 0, spouseCorrelation, diagnosticMode, robustMode, false, false),
    Pis(pis,piPar),
    numSnps(pis),
    Vgs(gamma),
    gamma(gamma, vector<string>(gamma.size())),
    vare(varResidual, data.numKeptInds, 0),
    varg(varGenotypic, data.numKeptInds),
    fixedEffects(data.fixedEffectNames),
    snpEffects(data.snpEffectNames),
    genVarPrior(varGenotypic),
    noscale(noscale),
    overdispersion(overdispersion),
    covg(spouseCorrelation, data.numKeptInds),
    scalePrior(sigmaSq.scale),
    originalModel(originalModel),
    estimateSigmaSq(estimateSigmaSq),
    vargBlk(data.ldblockNames, varGenotypic, data.numKeptInds),
    vareBlk(data.ldblockNames, data.varPhenotypic)
    {
        if (alg == "cg") algorithm = cg;
        else if (alg == "MH") algorithm = mh;
        else algorithm = gibbs;
        sparse = data.sparseLDM;
        // varg.value = varGenotypic; //// NOTE: write it into constructor!!!
        paramSetVec = {&snpEffects, &fixedEffects};
        // sigmaSq.value = varGenotypic/(data.snp2pq.array().sum()*(1-pis[0]));
        // scale.value = sigmaSq.scale = 0.5*sigmaSq.value;
        for (unsigned i=0; i<Pis.size(); ++i) { 
           Pis[i]->value=Pis.values[i];  
        }
        paramVec     = {&nnzSnp, &sigmaSq, &vare, &varg, &hsq};
        if (originalModel) paramVec.insert(paramVec.begin(), Vgs.begin(), Vgs.end());
        paramVec.insert(paramVec.begin(), numSnps.begin(), numSnps.end());
        paramToPrint = {&sigmaSq, &vare, &varg, &hsq};
        if (originalModel) paramToPrint.insert(paramToPrint.begin(), Vgs.begin(), Vgs.end());
        paramToPrint.insert(paramToPrint.begin(), numSnps.begin(), numSnps.end());
        if (lowRankModel) {
            paramSetVec.push_back(&vargBlk);
            paramSetVec.push_back(&vareBlk);
        }
        if (modelPS) {
            paramVec.push_back(&ps);
            paramToPrint.push_back(&ps);
        }
        if (diagnose) {
            nro.out.open((data.label+".diagnostics").c_str());
            paramVec.push_back(&nro);
            paramToPrint.push_back(&nro);
        }

        if (data.Z.size()) ghat.setZero(data.Z.rows());   // TMP_JZ
        else ghat.resize(0);                              // TMP_JZ
        
        if (!estimateSigmaSq) {
            sigmaSq.value = varg.value/(data.numIncdSnps*pis.dot(gamma));
            //sigmaSq.value = 0.1369346;
            sigmaSq.scale = 0.5*sigmaSq.value;
            LOGGER << "fixing sigmaSq to be " << sigmaSq.value << endl;
        }
        if (!estimateSigmaSq && (originalModel || robustMode)) {
            estimateHsq = false;
            hsq.value = varg.value/(varg.value + vare.value);
        } else estimateHsq = true;
        
        if (message) {
            LOGGER << "\nSBayesR" << endl;
            if (lowRankModel) {
                LOGGER << "Using the low-rank model" << endl;
            }
//            LOGGER << "scale factor: " << sigmaSq.scale << endl;
            LOGGER << "Gamma: " << gamma.transpose() << endl;
            if (noscale)
            {
               LOGGER << "Fitting model assuming unscaled genotypes " << endl; 
            } else
            {
               LOGGER << "Fitting model assuming scaled genotypes "  << endl;
            }
            if (robustMode) LOGGER << "Using a more robust parameterisation " << endl;
            if (algorithm == cg) LOGGER << "Conjugate gradient-adjusted Gibbs sampling" << endl;
        }
    }
    
    void sampleUnknowns(void);
};


#endif /* model_hpp */
