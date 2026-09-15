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


#ifndef mcmc_hpp
#define mcmc_hpp
#include <complex>
// #define lapack_complex_double std::complex <double> 
#define lapack_complex_double std::complex<double>

#include <cstdio>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <boost/format.hpp>
#include "Model.hpp"
#include "Gadgets.hpp"

using namespace std;
using namespace Eigen;


class McmcSamples {
    // rows: MCMC cycles, cols: model parameters
public:
    const string label;
    string filename;
    enum {dense, sparse, stream} storageMode;
    VectorXd onlineM2, positiveCount, enrichmentCount;
    
    unsigned chainLength;
    unsigned burnin;
    unsigned thin;
    
    unsigned nrow;
    unsigned ncol;
    unsigned nnz;  // number of non-zeros for sparse matrix
    
    // datMat stores various parameters; if parameter is single type, datMat is iteration/thin * 1 matrix,
    // If parameter is vector type, datMat is iteration/thin * ncol(npar) matrix.
    // If parameter is matrix type, datMat is iteration/thin *ncor matrix. In this situation, I have no idea
    // about how to design getSample function and will deal with this situation later.
    MatrixXd datMat; 
    SpMat datMatSp; // most of the snp effects will be zero if pi value is high
    
    VectorXd posteriorMean;
    VectorXd posteriorSqrMean;
    VectorXd pip; 
    VectorXd lastSample; // save the last sample of MCMC
    
    FILE *bout = nullptr;
    ofstream tout;
    
    McmcSamples(const string &label, const unsigned chainLength, const unsigned burnin, const unsigned thin,
                const unsigned npar, const string &storage_mode = "dense"):
    label(label), chainLength(chainLength), burnin(burnin), thin(thin) {
        if (!thin || burnin>=chainLength) throw std::invalid_argument("MCMC requires thin > 0 and burn-in < chain length");
        nrow = (uint64_t(chainLength)+thin-1)/thin - (uint64_t(burnin)+thin-1)/thin;
        if(!nrow) throw std::invalid_argument("No retained MCMC samples for these settings");
        ncol = npar;
        if (storage_mode == "stream") {
            storageMode=stream;onlineM2.setZero(ncol);positiveCount.setZero(ncol);enrichmentCount.setZero(ncol);
        } else if (storage_mode == "dense") {
            storageMode = dense;
            datMat.setZero(nrow, ncol);
        } else if (storage_mode == "sparse") {
            storageMode = sparse;
            datMatSp.resize(nrow,ncol);  // for faster filling the matrix
        } else {
            cerr << "Error: Unrecognized storage mode: " << storage_mode << endl;
        }
        posteriorMean.setZero(ncol);
        posteriorSqrMean.setZero(ncol);
        pip.setZero(ncol);
        lastSample.setZero(ncol);
    }
    
    McmcSamples(const string &label): label(label) {}
    // the function of getSample is to 
    // for single value
    void getSample(const unsigned iter, const double sample, const bool writeTxtPosterior, ofstream &out);
    // for vector
    void getSample(const unsigned iter, const VectorXd &sample, const bool writeBinPosterior, const bool writeTxtPosterior);
    // for matix
    void getSample(const unsigned iter, const MatrixXd &sample, const bool writeBinPosterior, const bool writeTxtPosterior);

    VectorXd mean(void);
    VectorXd sd(void);
    
    void initBinFile(const string &title);
    void writeBinSample(unsigned iter,const VectorXd &sample);
    void finishFiles();
    void initTxtFile(const string &title);
    void writeDataBin(const string &title);
    void writeDataTxt(const string &title);
    void readDataBin(const string &filename);
    void readDataTxt(const string &filename);
    void readDataTxt(const string &filename, const string &label);
};

class MCMC {
private:
    string outfilename;
    ofstream out;
    
    void initTxtFile(const vector<Parameter*> &paramVec, const string &title);
    vector<McmcSamples*> initMcmcSamples(const Model &model, const unsigned chainLength, const unsigned burnin,
                                         const unsigned thin, const string &title, const bool writeBinPosterior, const bool writeTxtPosterior);
    void collectSamples(const Model &model, vector<McmcSamples*> &mcmcSampleVec, const unsigned iteration, const bool writeBinPosterior, const bool writeTxtPosterior);
    void printStatus(const vector<Parameter*> &paramToPrint, const unsigned thisIter, const unsigned outputFreq, const string &timeLeft);
    void printSummary(const vector<Parameter*> &paramToPrint, const vector<McmcSamples*> &mcmcSampleVec, const string &filename);
    void printSetSummary(const vector<ParamSet*> &paramSetToPrint, const vector<McmcSamples*> &mcmcSampleVec, const string &filename);
public:
    vector<McmcSamples*> run(Model &model, const unsigned chainLength, const unsigned burnin, const unsigned thin, const bool print,
                             const unsigned outputFreq, const string &title, const bool writeBinPosterior, const bool writeTxtPosterior);
    void convergeDiagGelmanRubin(const Model &model, vector<vector<McmcSamples*> > &mcmcSampleVecChain, const string &filename);
};

#endif /* mcmc_hpp */
