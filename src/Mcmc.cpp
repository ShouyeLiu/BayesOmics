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


#include "Mcmc.hpp"

// for  single parameter
void McmcSamples::getSample(const unsigned iter, const double sample, const bool writeTxtPosterior, ofstream &out){
    if (writeTxtPosterior) out << boost::format("%12s ") %sample;
    if (iter % thin) return;
    unsigned thin_iter_post_burnin = iter/thin - burnin/thin;
    if (iter >= burnin) {
        datMat(thin_iter_post_burnin,0) = sample;
        // LOGGER << "thin: " << thin_iter_post_burnin << " sample: " << sample << " ";
        posteriorMean.array() += (sample - posteriorMean.array())/(thin_iter_post_burnin+1);
        posteriorSqrMean.array() += (sample*sample - posteriorSqrMean.array())/(thin_iter_post_burnin+1);
    }
}

//  for vector parameter
void McmcSamples::getSample(const unsigned iter, const VectorXd &sample, const bool writeBinPosterior, const bool writeTxtPosterior){
    if (storageMode == dense) {
        if (writeTxtPosterior) tout << sample.transpose() << endl;
    }

    if (iter % thin) return;
    unsigned thin_iter = iter/thin;
    unsigned thin_iter_post_burnin = thin_iter - burnin/thin;
    if (storageMode == dense) {
        ArrayXd delta;
        delta.setZero(sample.size());
        for(unsigned i =0; i< sample.size();i++){
            if(abs(sample[i]) > 5e-6) delta[i] = 1;
        }

        if (iter >= burnin) {
            datMat.row(thin_iter_post_burnin) = sample;
            pip.array() += (delta - pip.array())/(thin_iter_post_burnin+1);
            posteriorMean.array() += (sample - posteriorMean).array()/(thin_iter_post_burnin+1);
            posteriorSqrMean.array() += (sample.array().square() - posteriorSqrMean.array())/(thin_iter_post_burnin+1);
        }
    } else if (storageMode == sparse) {
        lastSample = sample;
        SparseVector <double>  spvec = sample.sparseView();
        if (writeBinPosterior) {
            for (SparseVector <double> ::InnerIterator it(spvec); it; ++it) {
                unsigned rc[2] = {thin_iter, (unsigned)it.index()};
                fwrite(rc, sizeof(unsigned), 2, bout);
                double val = it.value();
                fwrite(&val, sizeof(double), 1, bout);
            }
        }
        ArrayXd delta;
        delta.setZero(sample.size());
        for (SparseVector <double> ::InnerIterator it(spvec); it; ++it) {
            delta[it.index()] = 1;
        }
        if (iter >= burnin) {
            pip.array() += (delta - pip.array())/(thin_iter_post_burnin+1);
            posteriorMean.array() += (sample - posteriorMean).array()/(thin_iter_post_burnin+1);
            posteriorSqrMean.array() += (sample.array().square() - posteriorSqrMean.array())/(thin_iter_post_burnin+1);
        }
    }
}

// for matrix parameter
void McmcSamples::getSample(const unsigned iter, const MatrixXd &sample, const bool writeBinPosterior, const bool writeTxtPosterior){
    // if (storageMode == dense) {
    //     if (writeTxtPosterior) tout << sample.transpose() << endl;
    // }
    // if (iter % thin) return;
    // unsigned thin_iter = iter/thin;
    // unsigned thin_iter_post_burnin = thin_iter - burnin/thin;
    
    // if (storageMode == dense) {
    //     if (iter >= burnin) {
    //         datMat.row(thin_iter_post_burnin) = sample;
    //         posteriorMean.array() += (sample - posteriorMean).array()/(thin_iter_post_burnin+1);
    //         posteriorSqrMean.array() += (sample.array().square() - posteriorSqrMean.array())/(thin_iter_post_burnin+1);
    //     }
    // } else if (storageMode == sparse) {
    //     lastSample = sample;
    //     SparseVector <double>  spvec = sample.sparseView();
    //     if (writeBinPosterior) {
    //         for (SparseVector <double> ::InnerIterator it(spvec); it; ++it) {
    //             unsigned rc[2] = {thin_iter, (unsigned)it.index()};
    //             fwrite(rc, sizeof(unsigned), 2, bout);
    //             double val = it.value();
    //             fwrite(&val, sizeof(double), 1, bout);
    //         }
    //     }
    //     // ArrayXd delta;
    //     // delta.setZero(sample.size());
    //     // for (SparseVector <double> ::InnerIterator it(spvec); it; ++it) {
    //     //     delta[it.index()] = 1;
    //     // }
    //     // if (iter >= burnin) {
    //     //     pip.array() += (delta - pip.array())/(thin_iter_post_burnin+1);
    //     //     posteriorMean.array() += (sample - posteriorMean).array()/(thin_iter_post_burnin+1);
    //     //     posteriorSqrMean.array() += (sample.array().square() - posteriorSqrMean.array())/(thin_iter_post_burnin+1);
    //     // }
    // }
}


VectorXd McmcSamples::mean(){
    if (storageMode == dense) {
        return VectorXd::Ones(nrow).transpose()*datMat/nrow;
    } else {
        return VectorXd::Ones(nrow).transpose()*datMatSp/nrow;
    }
}

VectorXd McmcSamples::sd(){
    VectorXd res(ncol);
    if (storageMode == dense) {
        for (unsigned i=0; i<ncol; ++i) {
            res[i] = std::sqrt(Gadget::calcVariance(datMat.col(i)));
        }
    } else {
        for (unsigned i=0; i<ncol; ++i) {
            res[i] = std::sqrt(Gadget::calcVariance(datMatSp.col(i)));
        }
    }
    return res;
}

void McmcSamples::initBinFile(const string &title){
    string dirname = title + ".mcmcsamples";
    if (!Gadget::directoryExist(dirname)) {
        LOGGER.e(0," cannot find directory " + dirname);
    }
    filename = dirname + "/" + label + ".mcmcsamples.bin";
    bout = fopen(filename.c_str(), "wb");
    if (!bout) {
        LOGGER.e(0," cannot open file " + filename);
    }
    nnz = 0;
    unsigned xyn[3] = {chainLength/thin, ncol, nnz};
    fwrite(xyn, sizeof(unsigned), 3, bout);
}

void McmcSamples::initTxtFile(const string &title){
    string dirname = title + ".mcmcsamples";
    if (!Gadget::directoryExist(dirname)) {
        LOGGER.e(0," cannot find directory " + dirname);
    }
    filename = dirname + "/" + label + ".mcmcsamples.txt";
    tout.open(filename.c_str());
    if (!tout) {
        LOGGER.e(0," cannot open file " + filename);
    }
}

void McmcSamples::writeDataBin(const string &title){
    string dirname = title + ".mcmcsamples";
    if (!Gadget::directoryExist(dirname)) {
        LOGGER.e(0," cannot find directory " + dirname);
    }
    filename = dirname + "/" + label + ".mcmcsamples.bin";
    FILE *out = fopen(filename.c_str(), "wb");
    if (!out) {
        LOGGER.e(0," cannot open file " + filename);
    }
    
    int xyn[3] = {static_cast<int>(datMatSp.rows()), static_cast<int>(datMatSp.cols()), static_cast<int>(datMatSp.nonZeros())};
    fwrite(xyn, sizeof(unsigned), 3, out);
    
    for (int i=0; i < datMatSp.outerSize(); ++i) {
        SpMat::InnerIterator it(datMatSp, i);
        for (; it; ++it) {
            unsigned rc[2] = {(unsigned)it.row(), (unsigned)it.col()};
            fwrite(rc, sizeof(unsigned), 2, out);
            double v = it.value();
            fwrite(&v, sizeof(double), 1, out);
        }
    }
    fclose(out);
}

void McmcSamples::readDataBin(const string &title){
    string dirname = title + ".mcmcsamples";
    if (!Gadget::directoryExist(dirname)) {
        LOGGER.e(0," cannot find directory " + dirname);
    }
    filename = dirname + "/" + label + ".mcmcsamples.bin";
    FILE *in = fopen(filename.c_str(), "rb");
    if (!in) {
        LOGGER.e(0," cannot open file " + filename);
    }
    
    unsigned xyn[3];
    fread(xyn, sizeof(unsigned), 3, in);
    
    nrow = xyn[0];
    ncol = xyn[1];
        
    datMatSp.resize(xyn[0], xyn[1]);
//    vector<Triplet <double> > trips(xyn[2]);
    vector<Triplet <double>  > trips;
    
    //for (int i=0; i < trips.size(); ++i){
    while (!feof(in)) {
        unsigned rc[2];
        fread(rc, sizeof(unsigned), 2, in);
        double v;
        fread(&v, sizeof(double), 1, in);
        
        if(rc[0]>xyn[0] || rc[1]>xyn[1]) continue;
        
        //trips[i] = Triplet <double> (rc[0], rc[1], v);
        trips.push_back(Triplet <double> (rc[0], rc[1], v));
    }
    fclose(in);
    
    datMatSp.setFromTriplets(trips.begin(), trips.end());
    datMatSp.makeCompressed();
    
    //LOGGER << "nrow: " << nrow << " ncol: " << ncol << " nonzeros: " << datMatSp.nonZeros() << " " << nnz << endl;
    //LOGGER << MatrixXd(datMatSp) << endl;
    
    storageMode = sparse;
}

void McmcSamples::readDataTxt(const string &title){
    string dirname = title + ".mcmcsamples";
    if (!Gadget::directoryExist(dirname)) {
        LOGGER.e(0," cannot find directory " + dirname);
    }
    filename = dirname + "/" + label + ".mcmcsamples.txt";
    ifstream in(filename.c_str());
    string inputStr;
    vector <double>  tmp;
    while (in >> inputStr) {
        tmp.push_back(stof(inputStr));
    }
    in.close();
    nrow = tmp.size();
    datMat.resize(nrow, 1);
    datMat.col(0) = Eigen::Map<VectorXd>(&tmp[0], nrow);
    storageMode = dense;
}

void McmcSamples::readDataTxt(const string &title, const string &label){
    string dirname = title + ".mcmcsamples";
    if (!Gadget::directoryExist(dirname)) {
        LOGGER.e(0," cannot find directory " + dirname);
    }
    filename = dirname + "/" + label + ".mcmcsamples.txt";
    ifstream in(filename.c_str());
    Gadget::Tokenizer colData;
    Gadget::Tokenizer header;
    string inputStr;
    string sep(" \t");
    vector <double>  tmp;
    unsigned line = 0;
    
    std::getline(in, inputStr);
    header.getTokens(inputStr, sep);
    int idx = header.getIndex(label);
    
    if (idx==-1) LOGGER.e(0," Cannot find " + label + " in file [" + filename + "].");
    
    while (getline(in, inputStr)) {
        ++line;
        colData.getTokens(inputStr, sep);
        tmp.push_back(stof(colData[idx]));
    }
    in.close();
    nrow = tmp.size();
    datMat.resize(nrow, 1);
    datMat.col(0) = Eigen::Map<VectorXd>(&tmp[0], nrow);
    storageMode = dense;
}

void McmcSamples::writeDataTxt(const string &title){
    string dirname = title + ".mcmcsamples";
    if (!Gadget::directoryExist(dirname)) {
        LOGGER.e(0," cannot find directory " + dirname);
    }
    filename = dirname + "/" + label + ".mcmcsamples.txt";
    ofstream out(filename);
    out << datMat << endl;
    out.close();
}

void MCMC::initTxtFile(const vector<Parameter*> &paramVec, const string &title){
    string dirname = title + ".mcmcsamples";
    if (!Gadget::directoryExist(dirname)) {
        LOGGER.e(0," cannot find directory " + dirname);
    }
    outfilename = dirname + "/CoreParameters.mcmcsamples.txt";
    out.open(outfilename.c_str());
    if (!out) {
        LOGGER.e(0," cannot open file " + outfilename);
    }
    for (unsigned i=0; i<paramVec.size(); ++i) {
        Parameter *par = paramVec[i];
        out << boost::format("%12s ") %par->label;
    }
    out << endl;
}

vector<McmcSamples*> MCMC::initMcmcSamples(const Model &model, const unsigned chainLength, const unsigned burnin, const unsigned thin,
                                           const string &title, const bool writeBinPosterior, const bool writeTxtPosterior){
    vector<McmcSamples*> mcmcSampleVec;
    // here we will assign a McmcSamples for each parameter (single, vector and matrix)
    // Step 1. parameter with single value(double);
    // e.g., vare, pi; datMat is iteration/thin * 1 matrix
    for (unsigned i=0; i<model.paramVec.size(); ++i) {
        Parameter *par = model.paramVec[i];
        McmcSamples *mcmcSamples = new McmcSamples(par->label, chainLength, burnin, thin, 1);
        if (writeTxtPosterior) mcmcSamples->initTxtFile(title);
        mcmcSampleVec.push_back(mcmcSamples);
    }
    // Step 2. parameter with many values(vector): 
    // e.g., delta, beta 
    for (unsigned i=0; i<model.paramSetVec.size(); ++i) {
        ParamSet *parSet = model.paramSetVec[i];
        McmcSamples *mcmcSamples;
        if (parSet->label.find("SnpEffects") != string::npos) {
            mcmcSamples = new McmcSamples(parSet->label, chainLength, burnin, thin, parSet->size, "dense");
            if (writeBinPosterior) mcmcSamples->initBinFile(title);
            if (writeTxtPosterior) mcmcSamples->initTxtFile(title);
        } else if (parSet->label.find("GeneEffects") != string::npos) {
            mcmcSamples = new McmcSamples(parSet->label, chainLength, burnin, thin, parSet->size, "dense");
            if (writeBinPosterior) mcmcSamples->initBinFile(title);
            if (writeTxtPosterior) mcmcSamples->initTxtFile(title);
        } else if (parSet->label.find("EQTLJointVec") != string::npos) {
            mcmcSamples = new McmcSamples(parSet->label, chainLength, burnin, thin, parSet->size, "dense");
            if (writeBinPosterior) mcmcSamples->initBinFile(title);
            if (writeTxtPosterior) mcmcSamples->initTxtFile(title);
        }

        ////////////////////////
        // if (parSet->label.find("SnpEffects") != string::npos) {
        //     mcmcSamples = new McmcSamples(parSet->label, chainLength, burnin, thin, parSet->size, "sparse");
        //     if (writeBinPosterior) mcmcSamples->initBinFile(title);
        //     // mcmcSamples->initTxtFile(title);
        // } else if (parSet->label.find("GeneEffects") != string::npos) {
        //     mcmcSamples = new McmcSamples(parSet->label, chainLength, burnin, thin, parSet->size, "dense");
        //     // if (writeBinPosterior) mcmcSamples->initBinFile(title);
        //     mcmcSamples->initTxtFile(title);
        // // } else if (parSet->label.find("EQTLJointVec") != string::npos) {
        // //     mcmcSamples = new McmcSamples(parSet->label, chainLength, burnin, thin, parSet->size, "sparse");
        // //     if (writeBinPosterior) mcmcSamples->initBinFile(title);
        // //     // mcmcSamples->initTxtFile(title);
        // } 
        /////////////////////
        else if (parSet->label.find("Delta") != string::npos) {
            mcmcSamples = new McmcSamples(parSet->label, chainLength, burnin, thin, parSet->size, "sparse");
            if (writeBinPosterior) mcmcSamples->initBinFile(title);
        } else {
            mcmcSamples = new McmcSamples(parSet->label, chainLength, burnin, thin, parSet->size);
            if (writeTxtPosterior) mcmcSamples->initTxtFile(title);
        }
        mcmcSampleVec.push_back(mcmcSamples);
    }
    // Step 3. parameter with matrix values.
    for (unsigned i = 0; i < model.paramMatVec.size();++i){
        ParamMat *parMatSet = model.paramMatVec[i];
        McmcSamples *mcmcSamples;
        if (parMatSet->label.find("EQTLJointMat") != string::npos){
            mcmcSamples = new McmcSamples(parMatSet->label, chainLength, burnin, thin, parMatSet->ncol, "sparse");
            if (writeBinPosterior) mcmcSamples->initBinFile(title);
        } else if(parMatSet->label.find("SnpEffectMat") != string::npos){
            mcmcSamples = new McmcSamples(parMatSet->label, chainLength, burnin, thin, parMatSet->ncol, "sparse");
            if (writeBinPosterior) mcmcSamples->initBinFile(title);
        } else if (parMatSet->label.find("SigmaSqMats") != string::npos) {
            mcmcSamples = new McmcSamples(parMatSet->label, chainLength, burnin, thin, parMatSet->ncol, "sparse");
        }else {
            mcmcSamples = new McmcSamples(parMatSet->label, chainLength, burnin, thin, parMatSet->ncol, "sparse");
        }
        mcmcSampleVec.push_back(mcmcSamples);
    }
    if (writeTxtPosterior) initTxtFile(model.paramVec, title);
    return mcmcSampleVec;
}

void MCMC::collectSamples(const Model &model, vector<McmcSamples*> &mcmcSampleVec, const unsigned iteration, const bool writeBinPosterior, const bool writeTxtPosterior){
    unsigned i = 0;
    // Step 1. parameter with single value(double);
    // e.g., vare, pi 
    for (unsigned j=0; j<model.paramVec.size(); ++j) {
        McmcSamples *mcmcSamples = mcmcSampleVec[i++];
        Parameter *par = model.paramVec[j];
        mcmcSamples->getSample(iteration, par->value, writeTxtPosterior, out);
    }
    // Step 2. parameter with many values(vector): 
    // e.g., delta, beta
    for (unsigned j=0; j<model.paramSetVec.size(); ++j) {
        McmcSamples *mcmcSamples = mcmcSampleVec[i++];
        ParamSet *parSet = model.paramSetVec[j];
        mcmcSamples->getSample(iteration, parSet->values, writeBinPosterior, writeTxtPosterior);
    }
    // Step 3. parameter with matrix values.
    for(unsigned j= 0; j < model.paramMatVec.size(); ++j) {
        McmcSamples *mcmcSamples = mcmcSampleVec[i++];
        ParamMat * parMatSet = model.paramMatVec[j];
        mcmcSamples->getSample(iteration,parMatSet->values,writeBinPosterior, writeTxtPosterior);
    }
    out << endl;
    // LOGGER << endl;
}

void MCMC::printStatus(const vector<Parameter*> &paramToPrint, const unsigned thisIter, const unsigned outputFreq, const string &timeLeft){
    if (thisIter==outputFreq) {
        LOGGER << boost::format("%=10s ") % "Iter";
        for (unsigned i=0; i<paramToPrint.size(); ++i) {
            LOGGER << boost::format("%=12s ") % paramToPrint[i]->label;
        }
        LOGGER << boost::format("%=12s\n") % "TimeLeft";
    }
    LOGGER << boost::format("%=10s ") % thisIter;
    for (unsigned i=0; i<paramToPrint.size(); ++i) {
        Parameter *par = paramToPrint[i];
        if (par->label[0] == 'N'){
            LOGGER << boost::format("%=12.0f ") % par->value;
        }else if( paramToPrint[i]->value < 0.001 && paramToPrint[i]->value > 10e-15 ){
            LOGGER << boost::format("%=12.6e ") % paramToPrint[i]->value;
        }else if( paramToPrint[i]->value > 10e6 ){
            LOGGER << boost::format("%=12.6e ") % paramToPrint[i]->value;
        } else if((paramToPrint[i]->value <= 10e-15 )){ 
            LOGGER << boost::format("%=12.0f ") % 0;
        } else {
            LOGGER << boost::format("%=12.6f ") % paramToPrint[i]->value;
        }
    }
    LOGGER << boost::format("%=12s\n") % timeLeft;
    
    LOGGER.flush();
}

void MCMC::printSummary(const vector<Parameter*> &paramToPrint, const vector<McmcSamples*> &mcmcSampleVec, const string &filename){
    if (!paramToPrint.size()) return;
    ofstream out;
    out.open(filename.c_str());
    if (!out) {
        LOGGER.e(0," cannot open file " + filename);
    }

    LOGGER << "\nPosterior statistics from MCMC samples:\n\n";
    LOGGER << boost::format("%13s %-15s %-15s\n") %"Parameter" % "Mean" % "SD ";
    out << "Posterior statistics from MCMC samples:\n\n";
    out << boost::format("%13s %-15s %-15s\n") %"Parameter" % "Mean" % "SD ";
    for (unsigned i=0; i<paramToPrint.size(); ++i) {
        Parameter *par = paramToPrint[i];
        for (unsigned j=0; j<mcmcSampleVec.size(); ++j) {
            McmcSamples *mcmcSamples = mcmcSampleVec[j];
            if (mcmcSamples->label == par->label) {
                LOGGER << boost::format("%10s %2s %-12.20f %-12.20f\n")
                // LOGGER << boost::format("%10s %2s %-1% %-1%\n")
                % par->label
                % ""
                % mcmcSamples->mean()
                % mcmcSamples->sd();
                out << boost::format("%10s %2s %-12.29f %-12.20f\n")
                % par->label
                % ""
                % mcmcSamples->mean()
                % mcmcSamples->sd();
                break;
            }
        }
    }
    out.close();
}

void MCMC::printSetSummary(const vector<ParamSet*> &paramSetToPrint, const vector<McmcSamples*> &mcmcSampleVec, const string &filename){
    if (!paramSetToPrint.size()) return;
    ofstream out;
    out.open(filename.c_str());
    if (!out) {
        LOGGER.e(0," cannot open file " + filename);
    }
//    LOGGER << "\nPosterior statistics from MCMC samples:\n\n";
//    LOGGER << boost::format("%13s %-15s %-15s\n") %"" % "Mean" % "SD ";
//    out << "Posterior statistics from MCMC samples:\n\n";
//    out << boost::format("%13s %-15s %-15s\n") %"" % "Mean" % "SD ";
    for (unsigned i=0; i<paramSetToPrint.size(); ++i) {
        ParamSet *parset = paramSetToPrint[i];
        if (parset->label == "SnpAnnoMembershipDelta") continue;
        for (unsigned j=0; j<mcmcSampleVec.size(); ++j) {
            McmcSamples *mcmcSamples = mcmcSampleVec[j];
            if (mcmcSamples->label == parset->label) {
                for (unsigned col=0; col<parset->size; ++col) {
//                    LOGGER << boost::format("%20s %10s %2s %-15.6f %-15.6f\n")
//                    % parset->label
//                    % parset->header[col]
//                    % ""
//                    % mcmcSamples->posteriorMean[col]
//                    % sqrt(mcmcSamples->posteriorSqrMean[col]-mcmcSamples->posteriorMean[col]*mcmcSamples->posteriorMean[col]);
                    out << boost::format("%25s %20s %2s %-15.6f %-15.6f ")
                    % parset->label
                    % parset->header[col]
                    % ""
                    % mcmcSamples->posteriorMean[col]
                    % sqrt(mcmcSamples->posteriorSqrMean[col]-mcmcSamples->posteriorMean[col]*mcmcSamples->posteriorMean[col]);
                    Gadget::Tokenizer token;
                    token.getTokens(parset->label, "_");
                    double postprob = 0;
                    if (token.back() == "Enrichment") {
                        for (unsigned row=0; row<mcmcSamples->nrow; ++row) {
                            if (mcmcSamples->datMat(row, col) > 1) ++postprob;
                        }
                    } else {
                        for (unsigned row=0; row<mcmcSamples->nrow; ++row) {
                            if (mcmcSamples->datMat(row, col) > 0) ++postprob;
                        }
                    }
                    postprob /= double(mcmcSamples->nrow);
                    out << boost::format("%-15.6f\n") % postprob;
                }
                break;
            }
        }
    }
    out.close();
}

vector<McmcSamples*> MCMC::run(Model &model, const unsigned chainLength, const unsigned burnin, const unsigned thin, const bool print,
                               const unsigned outputFreq, const string &title, const bool writeBinPosterior, const bool writeTxtPosterior){
    if (print) {
        LOGGER << "MCMC launched ..." << endl;
        LOGGER << "  Chain length: " << chainLength << " iterations" << endl;
        LOGGER << "  Burn-in: " << burnin << " iterations" << endl << endl;
    }
    if (writeBinPosterior || writeTxtPosterior) {
        if (!Gadget::directoryExist(title + ".mcmcsamples")){
            Gadget::createDirectory(title + ".mcmcsamples");
            if (print) LOGGER << "  Created directory [" << title << ".mcmcsamples] to store MCMC samples.\n\n";
        }
    }
    // Step 1. initialize various parameters
    vector<McmcSamples*> mcmcSampleVec = initMcmcSamples(model, chainLength, burnin, thin, title, writeBinPosterior, writeTxtPosterior);
    Gadget::Timer timer;
    timer.setTime();
    // Step 2. use loops to find best values
    for (unsigned iteration=0; iteration<chainLength; ++iteration) {
        unsigned thisIter = iteration + 1;
        
        model.sampleUnknowns();
        collectSamples(model, mcmcSampleVec, iteration, writeBinPosterior, writeTxtPosterior);
        
        if (!(thisIter % outputFreq)) {
            timer.getTime();
            time_t timeToFinish = (chainLength-thisIter)*timer.getElapse()/thisIter; // remaining iterations multiplied by average time per iteration in seconds
            if (print) {
                printStatus(model.paramToPrint, thisIter, outputFreq, timer.format(timeToFinish));
            }
        }
    }
    // save the samples in the last iteration for potential continual run
    if (print) {
        LOGGER << "\nMCMC cycles completed." << endl;
        printSummary(model.paramToPrint, mcmcSampleVec, title + ".parRes");
        printSetSummary(model.paramSetToPrint, mcmcSampleVec, title + ".parSetRes");
        // printSnpAnnoMembership(model.paramSetToPrint, mcmcSampleVec, title + ".snpAnnoMembership");
    }
    return mcmcSampleVec;
}

void MCMC::convergeDiagGelmanRubin(const Model &model, vector<vector<McmcSamples *> > &mcmcSampleVecChain, const string &filename){
    if (!model.paramToPrint.size()) return;
    ofstream out;
    out.open((filename + ".parRes").c_str());
    LOGGER << "\nPosterior statistics from multiple chains:\n\n";
    LOGGER << boost::format("%13s %-15s %-15s %-12s\n") %"" % "Mean" % "SD " % "R_GelmanRubin ";
    out << "Posterior statistics from multiple chains:\n\n";
    out << boost::format("%13s %-15s %-15s %-12s\n") %"" % "Mean" % "SD " % "R_GelmanRubin ";
    long numChains = mcmcSampleVecChain.size();
    VectorXd meanVec(numChains);
    VectorXd varVec(numChains);
    for (unsigned i=0; i<model.paramToPrint.size(); ++i) {
        Parameter *par = model.paramToPrint[i];
        for (unsigned j=0; j<mcmcSampleVecChain[0].size(); ++j) {
            McmcSamples *mcmcSamples = mcmcSampleVecChain[0][j];
            if (mcmcSamples->label == par->label) {
                double nsample = mcmcSamples->nrow;
                for (unsigned m=0; m<numChains; ++m) {
                    mcmcSamples = mcmcSampleVecChain[m][j];
                    meanVec[m] = mcmcSamples->mean()[0];
                    varVec[m]  = mcmcSamples->sd()[0];
                    varVec[m] *= varVec[m];
                }
                double posteriorMean = meanVec.mean();
                double B = (meanVec.array() - posteriorMean).matrix().squaredNorm()*nsample/double(numChains-1);
                double W = varVec.mean();
                double posteriorVar = (nsample-1.0)*W/nsample + B/nsample;
                double R = sqrt(posteriorVar/W);
                
                LOGGER << boost::format("%10s %2s %-15.6f %-15.6f %-12.3f\n")
                % par->label
                % ""
                % posteriorMean
                % sqrt(posteriorVar)
                % R;
                out << boost::format("%10s %2s %-15.6f %-15.6f %-12.3f\n")
                % par->label
                % ""
                % posteriorMean
                % sqrt(posteriorVar)
                % R;
                break;
            }
        }
    }
    out.close();
    
    if (model.paramSetToPrint.size()) {
        ofstream out2;
        out2.open((filename + ".parSetRes").c_str());
        
        for (unsigned i=0; i<model.paramSetToPrint.size(); ++i) {
            ParamSet *parset = model.paramSetToPrint[i];
            for (unsigned j=0; j<mcmcSampleVecChain[0].size(); ++j) {
                McmcSamples *mcmcSamples = mcmcSampleVecChain[0][j];
                if (mcmcSamples->label == parset->label) {
                    MatrixXd meanMat(numChains, parset->size);
                    MatrixXd varMat(numChains, parset->size);
                    double nsample = mcmcSamples->nrow;
                    for (unsigned m=0; m<numChains; ++m) {
                        mcmcSamples = mcmcSampleVecChain[m][j];
                        meanMat.row(m) = mcmcSamples->mean();
                        varMat.row(m)  = mcmcSamples->sd();
                        varMat.row(m) *= varMat.row(m);
                    }
                    VectorXd posteriorMean = meanMat.colwise().mean();
                    VectorXd B = (meanMat.rowwise() - posteriorMean.transpose()).colwise().squaredNorm()*nsample/double(numChains-1);
                    VectorXd W = varMat.colwise().mean();
                    VectorXd posteriorVar = (nsample-1.0)*W/nsample + B/nsample;
                    VectorXd R = sqrt(posteriorVar.array()/W.array());
                    
                    for (unsigned col=0; col<parset->size; ++col) {
                        
                        out2 << boost::format("%25s %20s %2s %-15.6f %-15.6f ")
                        % parset->label
                        % parset->header[col]
                        % ""
                        % posteriorMean[col]
                        % sqrt(posteriorVar[col]);
                        Gadget::Tokenizer token;
                        token.getTokens(parset->label, "_");
                        double postprob = 0;
                        if (token.back() == "Enrichment") {
                            for (unsigned row=0; row<mcmcSamples->nrow; ++row) {
                                if (mcmcSamples->datMat(row, col) > 1) ++postprob;
                            }
                        } else {
                            for (unsigned row=0; row<mcmcSamples->nrow; ++row) {
                                if (mcmcSamples->datMat(row, col) > 0) ++postprob;
                            }
                        }
                        postprob /= double(mcmcSamples->nrow);
                        out2 << boost::format("%-15.6f %-12.3f\n") % postprob % R[col];
                    }
                    break;
                }
            }
        }
        
    }

}
