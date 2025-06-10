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

#include "Model.hpp"



void BayesC::FixedEffects::sampleFromFC(VectorXd &ycorr, const MatrixXd &X,
                                        const VectorXd &XPXdiag, const double vare){
    double rhs;
    for (unsigned i=0; i<size; ++i) {
        if (!XPXdiag[i]) continue;
        double oldSample = values[i];
        double rhs = X.col(i).dot(ycorr);
        rhs += XPXdiag[i]*oldSample;
        double invLhs = 1.0/XPXdiag[i];
        double bhat = invLhs*rhs;
        values[i] = Normal::sample(bhat, invLhs*vare);
        ycorr += X.col(i) * (oldSample - values[i]);
    }
}

void BayesC::RandomEffects::sampleFromFC(VectorXd &ycorr, const MatrixXd &W, const VectorXd &WPWdiag, const VectorXd &Rsqrt, const bool weightedRes, const double sigmaSqRand, const double vare, VectorXd &rhat){
    rhat.setZero(ycorr.size());
    double invVare = 1.0/vare;
    double invSigmaSqRand = 1.0/sigmaSqRand;
    double rhs = 0.0;
    ssq = 0.0;
    for (unsigned i=0; i<size; ++i) {
        if (!WPWdiag[i]) continue;
        double oldSample = values[i];
        double rhs = W.col(i).dot(ycorr) + WPWdiag[i]*oldSample;
        rhs *= invVare;
        double invLhs = 1.0/(WPWdiag[i]*invVare + invSigmaSqRand);
        double uhat = invLhs*rhs;
        values[i] = Normal::sample(uhat, invLhs);
        ssq += values[i]*values[i];
        if (weightedRes) rhat += W.col(i).cwiseProduct(Rsqrt) * values[i];
        else rhat  += W.col(i) * values[i];
        ycorr += W.col(i) * (oldSample - values[i]);
    }
}

void BayesC::VarRandomEffects::sampleFromFC(const double randEffSumSq, const unsigned int numRandEff){
    double dfTilde = df + numRandEff;
    double scaleTilde = randEffSumSq + df*scale;
    value = InvChiSq::sample(dfTilde, scaleTilde);    
}

void BayesC::SnpEffects::sampleFromFC(VectorXd &ycorr, const MatrixXd &Z, const VectorXd &ZPZdiag, const VectorXd &Rsqrt, const bool weightedRes,
const double sigmaSq, const double pi, const double vare, VectorXd &ghat){
    if (algorithm == gibbs) {
        gibbsSampler(ycorr, Z, ZPZdiag, Rsqrt, weightedRes, sigmaSq, pi, vare, ghat);
    } else if (algorithm == hmc) {
        hmcSampler(ycorr, Z, ZPZdiag, sigmaSq, pi, vare, ghat);
    }
}

void BayesC::SnpEffects::gibbsSampler(VectorXd &ycorr, const MatrixXd &Z, const VectorXd &ZPZdiag, const VectorXd &Rsqrt, const bool weightedRes,
                                      const double sigmaSq, const double pi, const double vare, VectorXd &ghat){
    sumSq = 0.0;
    numNonZeros = 0;
    
    pip.setZero(size);
    ghat.setZero(ycorr.size());
    
    double oldSample;
    double rhs, invLhs, uhat;
    double logDelta0, logDelta1, probDelta1;
    double logPi = log(pi);
    double logPiComp = log(1.0-pi);
    double logSigmaSq = log(sigmaSq);
    double invVare = 1.0/vare;
    double invSigmaSq = 1.0/sigmaSq;
    
    for (unsigned i=0; i<size; ++i) {
        oldSample = values[i];
        rhs = Z.col(i).dot(ycorr);
        rhs += ZPZdiag[i]*oldSample;
        rhs *= invVare;
        invLhs = 1.0/(ZPZdiag[i]*invVare + invSigmaSq);
        uhat = invLhs*rhs;
        logDelta1 = 0.5*(logf(invLhs) - logSigmaSq + uhat*rhs) + logPi;
        //logDelta1 = rhs*oldSample - 0.5*ZPZdiag[i]*oldSample*oldSample/vare + logPiComp;
        logDelta0 = logPiComp;
        probDelta1 = 1.0/(1.0 + expf(logDelta0-logDelta1));
        pip[i] = probDelta1;
        
        if (bernoulli.sample(probDelta1)) {
            values[i] = normal.sample(uhat, invLhs);
            ycorr += Z.col(i) * (oldSample - values[i]);
            if (weightedRes) ghat += Z.col(i).cwiseProduct(Rsqrt) * values[i];
            else ghat  += Z.col(i) * values[i];
            sumSq += values[i]*values[i];
            ++numNonZeros;
        } else {
            if (oldSample) ycorr += Z.col(i) * oldSample;
            values[i] = 0.0;
        }
    }
}

void BayesC::SnpEffects::hmcSampler(VectorXd &ycorr, const MatrixXd &Z, const VectorXd &ZPZdiag, const double sigmaSq, const double pi, const double vare, VectorXd &ghat){
    // Hamiltonian Monte Carlo
    // Only BayesC0 model available
    
    double stepSize = 0.1;
    unsigned numSteps = 10;
    
    ycorr += Z*values;
    
    static MatrixXd ZPZ;
    if (cnt==0) ZPZ = Z.transpose()*Z;
    VectorXd ypZ = ycorr.transpose()*Z;
    
    VectorXd curr = values;
    
    ArrayXd curr_p(size);
    for (unsigned i=0; i<size; ++i) {
        curr_p[i] = Stat::snorm();
    }
    
    VectorXd cand = curr;
    // Make a half step for momentum at the beginning
    ArrayXd cand_p = curr_p - 0.5*stepSize * gradientU(curr, ZPZ, ypZ, sigmaSq, vare);
    
    for (unsigned i=0; i<numSteps; ++i) {
        cand.array() += stepSize * cand_p;
        if (i < numSteps-1) {
            cand_p -= stepSize * gradientU(cand, ZPZ, ypZ, sigmaSq, vare);
        } else {
            cand_p -= 0.5*stepSize * gradientU(cand, ZPZ, ypZ, sigmaSq, vare);
        }
    }
    
    double curr_H = computeU(curr, ZPZ, ypZ, sigmaSq, vare) + 0.5*curr_p.matrix().squaredNorm();
    double cand_H = computeU(cand, ZPZ, ypZ, sigmaSq, vare) + 0.5*cand_p.matrix().squaredNorm();
    
    if (Stat::ranf() < exp(curr_H-cand_H)) {  // accept
        values = cand;
        ghat = Z*values;
        ++mhr;
    }
    
    if (!(++cnt % 100)) {
        double ar = mhr/double(cnt);
        if      (ar < 0.5) LOGGER << "Warning: acceptance rate for SNP effects is too low "  << ar << endl;
        else if (ar > 0.9) LOGGER << "Warning: acceptance rate for SNP effects is too high " << ar << endl;
    }
    
    numNonZeros = size;
    sumSq = values.squaredNorm();
    
    ycorr -= Z*values;
}

ArrayXd BayesC::SnpEffects::gradientU(const VectorXd &alpha, const MatrixXd &ZPZ, const VectorXd &ypZ, const double sigmaSq, const double vare){
    return 1.0/vare*(ZPZ*alpha - ypZ) + 1/sigmaSq*alpha;
}

double BayesC::SnpEffects::computeU(const VectorXd &alpha, const MatrixXd &ZPZ, const VectorXd &ypZ, const double sigmaSq, const double vare){
    return 0.5/vare*(alpha.transpose()*ZPZ*alpha + vare/sigmaSq*alpha.squaredNorm() - 2.0*ypZ.dot(alpha));
}

void BayesC::SnpEffects::sampleFromFC_omp(VectorXd &ycorr, const MatrixXd &Z, const VectorXd &ZPZdiag,
                                          const double sigmaSq, const double pi, const double vare, VectorXd &ghat){
    // speed-enhanced single site Gibbs sampling due to the use of parallel computing on SNPs with zero effect
    
    unsigned blockSize = 1; //omp_get_max_threads();
    //LOGGER << blockSize << endl;
    
    sumSq = 0.0;
    numNonZeros = 0;
    
    ghat.setZero(ycorr.size());
    
    double oldSample;
    double logPi = log(pi);
    double logPiComp = log(1.0-pi);
    double logSigmaSq = log(sigmaSq);
    double invVare = 1.0/vare;
    double invSigmaSq = 1.0/sigmaSq;
    
    vector<int> deltaVec(blockSize);
    vector <double>  invLhsVec(blockSize);
    vector <double>  uhatVec(blockSize);
    
    unsigned blocki;
    unsigned i, j;
    bool breakFlag;
    
    for (i=0; i<size; ) {
        
        if (blockSize + i < size) {
            blocki = blockSize;
        } else {
            blocki = size - i;
            deltaVec.resize(blocki);
            invLhsVec.resize(blocki);
            uhatVec.resize(blocki);
        }
        
        #pragma omp parallel for
        for (j=0; j<blocki; ++j) {
            double rhsj = (Z.col(i+j).dot(ycorr) + ZPZdiag[i+j]*values[i+j])*invVare;
            invLhsVec[j] = 1.0/(ZPZdiag[i+j]*invVare + invSigmaSq);
            uhatVec[j] = invLhsVec[j]*rhsj;
            double logDelta0minusDelta1j = logPiComp - (0.5f*(logf(invLhsVec[j]) - logSigmaSq + uhatVec[j]*rhsj) + logPi);
            deltaVec[j] = bernoulli.sample(1.0/(1.0 + expf(logDelta0minusDelta1j)));
        }
        
        breakFlag = false;
        for (j=0; j<blocki; ++j) {
            if (values[i+j] || deltaVec[j]) {   // need to update ycorr for the first snp who is in the model at either last or this iteration
                i += j;
                breakFlag = true;
                break;
            }
        }
        
        if (breakFlag) {
            oldSample = values[i];
            if (deltaVec[j]) {
                values[i] = normal.sample(uhatVec[j], invLhsVec[j]);
                ycorr += Z.col(i) * (oldSample - values[i]);
                ghat  += Z.col(i) * values[i];
                sumSq += values[i]*values[i];
                ++numNonZeros;
            } else {
                if (oldSample) ycorr += Z.col(i) * oldSample;
                values[i] = 0.0;
            }
            ++i;
        }
        else {
            i += blocki;
        }
    }
}

void BayesC::VarEffects::sampleFromFC(const double snpEffSumSq, const unsigned numSnpEff){
    double dfTilde = df + numSnpEff;
    double scaleTilde = snpEffSumSq + df*scale;
    value = InvChiSq::sample(dfTilde, scaleTilde);
    //LOGGER << "snpEffSumSq " << snpEffSumSq << " scale " << scale << " scaleTilde " << scaleTilde << " dfTilde " << dfTilde << " value " << value << endl;
}

void BayesC::VarEffects::sampleFromPrior(){
    value = InvChiSq::sample(df, scale);
}

void BayesC::VarEffects::computeScale(const double varg, const VectorXd &snp2pq, const double pi){
    if (noscale)
        scale = (df-2)/df * varg/(snp2pq.sum()*pi);
    else
        scale = (df-2)/df * varg/(snp2pq.size()*pi);
}

void BayesC::VarEffects::computeScale(const double varg, const double sum2pq){
        scale = (df-2)/df * varg/sum2pq;
}

void BayesC::VarEffects::compute(const double snpEffSumSq, const double numSnpEff){
    if (numSnpEff) value = snpEffSumSq/numSnpEff;
}

void BayesC::ScaleVar::sampleFromFC(const double sigmaSq, const double df, double &scaleVar){
    double shapeTilde = shape + 0.5*df;
    double scaleTilde = 1.0/(1.0/scale + 0.5*df/sigmaSq);
    value = Gamma::sample(shapeTilde, scaleTilde);
    scaleVar = value;
}

void BayesC::Pi::sampleFromFC(const unsigned numSnps, const unsigned numSnpEff){
    double alphaTilde = numSnpEff + alpha;
    double betaTilde  = numSnps - numSnpEff + beta;
    value = Beta::sample(alphaTilde, betaTilde);
}

void BayesC::Pi::sampleFromPrior(){
    value = Beta::sample(alpha, beta);
}

void BayesC::Pi::compute(const double numSnps, const double numSnpEff){
    value = numSnpEff/numSnps;
}

void BayesC::ResidualVar::sampleFromFC(VectorXd &ycorr){
    double sse = ycorr.squaredNorm();
    double dfTilde = df + nobs;
    double scaleTilde = sse + df*scale;
    value = InvChiSq::sample(dfTilde, scaleTilde);
}

void BayesC::GenotypicVar::compute(const VectorXd &ghat){
    //value = Gadget::calcVariance(ghat);
    double sum = ghat.sum();
    double ssq = ghat.squaredNorm();
    unsigned size = (unsigned)ghat.size();
    double mean = sum/size;
    value = ssq/size - mean*mean;
}

void BayesC::RandomVar::compute(const VectorXd &rhat){
    //value = Gadget::calcVariance(ghat);
    double sum = rhat.sum();
    double ssq = rhat.squaredNorm();
    unsigned size = (unsigned)rhat.size();
    double mean = sum/size;
    value = ssq/size - mean*mean;
}

void BayesC::Rounding::computeYcorr(const VectorXd &y, const MatrixXd &X, const MatrixXd &W, const MatrixXd &Z,
                                    const VectorXd &fixedEffects, const VectorXd &randomEffects, const VectorXd &snpEffects,
                                    VectorXd &ycorr){
    if (count++ % 100) return;
    VectorXd oldYcorr = ycorr;
    ycorr = y - X*fixedEffects;
    if (randomEffects.size()) ycorr -= W*randomEffects;
    for (unsigned i=0; i<snpEffects.size(); ++i) {
        if (snpEffects[i]) ycorr -= Z.col(i)*snpEffects[i];
    }
    double ss = (ycorr - oldYcorr).squaredNorm();
    value = sqrt(ss);
}

void BayesC::sampleUnknowns(){
    fixedEffects.sampleFromFC(ycorr, data.X, data.XPXdiag, vare.value);
    if (data.numRandomEffects) {
        randomEffects.sampleFromFC(ycorr, data.W, data.WPWdiag, data.Rsqrt, data.weightedRes, sigmaSqRand.value, vare.value, rhat);
        sigmaSqRand.sampleFromFC(randomEffects.ssq, data.numRandomEffects);
        varRand.compute(rhat);
    }
    unsigned cnt=0;
   do {
        snpEffects.sampleFromFC(ycorr, data.Z, data.ZPZdiag, data.Rsqrt, data.weightedRes, sigmaSq.value, pi.value, vare.value, ghat);
       if (++cnt == 100) LOGGER.e(0," Zero SNP effect in the model for 100 cycles of sampling");
   } while (snpEffects.numNonZeros == 0);
    snpPip.getValues(snpEffects.pip);
    sigmaSq.sampleFromFC(snpEffects.sumSq, snpEffects.numNonZeros);
    //scale.sampleFromFC(sigmaSq.value, sigmaSq.df, sigmaSq.scale);
    if (estimatePi) pi.sampleFromFC(snpEffects.size, snpEffects.numNonZeros);
    vare.sampleFromFC(ycorr);
    
    varg.compute(ghat);
    hsq.compute(varg.value, vare.value);
    
    rounding.computeYcorr(data.y, data.X, data.W, data.Z, fixedEffects.values, randomEffects.values, snpEffects.values, ycorr);
    nnzSnp.getValue(snpEffects.numNonZeros);
}

void BayesC::sampleStartVal(){
    sigmaSq.sampleFromPrior();
    if (estimatePi) pi.sampleFromPrior();
        LOGGER << "  Starting value for " << sigmaSq.label << ": " << sigmaSq.value << endl;
    if (estimatePi) LOGGER << "  Starting value for " << pi.label << ": " << pi.value << endl;
        LOGGER << endl;
}

// ----------------------------------------------------------------------------------------
// Bayes R
// ----------------------------------------------------------------------------------------

void BayesR::ProbMixComps::sampleFromFC(const VectorXd &snpStore) {
	VectorXd dirx;
	dirx = snpStore + alphaVec;
    values = Dirichlet::sample(ndist, dirx);
    for (unsigned i=0; i<ndist; ++i) {
      (*this)[i]->value=values[i];  
    }
}

void BayesR::NumSnpMixComps::getValues(const VectorXd &snpStore) {
    values = snpStore;
    for (unsigned i=0; i<ndist; ++i) {
        (*this)[i]->value=values[i];
    }
}

void BayesR::VgMixComps::compute(const VectorXd &snpEffects, const MatrixXd &Z, const vector<vector<unsigned> > snpset, const double varg) {
    values.setZero(ndist);
    long nobs = Z.rows();
//    for (unsigned k=0; k<ndist; ++k) {
//        if (k!=zeroIdx && k!=minIdx) {
//            long numSnps = snpset[k].size();
//            unsigned idx;
//            VectorXd ghat;
//            ghat.setZero(nobs);
//            for (unsigned i=0; i<numSnps; ++i) {
//                idx = snpset[k][i];
//                ghat += snpEffects[idx]*Z.col(idx);
//            }
//            (*this)[k]->value = values[k] = Gadget::calcVariance(ghat)/varg;
//        }
//    }
//    double sum = values.sum();
//    (*this)[minIdx]->value = values[minIdx] = 1.0 - sum;

    for (unsigned k=0; k<ndist; ++k) {
        if (k!=zeroIdx) {
            long numSnps = snpset[k].size();
            unsigned idx;
            VectorXd ghat;
            ghat.setZero(nobs);
            for (unsigned i=0; i<numSnps; ++i) {
                idx = snpset[k][i];
                ghat += snpEffects[idx]*Z.col(idx);
            }
            values[k] = Gadget::calcVariance(ghat);
        }
    }
    double sum = values.sum();
    for (unsigned k=0; k<ndist; ++k) {
        (*this)[k]->value = values[k] = values[k]/sum;
    }

}

void BayesR::SnpEffects::sampleFromFC(VectorXd &ycorr, const MatrixXd &Z, const VectorXd &ZPZdiag, const VectorXd &Rsqrt, const bool weightedRes,
                                      const double sigmaSq, const VectorXd &pis, const VectorXd &gamma,
                                      const double vare, VectorXd &ghat, VectorXd &snpStore,
                                      const double varg, const bool originalModel){
    sumSq = 0.0;
    numNonZeros = 0;
    
    ghat.setZero(ycorr.size());
    double oldSample;
    double rhs;
    // -----------------------------------------
    // Initialise the parameters in MCMC sampler
    // -----------------------------------------
    // ----------------
    // Bayes R specific
    // ----------------
    int ndist, indistflag;
    double v1,  b_ls, ssculm, r;
    VectorXd gp, ll, ll2, pll, snpindist, var_b_ls;
    ndist = pis.size();
    snpStore.setZero(pis.size());
    pll.setZero(pis.size());
    // --------------------------------------------------------------------------------
    // Scale the variances in each of the normal distributions by the genetic variance
    // and initialise the class membership probabilities
    // --------------------------------------------------------------------------------
    if (originalModel)
        gp = gamma * 0.01 * varg;
    else
        gp = gamma * sigmaSq;
//    cout << varg << " " << gp.transpose() << endl;
    snpset.resize(ndist);
    for (unsigned k=0; k<ndist; ++k) {
        snpset[k].resize(0);
    }
    
    for (unsigned i=0; i<size; ++i) {
        // ------------------------------
        // Derived Bayes R implementation
        // ------------------------------
        // ----------------------------------------------------
        // Add back the content for the corrected rhs for SNP k
        // ----------------------------------------------------
        rhs = Z.col(i).dot(ycorr);
        oldSample = values[i];
        rhs += ZPZdiag[i] * oldSample;
        // ------------------------------------------------------
        // Calculate the beta least squares updates and variances
        // ------------------------------------------------------
        b_ls = rhs / ZPZdiag[i];
        var_b_ls = gp.array() + vare / ZPZdiag[i];
        // ------------------------------------------------------
        // Calculate the likelihoods for each distribution
        // ------------------------------------------------------
        // ll  = (-1.0 / 2.0) * var_b_ls.array().log()  - (b_ls * b_ls)  / (2 * var_b_ls.array());
        ll = (-1.0 / 2.0) * var_b_ls.array().log()  - (b_ls * b_ls)  / (2 * var_b_ls.array()) + pis.array().log();
        // --------------------------------------------------------------
        // Calculate probability that snp is in each of the distributions
        // in this iteration
        // --------------------------------------------------------------
        // pll = (ll.array().exp().cwiseProduct(pis.array())) / ((ll.array().exp()).cwiseProduct(pis.array())).sum();
        for (unsigned k=0; k<pis.size(); ++k) {
            pll[k] = 1.0 / (exp(ll.array() - ll[k])).sum();
        }
        // --------------------------------------------------------------
        // Sample the group based on the calculated probabilities
        // --------------------------------------------------------------
        ssculm = 0.0;
        r = Stat::ranf();
        indistflag = 1;
        for (int kk = 0; kk < ndist; kk++)
        {
            ssculm += pll(kk);
            if (r < ssculm)
            {
                indistflag = kk + 1;
                snpStore(kk) = snpStore(kk) + 1; 
                break;
            }
        }
        snpset[indistflag-1].push_back(i);
        // --------------------------------------------------------------
        // Sample the effect given the group and adjust the rhs
        // --------------------------------------------------------------
        if (indistflag != 1)
        {
            v1 = ZPZdiag[i] + vare / gp((indistflag - 1));
            values[i] = normal.sample(rhs / v1, vare / v1);
            ycorr += Z.col(i) * (oldSample - values[i]);
            if (weightedRes) ghat += Z.col(i).cwiseProduct(Rsqrt) * values[i];
            else ghat  += Z.col(i) * values[i];
            sumSq += (values[i] * values[i]) / gamma[indistflag - 1];
            ++numNonZeros;
        } else {
            if (oldSample) ycorr += Z.col(i) * oldSample;
            values[i] = 0.0;
        }
    }
}

void BayesR::VarEffects::computeScale(const double varg, const VectorXd &snp2pq, const VectorXd &gamma, const VectorXd &pi){
    if (noscale)
        scale = (df-2)/df * varg/(snp2pq.sum()*gamma.dot(pi));
    else
        scale = (df-2)/df * varg/(snp2pq.size()*gamma.dot(pi));
}

void BayesR::sampleUnknowns(){
    fixedEffects.sampleFromFC(ycorr, data.X, data.XPXdiag, vare.value);
    if (data.numRandomEffects) {
        randomEffects.sampleFromFC(ycorr, data.W, data.WPWdiag, data.Rsqrt, data.weightedRes, sigmaSqRand.value, vare.value, rhat);
        sigmaSqRand.sampleFromFC(randomEffects.ssq, data.numRandomEffects);
        varRand.compute(rhat);
    }
    unsigned cnt=0;
    do {
        snpEffects.sampleFromFC(ycorr, data.Z, data.ZPZdiag, data.Rsqrt, data.weightedRes, sigmaSq.value, Pis.values, gamma.values, vare.value, ghat, snpStore, varg.value, originalModel);
        if (++cnt == 100) LOGGER.e(0,"Error: Zero SNP effect in the model for 100 cycles of sampling");
    } while (snpEffects.numNonZeros == 0);  
    sigmaSq.sampleFromFC(snpEffects.sumSq, snpEffects.numNonZeros);
    vare.sampleFromFC(ycorr);
    Pis.sampleFromFC(snpStore);
    numSnps.getValues(snpStore);
    varg.compute(ghat);
    hsq.compute(varg.value, vare.value);
    if (originalModel) Vgs.compute(snpEffects.values, data.Z, snpEffects.snpset, varg.value);
    rounding.computeYcorr(data.y, data.X, data.W, data.Z, fixedEffects.values, randomEffects.values, snpEffects.values, ycorr);
    nnzSnp.getValue(snpEffects.numNonZeros);
}


void ApproxBayesC::FixedEffects::sampleFromFC(const MatrixXd &XPX, const VectorXd &XPXdiag,
                                              const MatrixXd &ZPX, const VectorXd &XPy,
                                              const VectorXd &snpEffects, const double vare,
                                              VectorXd &rcorr){
    for (unsigned i=0; i<size; ++i) {
        double oldSample = values[i];
        double XPZa = ZPX.col(i).dot(snpEffects);
        double rhs = XPy[i] - XPZa - XPX.row(i).dot(values) + XPXdiag[i]*values[i];
        double invLhs = 1.0/XPXdiag[i];
        double bhat = invLhs*rhs;
        values[i] = Normal::sample(bhat, invLhs*vare);
        //rcorr += ZPX.col(i) * (oldSample - values[i]);
    }

}

void ApproxBayesC::SnpEffects::sampleFromFC(VectorXd &rcorr,const vector<SparseVector <double>  > &ZPZ, const VectorXd &ZPZdiag, const VectorXd &ZPy,
                                            const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                                            const VectorXd &se, const VectorXd &tss, VectorXd &varei, const VectorXd &n, const VectorXd &snp2pq, const VectorXd &LDsamplVar,
                                            const double sigmaSq, const double pi, const double vare, const double varg, const double ps, const double overdispersion){
    
    long numChr = chromInfoVec.size();

    double ssq[numChr], s2pq[numChr], nnz[numChr];
    memset(ssq,0,sizeof(double)*numChr);
    memset(s2pq,0,sizeof(double)*numChr);
    memset(nnz,0, sizeof(double)*numChr);

//    for (unsigned chr=0; chr<numChr; ++chr) {
//        ChromInfo *chromInfo = chromInfoVec[chr];
//        unsigned chrStart = chromInfo->startSnpIdx;
//        unsigned chrEnd   = chromInfo->endSnpIdx;
//        if (iter==0) {
//            LOGGER << "chr " << chr+1 << " start " << chrStart << " end " << chrEnd << endl;
//        }
//    }
//    if (iter==0) LOGGER << endl; 

    double *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking precents the writing from other threads

    vector <double>  urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }
    
#pragma omp parallel for
    for (unsigned chr=0; chr<numChr; ++chr) {
        //LOGGER << " thread " << omp_get_thread_num() << " chr " << chr << endl;
        
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned windEnd, j;
        
        double oldSample;
        double rhs, invLhs, uhat;
        double logDelta0, logDelta1, probDelta1;
        double logPi = log(pi);
        double logPiComp = log(1.0-pi);
        double logSigmaSq = log(sigmaSq);
        double invSigmaSq = 1.0/sigmaSq;
        double varei;
        
        for (unsigned i=chrStart; i<=chrEnd; ++i) {

            oldSample = valuesPtr[i];

            if (leaveout[i]) {
                probDelta1 = 0;
            }
            else {

            //LOGGER << i << " " << chrStart << " " << chrEnd << " " << ssq << " " << nnz << endl;
//            if (!(iter % 100)) {
//                //double varei = (sse[i] - values.segment(windStart[i], windSize[i]).dot(ZPy.segment(windStart[i], windSize[i]) + rcorr.segment(windStart[i], windSize[i])))/n[i];
//                windEnd = windStart[i] + windSize[i];
//                varei[i] = tss[i];
//                for (j=windStart[i]; j<windEnd; ++j) {
//                    if (valuesPtr[j]) varei[i] -= valuesPtr[j]*(ZPy[j] + rcorr[j]);
//                }
//                varei[i] /= n[i];
//            }
                
            varei = LDsamplVar[i]*varg + vare + ps + overdispersion;
                
            rhs = rcorr[i] + ZPZdiag[i]*oldSample;
            rhs /= varei;
            invLhs = 1.0/(ZPZdiag[i]/varei + invSigmaSq);
            uhat = invLhs*rhs;
            logDelta1 = 0.5*(logf(invLhs) - logSigmaSq + uhat*rhs) + logPi;
            logDelta0 = logPiComp;
            probDelta1 = 1.0/(1.0 + expf(logDelta0-logDelta1));
                
            }
            
//            LOGGER << i << " " << rhs << " " << invLhs << " " << uhat << " " << probDelta1 << " " << sigmaSq << endl;
//            int tmp;
//            cin >> tmp;

//            if (bernoulli.sample(probDelta1)) {
            if (urnd[i] < probDelta1) {
//                valuesPtr[i] = normal.sample(uhat, invLhs);
                valuesPtr[i] = uhat + nrnd[i]*sqrtf(invLhs);
//                rcorr.segment(windStart[i], windSize[i]) += ZPZ[i]*(oldSample - values[i]);
                double sampleDiff = oldSample - valuesPtr[i];
                for (SparseVector <double> ::InnerIterator it(ZPZ[i]); it; ++it) {
                    rcorr[it.index()] += it.value() * sampleDiff;
                }
                ssq[chr]  += valuesPtr[i]*valuesPtr[i];
                s2pq[chr] += snp2pq[i];
                ++nnz[chr];
            } else {
                if (oldSample) {
//                    rcorr.segment(windStart[i], windSize[i]) += ZPZ[i]*oldSample;
                    for (SparseVector <double> ::InnerIterator it(ZPZ[i]); it; ++it) {
                        rcorr[it.index()] += it.value() * oldSample;
                    }
                }
                valuesPtr[i] = 0.0;
            }
        }
    }
    
    //LOGGER << ssq << " " << nnz << endl;

    sumSq = 0.0;
    sum2pq = 0.0;
    numNonZeros = 0;
    nnzPerChr.setZero(numChr);
    for (unsigned i=0; i<numChr; ++i) {
        sumSq += ssq[i];
        sum2pq += s2pq[i];
        numNonZeros += nnz[i];
        nnzPerChr[i] = nnz[i];
    }

    values = VectorXd::Map(valuesPtr, size);
}

void ApproxBayesC::SnpEffects::sampleFromFC(VectorXd &rcorr,const vector<VectorXd> &ZPZ, const VectorXd &ZPZdiag, const VectorXd &ZPy,
                                            const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                                            const VectorXd &se, const VectorXd &tss, VectorXd &varei, const VectorXd &n, const VectorXd &snp2pq, const VectorXd &LDsamplVar,
                                            const double sigmaSq, const double pi, const double vare, const double varg, const double ps, const double overdispersion){
    
    long numChr = chromInfoVec.size();
    
    double ssq[numChr], nnz[numChr], s2pq[numChr];
    memset(ssq,0,sizeof(double)*numChr);
    memset(nnz,0,sizeof(double)*numChr);
    memset(s2pq,0,sizeof(double)*numChr);

    double *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking precents the writing from other threads
  
    vector <double>  urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }
    
#pragma omp parallel for
    for (unsigned chr=0; chr<numChr; ++chr) {
        //LOGGER << " thread " << omp_get_thread_num() << " chr " << chr << endl;
        
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned windEnd, j;
        
        double oldSample;
        double rhs, invLhs, uhat;
        double logDelta0, logDelta1, probDelta1;
        double logPi = log(pi);
        double logPiComp = log(1.0-pi);
        double logSigmaSq = log(sigmaSq);
        double invSigmaSq = 1.0/sigmaSq;
        double varei;
        
        for (unsigned i=chrStart; i<=chrEnd; ++i) {

            oldSample = valuesPtr[i];

            if (leaveout[i]) {
                probDelta1 = 0;
            }
            else {
            
            //LOGGER << i << " " << chrStart << " " << chrEnd << " " << ssq << " " << nnz << endl;
            // LOGGER << "Varei i " << varei[i] << endl;
//            if (!(iter % 100)) {
//                //double varei = (sse[i] - values.segment(windStart[i], windSize[i]).dot(ZPy.segment(windStart[i], windSize[i]) + rcorr.segment(windStart[i], windSize[i])))/n[i];
//                windEnd = windStart[i] + windSize[i];
//                varei[i] = tss[i];
//                for (j=windStart[i]; j<windEnd; ++j) {
//                    if (valuesPtr[j]) varei[i] -= valuesPtr[j]*(ZPy[j] + rcorr[j]);
//                }
//                varei[i] /= n[i];
//                // LOGGER << "Varei 100 " << varei[i] << endl;
//            }
            // varei[i] = tss[i] / n[i];
            //            varei = se[i]*se[i]*ZPZdiag[i];
            
            varei = LDsamplVar[i]*varg + vare + ps + overdispersion;

            rhs = rcorr[i] + ZPZdiag[i]*oldSample;
            rhs /= varei;
            invLhs = 1.0/(ZPZdiag[i]/varei + invSigmaSq);
            uhat = invLhs*rhs;
            logDelta1 = 0.5*(logf(invLhs) - logSigmaSq + uhat*rhs) + logPi;
            logDelta0 = logPiComp;
            probDelta1 = 1.0/(1.0 + expf(logDelta0-logDelta1));
            //LOGGER << rhs << " " << invLhs << " " << logDelta1 << " " << logSigmaSq << " " << sigmaSq << endl;
                
            }
            
//            if (bernoulli.sample(probDelta1)) {
            if (urnd[i] < probDelta1) {
//                valuesPtr[i] = normal.sample(uhat, invLhs);
                valuesPtr[i] = uhat + nrnd[i]*sqrtf(invLhs);
                rcorr.segment(windStart[i], windSize[i]) += ZPZ[i]*(oldSample - valuesPtr[i]);
                ssq[chr] += valuesPtr[i]*valuesPtr[i];
                s2pq[chr] += snp2pq[i];
                ++nnz[chr];
            } else {
                if (oldSample) rcorr.segment(windStart[i], windSize[i]) += ZPZ[i]*oldSample;
                valuesPtr[i] = 0.0;
            }
        }
    }
    // LOGGER << "Varei 1 max" << varei.maxCoeff() << endl;
    //LOGGER << ssq << " " << nnz << endl;
    
    sumSq = 0.0;
    sum2pq = 0.0;
    numNonZeros = 0.0;
    nnzPerChr.setZero(numChr);
    for (unsigned i=0; i<numChr; ++i) {
        sumSq += ssq[i];
        sum2pq += s2pq[i];
        numNonZeros += nnz[i];
        nnzPerChr[i] = nnz[i];
    }
    
    values = VectorXd::Map(valuesPtr, size);
}

void ApproxBayesC::SnpEffects::hmcSampler(VectorXd &rcorr, const VectorXd &ZPy, const vector<VectorXd> &ZPZ,
                                            const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                                            const double sigmaSq, const double pi, const double vare){
    
    double stepSize = 0.001;
    unsigned numSteps = 1;
    
    
    //#pragma omp parallel for   // this multi-thread may not work due to vector locking when write to the vector
    for (unsigned chr=0; chr<chromInfoVec.size(); ++chr) {
        //LOGGER << " thread " << omp_get_thread_num() << " chr " << chr << endl;
        
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned chrSize  = chromInfo->size;
        
        VectorXd chrZPy = ZPy.segment(chrStart, chrSize);
        VectorXi chrWindStart = windStart.segment(chrStart, chrSize);
        VectorXi chrWindSize = windSize.segment(chrStart, chrSize);
        chrWindStart.array() -= chrStart;
        

        VectorXd delta;
        delta.setZero(chrSize);
        for (unsigned i=chrStart, j=0; i<=chrEnd; ++i) {
            if (values[i]) {
                delta[j++] = 1;
            }
        }
        
        
        VectorXd curr = values.segment(chrStart, chrSize);
        VectorXd curr_p(chrSize);
        
        for (unsigned i=0; i<chrSize; ++i) {
            curr_p[i] = Stat::snorm();
        }
        
        VectorXd cand = curr.cwiseProduct(delta);
        // Make a half step for momentum at the beginning
        VectorXd rc = chrZPy;
        VectorXd cand_p = curr_p.cwiseProduct(delta) - 0.5*stepSize * gradientU(curr, rc, chrZPy, ZPZ, chrWindStart, chrWindSize, chrStart, chrSize, sigmaSq, vare).cwiseProduct(delta);
        
        for (unsigned i=0; i<numSteps; ++i) {
            cand += stepSize * cand_p.cwiseProduct(delta);
            if (i < numSteps-1) {
                cand_p -= stepSize * gradientU(cand, rc, chrZPy, ZPZ, chrWindStart, chrWindSize, chrStart, chrSize, sigmaSq, vare).cwiseProduct(delta);
            } else {
                cand_p -= 0.5* stepSize * gradientU(cand, rc, chrZPy, ZPZ, chrWindStart, chrWindSize, chrStart, chrSize, sigmaSq, vare).cwiseProduct(delta);
            }
        }
        
        double curr_H = computeU(curr, rcorr.segment(chrStart, chrSize), chrZPy, sigmaSq, vare) + 0.5*curr_p.squaredNorm();
        double cand_H = computeU(cand, rc, chrZPy, sigmaSq, vare) + 0.5*cand_p.squaredNorm();
        
        if (Stat::ranf() < exp(curr_H-cand_H)) {  // accept
            values.segment(chrStart, chrSize) = cand;
            rcorr.segment(chrStart, chrSize) = rc;
            ++mhr;
        }
    }
    
    sumSq = values.squaredNorm();
    //numNonZeros = size;
    
    for (unsigned i=0; i<size; ++i) {
        if(values[i]) ++numNonZeros;
    }
    //LOGGER << sumSq << " " << nnz << " " << numNonZeros << endl;
    
    //LOGGER << values.head(10).transpose() << endl;
    
//    if (!(++cnt % 100) && myMPI::rank==0) {
//        double ar = mhr/double(cnt*22);
//        if      (ar < 0.5) LOGGER << "Warning: acceptance rate for SNP effects is too low "  << ar << endl;
//        else if (ar > 0.9) LOGGER << "Warning: acceptance rate for SNP effects is too high " << ar << endl;
//    }

}

VectorXd ApproxBayesC::SnpEffects::gradientU(const VectorXd &effects, VectorXd &rcorr, const VectorXd &ZPy, const vector<VectorXd> &ZPZ, 
                                             const VectorXi &windStart, const VectorXi &windSize, const unsigned chrStart, const unsigned chrSize,
                                             const double sigmaSq, const double vare){
    rcorr = ZPy;
    for (unsigned i=0; i<chrSize; ++i) {
        if (effects[i]) {
            rcorr.segment(windStart[i], windSize[i]) -= ZPZ[chrStart+i]*effects[i];
        }
    }
    return -rcorr/vare + effects/sigmaSq;
}

double ApproxBayesC::SnpEffects::computeU(const VectorXd &effects, const VectorXd &rcorr, const VectorXd &ZPy,                                             const double sigmaSq, const double vare){
    return -0.5f/vare*effects.dot(ZPy+rcorr) + 0.5/sigmaSq*effects.squaredNorm();
}

void ApproxBayesC::SnpEffects::sampleFromFC(vector<VectorXd> &wcorrBlocks, const vector<MatrixXd> &Qblocks, vector<VectorXd> &whatBlocks,
                                            const vector<LDBlockInfo*> keptLdBlockInfoVec, const VectorXd &nGWASblocks, const VectorXd &vareBlocks,
                                            const double sigmaSq, const double pi, const double varg, const VectorXd &snp2pq){
    
    long nBlocks = keptLdBlockInfoVec.size();

    whatBlocks.resize(nBlocks);
    ssqBlocks.resize(nBlocks);
    for (unsigned i=0; i<nBlocks; ++i) {
        whatBlocks[i].resize(wcorrBlocks[i].size());
    }

    double ssq[nBlocks], s2pq[nBlocks], nnz[nBlocks];
    memset(ssq,0, sizeof(double)*nBlocks);
    memset(s2pq,0,sizeof(double)*nBlocks);
    memset(nnz,0, sizeof(double)*nBlocks);

    double *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking precents the writing from other threads
  
    vector <double>  urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }
    
    double logPi = log(pi);
    double logPiComp = log(1.0-pi);
    double logSigmaSq = log(sigmaSq);
    double invSigmaSq = 1.0/sigmaSq;

    
#pragma omp parallel for schedule(dynamic)
    for(unsigned blk = 0; blk < nBlocks; blk++){
        Ref<const MatrixXd> Q = Qblocks[blk];
        Ref<VectorXd> wcorr = wcorrBlocks[blk];
        Ref<VectorXd> what = whatBlocks[blk];

        what.setZero();
        
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[blk];
        
        unsigned blockStart = blockInfo->startSnpIdx;
        unsigned blockEnd   = blockInfo->endSnpIdx;
        
        double invVareDn = nGWASblocks[blk] / vareBlocks[blk];

        double invLhs = 1.0/(invVareDn + invSigmaSq);
        double logInvLhsMsigma = logf(invLhs) - logSigmaSq;

        for(unsigned i = blockStart; i <= blockEnd; i++){
            double oldSample = valuesPtr[i];
            Ref<const VectorXd> Qi = Q.col(i - blockStart);
            double rhs = (Qi.dot(wcorr) + oldSample)*invVareDn;
            double uhat = invLhs * rhs;
            double logDelta1 = 0.5*(logInvLhsMsigma + uhat*rhs) + logPi;
            double logDelta0 = logPiComp;
            double probDelta1 = 1.0/(1.0 + expf(logDelta0-logDelta1));
            
//            LOGGER << i << " " << rhs << " " << invLhs << " " << uhat << " " << probDelta1 << " " << sigmaSq << endl;
//            int tmp;
//            cin >> tmp;
            
//            if (bernoulli.sample(probDelta1)) {
            if (urnd[i] < probDelta1) {
//                valuesPtr[i] = normal.sample(uhat, invLhs);
                valuesPtr[i] = uhat + nrnd[i]*sqrtf(invLhs);
                wcorr += Qi*(oldSample - valuesPtr[i]);
                what  += Qi* valuesPtr[i];
                ssq[blk] += (valuesPtr[i] * valuesPtr[i]);
                s2pq[blk] += snp2pq[i];
                ++nnz[blk];
            } else {
                if (oldSample) wcorr += Qi * oldSample;
                valuesPtr[i] = 0.0;
            }
        }
    }
    // LOGGER << "Varei 1 max" << varei.maxCoeff() << endl;
    //LOGGER << ssq << " " << nnz << endl;
    
    sumSq = 0.0;
    sum2pq = 0.0;
    numNonZeros = 0.0;
    nnzPerBlk.setZero(nBlocks);
    for (unsigned blk=0; blk<nBlocks; ++blk) {
        sumSq += ssq[blk];
        sum2pq += s2pq[blk];
        numNonZeros += nnz[blk];
        nnzPerBlk[blk] = nnz[blk];
        ssqBlocks[blk] = ssq[blk];
    }
    
    values = VectorXd::Map(valuesPtr, size);
}

void ApproxBayesC::ResidualVar::sampleFromFC(const double ypy, const VectorXd &effects, const VectorXd &ZPy, const VectorXd &rcorr, const double covg) {
    double sse = ypy - effects.dot(ZPy) - effects.dot(rcorr) + nobs*covg;
    if (sse < 0) {
        string vare_str = to_string(static_cast<long double>(sse/nobs));
        LOGGER.e(0,"\nError: Residual variance is negative (" + vare_str + "). This may indicate that effect sizes are \"blowing up\" likely due to a convergence problem. If SigmaSq variable is increasing with MCMC iterations, then this further indicates MCMC may not converge.");
    }
//    if (sse < 0) sse = 0; //-sse;
//    if (sse > ypy) sse = ypy;
    double dfTilde = df + nobs;
    double scaleTilde = sse + df*scale;
    value = InvChiSq::sample(dfTilde, scaleTilde);
}


void ApproxBayesC::GenotypicVar::compute(const VectorXd &effects, const VectorXd &ZPy, const VectorXd &rcorr, const double covg){
    double modelSS = effects.dot(ZPy) - effects.dot(rcorr) + nobs*covg;
    if (modelSS < 0) modelSS = 0; // -modelSS;
    value = modelSS/nobs;
}

void ApproxBayesC::BlockGenotypicVar::compute(const vector<VectorXd> &whatBlocks){
    for (unsigned i=0; i<numBlocks; ++i) {
        values[i] = whatBlocks[i].squaredNorm();
        //LOGGER << "varg " << i << " " << values[i] << endl;
    }
    total = values.sum();
}

void ApproxBayesC::BlockResidualVar::sampleFromFC(vector<VectorXd> &wcorrBlocks, VectorXd &ssqBlocks, const VectorXd &nGWASblocks, const VectorXd &numEigenvalBlock){
    for (unsigned i=0; i<numBlocks; ++i) {
        double sse = wcorrBlocks[i].squaredNorm() * nGWASblocks[i];
        double dfTilde = df + numEigenvalBlock[i];
        double scaleTilde = sse + df*scale;
        double sample = InvChiSq::sample(dfTilde, scaleTilde);
        if (ssqBlocks[i]/sample > threshold) {
            values[i] = sample;
        } else {
            values[i] = vary;
        }
        //LOGGER << "vare " << i << " " << values[i] << endl;
    }
    mean = values.mean();
}




void ApproxBayesC::Rounding::computeRcorr(const VectorXd &ZPy, const vector<SparseVector <double>  > &ZPZ,
                                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                                          const VectorXd &snpEffects, VectorXd &rcorr){
    if (count++ % 100) return;
    VectorXd rcorrOld = rcorr;
    rcorr = ZPy;
#pragma omp parallel for
    for (unsigned chr=0; chr<chromInfoVec.size(); ++chr) {
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        for (unsigned i=chrStart; i<=chrEnd; ++i) {
            for (SparseVector <double> ::InnerIterator it(ZPZ[i]); it; ++it) {
                //rcorr[windStart[i]+it.index()] -= it.value() * snpEffects[i];
                rcorr[it.index()] -= it.value() * snpEffects[i];
            }
//            rcorr.segment(windStart[i], windSize[i]) -= ZPZ[i]*snpEffects[i];
        }
    }
    value = sqrt(Gadget::calcVariance(rcorrOld-rcorr));
}

void ApproxBayesC::Rounding::computeRcorr(const VectorXd &ZPy, const vector<VectorXd> &ZPZ,
                                          const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                                          const VectorXd &snpEffects, VectorXd &rcorr){
    if (count++ % 100) return;
    VectorXd rcorrOld = rcorr;
    rcorr = ZPy;
#pragma omp parallel for
    for (unsigned chr=0; chr<chromInfoVec.size(); ++chr) {
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        for (unsigned i=chrStart; i<=chrEnd; ++i) {
            rcorr.segment(windStart[i], windSize[i]) -= ZPZ[i]*snpEffects[i];
        }
    }
    value = sqrt(Gadget::calcVariance(rcorrOld-rcorr));
}

void ApproxBayesC::Rounding::computeGhat(const MatrixXd &Z, const VectorXd &snpEffects, VectorXd &ghat){
    if (count++ % 100) return;
    VectorXd ghatOld = ghat;
    ghat.setZero(ghat.size());
    for (unsigned i=0; i<snpEffects.size(); ++i) {
        if (snpEffects[i]) ghat += Z.col(i)*snpEffects[i];
    }
    value = sqrt(Gadget::calcVariance(ghatOld-ghat));
}

void ApproxBayesC::PopulationStratification::compute(const VectorXd &rcorr, const VectorXd &ZPZdiag, const VectorXd &LDsamplVar, const double varg, const double vare, const VectorXd &chisq){
    
    value = (rcorr.array().square()/(ZPZdiag.array() * (LDsamplVar.array()*varg + value + vare))).mean() - 1.0;
    value = value < -0.01 ? -0.01 : value;
    
//    VectorXd varEta = ZPZdiag.array() * (LDsamplVar.array()*varg + value + vare);
////    VectorXd wt = varEta.array().square().inverse();
//    VectorXd wt = (rcorr.array().square()/ZPZdiag.array().square()).square().inverse();
////    VectorXd zsq = rcorr.array().square()/varEta.array();
//    VectorXd zsq = rcorr.array().square()/ZPZdiag.array() - LDsamplVar.array()*varg - vare;
//    value = zsq.cwiseProduct(wt).sum()/wt.sum();

    
//    VectorXd tmp = rcorr.array().square()/ZPZdiag.array() - LDsamplVar.array()*varg - vare;
//    double ssq = 0.0;
//    long size = rcorr.size();
//    long cnt = 0;
//    for (unsigned i=0; i<size; ++i) {
////        if (chisq[i] < 20 && !(i%20)) {
//            ssq += tmp[i];
//            ++cnt;
////        }
//    }
////    ssq /= double(cnt);
//    if (ssq < 0) ssq = 0.0;
//    double dfTilde = df + cnt;
//    double scaleTilde = ssq + df*scale;
//    value = InvChiSq::sample(dfTilde, scaleTilde);

    
//        ofstream out("rcorr.txt");
//        out << rcorr.array().square()/(ZPZdiag.array() * (LDsamplVar.array()*varg + value + vare)) << endl;
//        out.close();
    
}

void ApproxBayesC::PopulationStratification::compute(const VectorXd &rcorr, const VectorXd &ZPZdiag, const VectorXd &LDsamplVar, const double varg, const double vare, const vector<ChromInfo*> chromInfoVec){
    
    for (unsigned i=0; i<22; ++i) {
        unsigned start = chromInfoVec[i]->startSnpIdx;
        unsigned end = chromInfoVec[i]->endSnpIdx;
        unsigned size = end - start + 1;
        chrSpecific[i] = (rcorr.segment(start,size).array().square()/(ZPZdiag.segment(start,size).array() * (LDsamplVar.segment(start,size).array()*varg + chrSpecific[i] + vare))).mean() - 1.0;
    }    
}

void ApproxBayesC::NumResidualOutlier::compute(const VectorXd &rcorr, const VectorXd &ZPZdiag, const VectorXd &LDsamplVar, const double varg, const double vare, const vector<string> &snpName, VectorXi &leaveout, const vector<SparseVector <double>  > &ZPZ, const VectorXd &ZPy, const VectorXd &snpEffects) {
    iter++;
    //if (iter<10) return;
    
    VectorXd tss = ZPy.array().square()/ZPZdiag.array();
    VectorXd sse = rcorr.array().square()/ZPZdiag.array();
    value = 0;
    long size = tss.size();
    for (unsigned i=0; i<size; ++i) {
        if (sse[i] > 10 && tss[i] < 10) ++value;
        //if (sse[i] > 30 && sse[i] > tss[i]) ++value;
    }
    
//    MatrixXd X(tss.size(), 2);
//    X.col(0) = VectorXd::Ones(tss.size());
//    X.col(1) = snpEffects;
//    VectorXd bhat = ZPy.array()/ZPZdiag.array();
//    VectorXd b = X.householderQr().solve(bhat);
//    value = b[1];
    
//    LOGGER << sse.array().maxCoeff() << endl;
    
//    LOGGER << tss.mean() << " " << sse.mean() << endl;
    
//    value = 0;
//    
//    VectorXd tmp = rcorr.array().square()/(ZPZdiag.array() * (LDsamplVar.array()*varg + value + vare));
//    //VectorXd tmp = rcorr.array().square()/ZPZdiag.array() - LDsamplVar.array()*varg;
//    
//    long size = tmp.size();
//    stringstream ss;
//    for (unsigned i=0; i<size; ++i) {
//        if (tmp[i]>20) {
//            ++value;
//            //leaveout[i] = 1;
//            ss << " " << snpName[i];
//            
////            for (SparseVector <double> ::InnerIterator it(ZPZ[i]); it; ++it) {
////                leaveout[it.index()] = 1;
////            }
//            
//        }
//    }
//    if (value) out << iter << ss.str() << endl;
}

void ApproxBayesC::ldScoreReg(const VectorXd &chisq, const VectorXd &LDscore, const VectorXd &LDsamplVar,
                              const double varg, const double vare, double &ps){
    
    long nrow = chisq.size();
    
//    ps = (chisq - LDscore - LDsamplVar*varg - VectorXd::Ones(nrow)*vare).mean();
//    ps = (chisq - LDscore*0.17/double(nrow) - VectorXd::Ones(nrow)).mean();

//    return;
    
    VectorXd y = chisq - LDsamplVar*varg - VectorXd::Ones(nrow)*vare;
//    VectorXd y = chisq;
//    VectorXd weight = 2.0*(LDscore*vargj + LDsamplVar*varg + VectorXd::Ones(nrow)*vare).array().square();
//    VectorXd weight = 2.0*(LDscore*vargj + VectorXd::Ones(nrow)).array().square();
//    VectorXd weightInv = weight.cwiseInverse();
    
    MatrixXd X(nrow, 2);
//    X.col(0) = weight;
//    X.col(1) = LDscore.cwiseProduct(weight);
//    y.array() *= weight.array();

    X.col(0) = VectorXd::Ones(nrow);
    X.col(1) = LDscore;
    
//    unsigned m = 0;
//    for (unsigned i=0; i<nrow; ++i) {
//        if (chisq[i] < 30) ++m;
//    }
//    VectorXd ysub(m);
//    MatrixXd Xsub(m, 2);
//    unsigned j=0;
//    for (unsigned i=0; i<nrow; ++i) {
//        if (chisq[i] < 30) {
//            ysub[j] = y[i];
//            Xsub.row(j) = X.row(i);
//            ++j;
//        }
//    }
//    VectorXd b = Xsub.householderQr().solve(ysub);

    
    VectorXd b = X.householderQr().solve(y);
//    VectorXd b = (X.transpose()*weightInv.asDiagonal()*X).inverse()*X.transpose()*weightInv.asDiagonal()*y;

    ps = b[0];
    
//    LOGGER << b.transpose() << endl;
    
}

void ApproxBayesC::InterChrGenetCov::compute(const double ypy, const VectorXd &effects, const VectorXd &ZPy, const VectorXd &rcorr) {
    if (!spouseCorrelation) return;
    double bZPy = effects.dot(ZPy);
    double brcorr = effects.dot(rcorr);
    double varg = (bZPy - brcorr)/nobs;
//    double vare = (ypy - bZPy - brcorr)/nobs;
    double varp = ypy/nobs;
    double hsq = varg/varp;
    double R = spouseCorrelation*hsq / (1 - spouseCorrelation*hsq);
    value = varg * R; // * 0.95;
}

void ApproxBayesC::NnzGwas::compute(const VectorXd &effects, const vector<SparseVector <double>  > &ZPZ, const VectorXd &ZPZdiag) {
    if (iter++ % 100) return;
    value = 0;
    long numSnps = effects.size();
    unsigned i, j;
    for (i=0; i<numSnps; ++i) {
        for (SparseVector <double> ::InnerIterator it(ZPZ[i]); it; ++it) {
            j = it.index();
            if (effects[j]){
                if (it.value()*it.value() > 0.1*ZPZdiag[i]*ZPZdiag[j]) {
                    ++value;
                    break;
                }
            }
        }
    }
}

void ApproxBayesC::PiGwas::compute(const double nnzGwas, const unsigned int numSnps) {
    if (iter++ % 100) return;
    value = nnzGwas/double(numSnps);
}

void ApproxBayesC::checkHsq(vector <double>  &hsqMCMC) {
    long niter = hsqMCMC.size();
    VectorXd y(niter);
    MatrixXd X(niter, 2);
    for (unsigned i=0; i<niter; ++i) {
        y[i] = hsqMCMC[i];
        X(i,0) = 1;
        X(i,1) = i;
    }
    VectorXd b = X.householderQr().solve(y);
    double slope = b[1];
    double vare = (y.squaredNorm() - (X*b).squaredNorm())/double(niter);
    double se = sqrt(vare/X.col(1).squaredNorm());
    if ((slope/se) > 3) {   // 3 corresponds to P = 0.001 at one-way test
        string slope_str = to_string(static_cast <double> (slope));
        string se_str = to_string(static_cast <double> (se));
        LOGGER.e(0,"\nError: The SNP-heritability is increasing over MCMC iterations (slope: " + slope_str + "; se: " + se_str + "). This may indicate that effect sizes are \"blowing up\" likely due to a convergence problem.");
    }

}

void ApproxBayesC::sampleUnknowns(){
    static int iter = 0;
//    fixedEffects.sampleFromFC(data.XPX, data.XPXdiag, data.ZPX, data.XPy, snpEffects.values, vare.value, rcorr);
    unsigned cnt=0;
    //do {
        //snpEffects.sampleFromFC(rcorr, data.ZPZ, data.ZPZdiag, data.ZPy, data.windStart, data.windSize, data.chromInfoVec, data.se, tss, data.n, data.snp2pq, sigmaSq.value, pi.value, vare.value);
        //snpEffects.hmcSampler(rcorr, data.ZPy, data.ZPZ, data.windStart, data.windSize, data.chromInfoVec, sigmaSq.value, pi.value, vare.value);
    if (lowRankModel)
        snpEffects.sampleFromFC(wcorrBlocks, data.Qblocks, whatBlocks, data.keptLdBlockInfoVec, data.nGWASblock, vareBlk.values, sigmaSq.value, pi.value, varg.value, data.snp2pq);
    else if (sparse)
        snpEffects.sampleFromFC(rcorr, data.ZPZsp, data.ZPZdiag, data.ZPy, data.windStart, data.windSize, data.chromInfoVec, data.se, data.tss, varei,
                                data.n, data.snp2pq, data.LDsamplVar, sigmaSq.value, pi.value, vare.value, varg.value, ps.value, overdispersion);
    else
        snpEffects.sampleFromFC(rcorr, data.ZPZ, data.ZPZdiag, data.ZPy, data.windStart, data.windSize, data.chromInfoVec, data.se, data.tss, varei,
                                data.n, data.snp2pq, data.LDsamplVar, sigmaSq.value, pi.value, vare.value, varg.value, ps.value, overdispersion);
    //    if (++cnt == 100) LOGGER.e(0," Zero SNP effect in the model for 100 cycles of sampling");
    //} while (snpEffects.numNonZeros == 0);
    if (diagnose) nro.compute(rcorr, data.ZPZdiag, data.LDsamplVar, varg.value, vare.value, snpEffects.header, snpEffects.leaveout, data.ZPZsp, data.ZPy, snpEffects.values);
    
    if (robustMode) {
        if (noscale) {
            sigmaSq.value = varg.value/(data.snp2pq.array().sum()*pi.value);
        } else {
            sigmaSq.value = varg.value/(data.numIncdSnps*pi.value);  // LDpred2's parameterisation
        }
    } else {
        sigmaSq.sampleFromFC(snpEffects.sumSq, snpEffects.numNonZeros);
    }

    if (estimatePi) pi.sampleFromFC(data.numIncdSnps, snpEffects.numNonZeros);
    nnzSnp.getValue(snpEffects.numNonZeros);
    if (noscale) {
        sigmaSqG.compute(sigmaSq.value, snpEffects.sum2pq);
    } else {
        sigmaSqG.value = sigmaSq.value * snpEffects.numNonZeros;
    }

    if (lowRankModel) {
        vargBlk.compute(whatBlocks);
        vareBlk.sampleFromFC(wcorrBlocks, snpEffects.ssqBlocks, data.nGWASblock, data.numEigenvalBlock);
        varg.value = vargBlk.total;
        vare.value = vareBlk.mean;
    }
    else {
        covg.compute(data.ypy, snpEffects.values, data.ZPy, rcorr);
        varg.compute(snpEffects.values, data.ZPy, rcorr, covg.value);
        //    varg.value = sigmaSqG.value;
        vare.sampleFromFC(data.ypy, snpEffects.values, data.ZPy, rcorr, covg.value);
    }
    //hsq.compute(varg.value, vare.value);
    hsq.value = varg.value / data.varPhenotypic;

    if (iter >= 2000) sigmaSq.scale = scalePrior;
    //sigmaSq.scale = scalePrior;
    scale.getValue(sigmaSq.scale);
        
//    if (sparse)
//        rounding.computeRcorr(data.ZPy, data.ZPZsp, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
//    else
//        rounding.computeRcorr(data.ZPy, data.ZPZ, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
    if (modelPS) ps.compute(rcorr, data.ZPZdiag, data.LDsamplVar, varg.value, vare.value, data.chisq);
//    if (modelPS) ldScoreReg(data.chisq, data.LDscore, data.LDsamplVar, varg.value, vare.value, ps.value);
    
//    if (sparse) {
//        nnzgwas.compute(snpEffects.values, data.ZPZsp, data.ZPZdiag);
//        pigwas.compute(nnzgwas.value, data.numIncdSnps);
//    }
    
    double scaleIteri = 0;
    if (++iter < 2000) {
        if (noscale)
        {
            scaleIteri = 0.5f * varg.value / (data.snp2pq.array().sum()*(pi.value));
        } else
        {
            scaleIteri = 0.5f * varg.value / (data.snp2pq.size()*(pi.value));
        }
        genVarPrior += (varg.value - genVarPrior)/iter;
        scalePrior += (scaleIteri - scalePrior)/iter;
    }
    //if(iter>1990 && iter < 2010){
    //    LOGGER << "iter " << iter << " , scalePrior " << scalePrior << " , sigmaSq.scale " << sigmaSq.scale << " , genVarPrior " << genVarPrior << " , varg " << varg.value << " , pi " << pi.value << endl;
    //}
    
//    if (iter > 100 & !(iter % 10)) {
//        hsqMCMC.push_back(hsq.value);
//        if (!(iter % 1000)) checkHsq(hsqMCMC);
//    }
}


// *******************************************************
// Bayes R - Approximate
// *******************************************************

void ApproxBayesR::sampleUnknowns(){
    //sigmaSq.value = 0.000275;   // TMP_JZ
    static int iter = 0;
//    fixedEffects.sampleFromFC(data.XPX, data.XPXdiag, data.ZPX, data.XPy, snpEffects.values, vare.value, rcorr);
    unsigned cnt=0;
    //do {
        if (data.Z.size()) {
            snpEffects.sampleFromFC(data.ZPy, data.ZPZdiag, data.Z, data.Z.rows(), data.numKeptInds, sigmaSq.value, Pis.values, gamma.values, vare.value, snpStore, ghat, varg.value, originalModel);
        }
        else {
            if (lowRankModel) {
                snpEffects.sampleFromFC(wcorrBlocks, data.Qblocks, whatBlocks, data.keptLdBlockInfoVec, data.nGWASblock, vareBlk.values, sigmaSq.value, Pis.values, gamma.values, snpStore, varg.value, originalModel);
            }
            else if (sparse)
                snpEffects.sampleFromFC(rcorr, data.ZPZsp, data.ZPZdiag, data.ZPy, data.windStart, data.windSize, data.chromInfoVec, data.se, data.tss, varei, data.n, data.snp2pq, data.LDsamplVar, sigmaSq.value, Pis.values, gamma.values, vare.value, snpStore,
                                        varg.value, ps.value, overdispersion, originalModel);
            else
                snpEffects.sampleFromFC(rcorr, data.ZPZ, data.ZPZdiag, data.ZPy, data.windStart, data.windSize, data.chromInfoVec, data.se, data.tss, varei, data.n, data.snp2pq, data.LDsamplVar, sigmaSq.value, Pis.values, gamma.values, vare.value, snpStore,
                                        varg.value, ps.value, overdispersion, originalModel);
        }
    //    if (++cnt == 100) LOGGER.e(0,"Error: Zero SNP effect in the model for 100 cycles of sampling");
    //} while (snpEffects.numNonZeros == 0);
        
    if (algorithm == cg) {
        snpEffects.adjustByCG(data.ZPy, data.ZPZsp, rcorr);
    }
    
    if (diagnose) nro.compute(rcorr, data.ZPZdiag, data.LDsamplVar, varg.value, vare.value, snpEffects.header, snpEffects.leaveout, data.ZPZsp, data.ZPy, snpEffects.values);
    
    if (robustMode) {
        if (noscale) {
            sigmaSq.value = varg.value/(data.snp2pq.array().sum()*gamma.values.dot(Pis.values));
        } else {
            sigmaSq.value = varg.value/(data.numIncdSnps*gamma.values.dot(Pis.values));  // LDpred2's parameterisation
        }
    } else if (originalModel) {
        sigmaSq.sampleFromFC(snpEffects.sumSq, snpEffects.numNonZeros);
    } else {
        if (estimateSigmaSq) sigmaSq.sampleFromFC(snpEffects.sumSq, snpEffects.numNonZeros);
    }
        
    if (estimatePi) Pis.sampleFromFC(snpStore);
    numSnps.getValues(snpStore);
    nnzSnp.getValue(snpEffects.numNonZeros);
    sigmaSqG.compute(sigmaSq.value, snpEffects.sum2pq);

    if (estimateHsq) {
        if (data.Z.size()) {   // TMP_JZ
            //        ghat.setZero(data.Z.rows());
            //        for (unsigned i=0; i<snpEffects.size; ++i) {
            //            if (snpEffects.values[i]) ghat += data.Z.col(i)*snpEffects.values[i];
            //        }
            varg.value = Gadget::calcVariance(ghat);
            double n_ref = data.Z.rows();
            double n_ratio = n_ref/double(data.numKeptInds);
            vare.value = (data.ypy*n_ratio - 2.0f*snpEffects.values.dot(data.ZPy)*n_ratio + ghat.dot(ghat))/n_ref;
            //cout << "varg " << varg.value << " vare " << vare.value << endl;
            //cout << "n_ratio " << n_ratio << " ypy " << data.ypy*n_ratio << " ypg " << 2.0f*snpEffects.values.dot(data.ZPy)*n_ratio << " gpg " << ghat.dot(ghat) << " n_ref " << n_ref << endl;
        }
        else if (lowRankModel) {
            vargBlk.compute(whatBlocks);
            vareBlk.sampleFromFC(wcorrBlocks, snpEffects.ssqBlocks, data.nGWASblock, data.numEigenvalBlock);
            varg.value = vargBlk.total;
            vare.value = vareBlk.mean;
        }
        else {
            covg.compute(data.ypy, snpEffects.values, data.ZPy, rcorr);
            varg.compute(snpEffects.values, data.ZPy, rcorr, covg.value);
            vare.sampleFromFC(data.ypy, snpEffects.values, data.ZPy, rcorr, covg.value);
        }
        
        //hsq.compute(varg.value, vare.value);
        hsq.value = varg.value / data.varPhenotypic;
    }
    
    if (iter >= 2000) sigmaSq.scale = scalePrior;
    scale.getValue(sigmaSq.scale);
    // cout << "iter " << iter << " scalePrior " << scalePrior << "sigmaSq.scale " << sigmaSq.scale << endl;

//    if (sparse)
//        rounding.computeRcorr(data.ZPy, data.ZPZsp, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
//    else
//        rounding.computeRcorr(data.ZPy, data.ZPZ, data.windStart, data.windSize, data.chromInfoVec, snpEffects.values, rcorr);
    if (modelPS) ps.compute(rcorr, data.ZPZdiag, data.LDsamplVar, varg.value, vare.value, data.chisq);

    nnzSnp.getValue(snpEffects.numNonZeros);
    if (noscale) {
        sigmaSqG.compute(sigmaSq.value, snpEffects.sum2pq);
    } else {
        sigmaSqG.value = sigmaSq.value * snpEffects.numNonZeros;
    }

//    numSnpVg.compute(snpEffects.values, data.ZPZdiag, varg.value, vare.nobs);
    if (originalModel) {
        Vgs.compute(snpEffects.values, data.ZPy, rcorr, snpEffects.snpset, varg.value, vare.nobs);
//        if (sparse)
//            Vgs.compute(snpEffects.values, data.ZPZsp, snpEffects.snpset, varg.value, vare.nobs);
//        else
//            Vgs.compute(snpEffects.values, data.ZPZ, snpEffects.snpset, varg.value, vare.nobs);
    }

    double scaleIteri = 0;
    if (++iter < 2000) {
        if (noscale)
        {
            scaleIteri = 0.5f * varg.value / (data.snp2pq.array().sum()*gamma.values.dot(Pis.values));
        } else
        {
            scaleIteri = 0.5f * varg.value / (data.snp2pq.size()*gamma.values.dot(Pis.values));
        }
        genVarPrior += (varg.value - genVarPrior)/iter;
        scalePrior += (scaleIteri - scalePrior)/iter;
    }
    
//    if (iter > 100 & !(iter % 10)) {
//        hsqMCMC.push_back(hsq.value);
//        if (!(iter % 1000)) checkHsq(hsqMCMC);
//    }
}

void ApproxBayesR::VgMixComps::compute(const VectorXd &snpEffects, const VectorXd &ZPy, const VectorXd &rcorr, const vector<vector<unsigned> > snpset, const double varg, const double nobs) {
    values.setZero(ndist);
    for (unsigned k=1; k<ndist; ++k) {
        unsigned size = snpset[k].size();
        for (unsigned j=0; j<size; ++j) {
            unsigned snpIdx = snpset[k][j];
            double varj = snpEffects[snpIdx] * (ZPy[snpIdx] - rcorr[snpIdx]);
            values[k] += varj;
        }
        values[k] /= varg * nobs;
        (*this)[k]->value = values[k];
    }
}


//void ApproxBayesR::VgMixComps::compute(const VectorXd &snpEffects, const vector<SparseVector<double> > &ZPZsp, const vector<vector<unsigned> > snpset, const double varg, const double nobs) {
//    values.setZero(ndist);
//    for (unsigned k=0; k<ndist; ++k) {
//        if (k!=zeroIdx && k!=minIdx) {
//            long numSnps = snpset[k].size();
//            double vargk = 0.0;
//            for (unsigned i=0, j=0; i<numSnps; ++i) {
//                unsigned ki = snpset[k][i];
//                unsigned kj = snpset[k][j];
//                for (SparseVector<double>::InnerIterator it(ZPZsp[ki]); it; ++it) {
//                    if (it.index() == kj) {
//                        vargk += snpEffects[ki]*snpEffects[kj]*it.value();
//                        kj = snpset[k][++j];
//                        if (j==numSnps) break;
//                    }
//                }
//            }
//            (*this)[k]->value = values[k] = vargk/(varg*nobs);
//        }
//    }
//    double sum = values.sum();
//    (*this)[minIdx]->value = values[minIdx] = 1.0 - sum;
//}
//
//void ApproxBayesR::VgMixComps::compute(const VectorXd &snpEffects, const vector<VectorXd> &ZPZ, const vector<vector<unsigned> > snpset, const double varg, const double nobs) {
//    values.setZero(ndist);
//    for (unsigned k=0; k<ndist; ++k) {
//        if (k!=zeroIdx && k!=minIdx) {
//            long numSnps = snpset[k].size();
//            double vargk = 0.0;
//            for (unsigned i=0; i<numSnps; ++i) {
//                unsigned ki = snpset[k][i];
//                for (unsigned j=0; j<numSnps; ++j) {
//                    unsigned kj = snpset[k][j];
//                    vargk += snpEffects[ki]*snpEffects[kj]*ZPZ[ki][kj];
//                }
//            }
//            (*this)[k]->value = values[k] = vargk/(varg*nobs);
//        }
//    }
//    double sum = values.sum();
//    (*this)[minIdx]->value = values[minIdx] = 1.0 - sum;
//}

// ==============================================================
// Sparse vector version
// ==============================================================

void ApproxBayesR::SnpEffects::sampleFromFC(VectorXd &rcorr, const vector<SparseVector<double>> &ZPZ, const VectorXd &ZPZdiag, const VectorXd &ZPy,
                                            const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                                            const VectorXd &se, const VectorXd &tss, VectorXd &varei, const VectorXd &n, const VectorXd &snp2pq, const VectorXd &LDsamplVar,
                                            const double sigmaSq, const VectorXd &pis, const VectorXd &gamma, const double vare, VectorXd &snpStore,
                                            const double varg, const double ps, const double overdispersion,
                                            const bool originalModel){
    // -----------------------------------------
    // Initialise the parameters in MCMC sampler
    // -----------------------------------------
    static unsigned iter = 0;
    long numChr = chromInfoVec.size();

    double ssq[numChr], s2pq[numChr], nnz[numChr];
    memset(ssq,0,sizeof(double)*numChr);
    memset(s2pq,0,sizeof(double)*numChr);
    memset(nnz,0, sizeof(double)*numChr);

    double *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking prevents the writing from other threads

    vector<double> urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }
    
    // R specific parameters
    int ndist;
    VectorXd gp;
    snpStore.setZero(pis.size());
    // --------------------------------------------------------------------------------
    // Scale the variances in each of the normal distributions by the genetic variance
    // and initialise the class membership probabilities
    // --------------------------------------------------------------------------------
    ndist = pis.size();
    if (originalModel)
        gp = gamma * 0.01 * varg;
    else
        gp = gamma * sigmaSq;
    
    vector<vector<vector<unsigned> > > snpsetChr(numChr);
    for (unsigned i=0; i<numChr; ++i) {
        snpsetChr[i].resize(ndist);
        for (unsigned k=0; k<ndist; ++k) {
            snpsetChr[i][k].resize(0);
        }
    }

    VectorXd invGamma = gamma.array().inverse();
    invGamma[0] = 0.0;
    
    lambdaVec.setZero(size);
    uhatVec.setZero(size);
    invGammaVec.setZero(size);
    deltaNZ.setZero(size);
    
    deltaNzIdx.clear();
    deltaNzIdx.reserve(size);
    
    // --------------------------------------------------------------------------------
    // Cycle over all variants in the window and sample the genetics effects
    // --------------------------------------------------------------------------------

#pragma omp parallel for schedule(dynamic)
    for (unsigned chr=0; chr<numChr; ++chr) 
    {
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned windEnd, j;
        
        // R specific parameters
        int indistflag;
        double rhs, v1,  b_ls, ssculm, r;
        VectorXd ll, pll, snpindist, var_b_ls;
        ll.setZero(pis.size());
        pll.setZero(pis.size());

        double oldSample, varei;
        
        for (unsigned i=chrStart; i<=chrEnd; ++i) {
            oldSample = valuesPtr[i]; 
            varei = LDsamplVar[i]*varg + vare + ps + overdispersion;
                        
//            varei = (tss[i] + oldSample*oldSample*ZPZdiag[i])/n[i];
//            double ssei = rcorr[i]*rcorr[i]/ZPZdiag[i];

//            double dfTilde = 10 + 1;
//            double scaleTilde = ssei + 10*(LDsamplVar[i]*varg + vare + ps + overdispersion);
//            Stat::InvChiSq invchisq;
//            varei = invchisq.sample(dfTilde, scaleTilde);
//
//            cout << i << " " << varei << endl;
            
            // ------------------------------
            // Derived Bayes R implementation
            // ------------------------------
            // ----------------------------------------------------
            // Add back the content for the corrected rhs for SNP k
            // ----------------------------------------------------
            rhs = rcorr[i] + ZPZdiag[i] * oldSample;
            // ------------------------------------------------------
            // Calculate the beta least squares updates and variances
            // ------------------------------------------------------
            b_ls = rhs / ZPZdiag[i];
            var_b_ls = gp.array() + varei / ZPZdiag[i];
            // ------------------------------------------------------
            // Calculate the likelihoods for each distribution
            // ------------------------------------------------------
            ll = (-1.0 / 2.0) * var_b_ls.array().log()  - (b_ls * b_ls)  / (2 * var_b_ls.array()) + pis.array().log();
            // --------------------------------------------------------------
            // Calculate probability that snp is in each of the distributions
            // in this iteration
            // --------------------------------------------------------------
            // pll = (ll.array().exp().cwiseProduct(pis.array())) / ((ll.array().exp()).cwiseProduct(pis.array())).sum();
            for (unsigned k=0; k<pis.size(); ++k) {
              pll[k] = 1.0 / (exp(ll.array() - ll[k])).sum();
            }
            // --------------------------------------------------------------
            // Sample the group based on the calculated probabilities
            // --------------------------------------------------------------
            ssculm = 0.0;
            r = urnd[i];
            indistflag = 1;
            for (int kk = 0; kk < ndist; kk++)
            {
                ssculm += pll(kk);
                if (r < ssculm)
                {
                    indistflag = kk + 1;
                    break;
                }
            }
            snpsetChr[chr][indistflag-1].push_back(i);
            // --------------------------------------------------------------
            // Sample the effect given the group and adjust the rhs                                                                                                                 
            // --------------------------------------------------------------                                                                                                       
            if (indistflag != 1)                                                                                                                                                    
            {                                                                                                                                                                       
                v1 = ZPZdiag[i] + varei / gp((indistflag - 1));                                                                                                                     
//                valuesPtr[i] = normal.sample(rhs / v1, varei / v1);                                                                                                                 
                valuesPtr[i] = rhs / v1 + nrnd[i]*sqrtf(varei / v1);
                double sampleDiff = oldSample - valuesPtr[i];
                for (SparseVector<double>::InnerIterator it(ZPZ[i]); it; ++it) {                                                                                                     
                    rcorr[it.index()] += it.value() * sampleDiff;                                                                                                                   
                }                                                                                                                                                                   
                ssq[chr]  += (valuesPtr[i]*valuesPtr[i]) / gamma[indistflag - 1];
                s2pq[chr] += snp2pq[i];
                deltaNZ[i] = 1;
                ++nnz[chr];
                deltaNzIdx.push_back(i);
            } else {                                                                                                                                                                
                if (oldSample) {                                                                                                                                                    
                    for (SparseVector<double>::InnerIterator it(ZPZ[i]); it; ++it) {                                                                                                 
                        rcorr[it.index()] += it.value() * oldSample;                                                                                                                
                    }
                }                                                                                                                                                                   
                valuesPtr[i] = 0.0;                                                                                                                                                 
            }
            
            uhatVec[i] = rhs/v1;
            lambdaVec[i] = vare/gp[indistflag-1];
            invGammaVec[i] = invGamma[indistflag-1];
        }
    }
    // ---------------------------------------------------------------------
    // Tally up the effect sum of squares and the number of non-zero effects
    // ---------------------------------------------------------------------
    sumSq = 0.0;                                                                                                                                                                    
    sum2pq = 0.0;                                                                                                                                                                   
    numNonZeros = 0;                                                                                                                                                                
    nnzPerChr.setZero(numChr);                                                                                                                                                      
    snpStore.setZero(ndist);
    snpset.resize(ndist);
    for (unsigned k=0; k<ndist; ++k) {
        snpset[k].resize(0);
    }
    for (unsigned i=0; i<numChr; ++i) {
        sumSq += ssq[i];
        sum2pq += s2pq[i];                                                                                                                                                          
        numNonZeros += nnz[i];                                                                                                                                                      
        nnzPerChr[i] = nnz[i];                                                                                                                                                      
        for (unsigned k=0; k<ndist; ++k) {
            for (unsigned j=0; j<snpsetChr[i][k].size(); ++j) {
                snpset[k].push_back(snpsetChr[i][k][j]);
                snpStore[k]++;
            }
        }
    }
    ++iter;

    values = VectorXd::Map(valuesPtr, size);
}

// ==============================================================
// Vector of vectors version
// ==============================================================

void ApproxBayesR::SnpEffects::sampleFromFC(VectorXd &rcorr, const vector<VectorXd> &ZPZ, const VectorXd &ZPZdiag, const VectorXd &ZPy,
                                            const VectorXi &windStart, const VectorXi &windSize, const vector<ChromInfo*> &chromInfoVec,
                                            const VectorXd &se, const VectorXd &tss, VectorXd &varei, const VectorXd &n, const VectorXd &snp2pq, const VectorXd &LDsamplVar,
                                            const double sigmaSq, const VectorXd &pis, const VectorXd &gamma, const double vare, VectorXd &snpStore,
                                            const double varg, const double ps, const double overdispersion,
                                            const bool originalModel){
    // -----------------------------------------
    // Initialise the parameters in MCMC sampler
    // -----------------------------------------
    static unsigned iter = 0;                                                                                                                                                       
    long numChr = chromInfoVec.size();                                                                                                                                              
    
    double ssq[numChr], nnz[numChr], s2pq[numChr];                                                                                                                                   
    memset(ssq,0,sizeof(double)*numChr);                                                                                                                                             
    memset(nnz,0,sizeof(double)*numChr);                                                                                                                                             
    memset(s2pq,0,sizeof(double)*numChr);                                                                                                                                            

    double *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking precents the writing from other threads                     
  
    vector<double> urnd(size), nrnd(size);                                                                                                                                           
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work                                                                                                               
        urnd[i] = Stat::ranf();                                                                                                                                                     
        nrnd[i] = Stat::snorm();                                                                                                                                                    
    }

    // R specific parameters
    int ndist;
    VectorXd gp;
    snpStore.setZero(pis.size());
    // --------------------------------------------------------------------------------
    // Scale the variances in each of the normal distributions by the genetic variance
    // and initialise the class membership probabilities
    // --------------------------------------------------------------------------------
    ndist = pis.size();
    if (originalModel)
        gp = gamma * 0.01 * varg;
    else
        gp = gamma * sigmaSq;

    vector<vector<vector<unsigned> > > snpsetChr(numChr);
    for (unsigned i=0; i<numChr; ++i) {
        snpsetChr[i].resize(ndist);
        for (unsigned k=0; k<ndist; ++k) {
            snpsetChr[i][k].resize(0);
        }
    }

    // --------------------------------------------------------------------------------
    // Cycle over all variants in the window and sample the genetics effects
    // --------------------------------------------------------------------------------

#pragma omp parallel for schedule(dynamic)
    for (unsigned chr=0; chr<numChr; ++chr) 
    {
        ChromInfo *chromInfo = chromInfoVec[chr];
        unsigned chrStart = chromInfo->startSnpIdx;
        unsigned chrEnd   = chromInfo->endSnpIdx;
        unsigned windEnd, j;
        double oldSample, varei;
        double rhs, invLhs, uhat;
        
        int indistflag;
        double v1,  b_ls, ssculm, r;
        VectorXd ll, pll, snpindist, var_b_ls;
        ll.setZero(pis.size());
        pll.setZero(pis.size());

        for (unsigned i=chrStart; i<=chrEnd; ++i) {
            oldSample = valuesPtr[i];
            // ---------------------------------------------
            // Calculate residual variance including a 
            // correction for the sampling variation and
            // LD ignored
            // ---------------------------------------------
            varei = LDsamplVar[i]*varg + vare + ps + overdispersion;
            // ------------------------------
            // Derived Bayes R implementation
            // ------------------------------
            // ----------------------------------------------------
            // Add back the content for the corrected rhs for SNP k
            // ----------------------------------------------------
            rhs = rcorr[i] + ZPZdiag[i] * oldSample;
            // ------------------------------------------------------
            // Calculate the beta least squares updates and variances
            // ------------------------------------------------------
            b_ls = rhs / ZPZdiag[i];
            var_b_ls = gp.array() + varei / ZPZdiag[i];
            // ------------------------------------------------------
            // Calculate the likelihoods for each distribution
            // ------------------------------------------------------
            // ll  = (-1.0 / 2.0) * var_b_ls.array().log()  - (b_ls * b_ls)  / (2 * var_b_ls.array());
            ll = (-1.0 / 2.0) * var_b_ls.array().log()  - (b_ls * b_ls)  / (2 * var_b_ls.array()) + pis.array().log();
            // --------------------------------------------------------------
            // Calculate probability that snp is in each of the distributions
            // in this iteration
            // --------------------------------------------------------------
            // pll = (ll.array().exp().cwiseProduct(pis.array())) / ((ll.array().exp()).cwiseProduct(pis.array())).sum();
            for (unsigned k=0; k<pis.size(); ++k) {
              pll[k] = 1.0 / (exp(ll.array() - ll[k])).sum();
            }
            // if (i < 10) {
            //   cout << "P likelihood 1 " << pll << endl;
            //   cout << "P likelihood 2 " << pll2 << endl;
            // }
            // --------------------------------------------------------------
            // Sample the group based on the calculated probabilities
            // --------------------------------------------------------------
            ssculm = 0.0;
            r = urnd[i];
            indistflag = 1;
            for (int kk = 0; kk < ndist; kk++)
            {
                ssculm += pll(kk);
                if (r < ssculm)
                {
                    indistflag = kk + 1;
                    break;
                }
            }
            snpsetChr[chr][indistflag-1].push_back(i);
            // --------------------------------------------------------------
            // Sample the effect given the group and adjust the rhs
            // --------------------------------------------------------------
            if (indistflag != 1)                                                                                                                                                    
            {                                                                                                                                                                       
                v1 = ZPZdiag[i] + varei / gp((indistflag - 1));                                                                                                                     
//                valuesPtr[i] = normal.sample(rhs / v1, varei / v1);                                                                                                                 
                valuesPtr[i] = rhs / v1 + nrnd[i]*sqrtf(varei / v1);
                rcorr.segment(windStart[i], windSize[i]) += ZPZ[i] * (oldSample - valuesPtr[i]);
                ssq[chr] += (valuesPtr[i] * valuesPtr[i]) / gamma[indistflag - 1];                                                                                                  
                s2pq[chr] += snp2pq[i];                                                                                                                                             
                ++nnz[chr];                                                                                                                                                         
            } else {                                                                                                                                                                
                if (oldSample) rcorr.segment(windStart[i], windSize[i]) += ZPZ[i] * oldSample;                                                                                      
                valuesPtr[i] = 0.0;                                                                                                                                                 
            }  
        }
    }
    // ---------------------------------------------------------------------                                                                                                        
    // Tally up the effect sum of squares and the number of non-zero effects                                                                                                        
    // ---------------------------------------------------------------------                                                                                                        
    sumSq = 0.0;                                                                                                                                                                    
    sum2pq = 0.0;                                                                                                                                                                   
    numNonZeros = 0.0;                                                                                                                                                              
    nnzPerChr.setZero(numChr);                                                                                                                                                      
    snpStore.setZero(ndist);
    snpset.resize(ndist);
    for (unsigned k=0; k<ndist; ++k) {
        snpset[k].resize(0);
    }
    for (unsigned i=0; i<numChr; ++i) {
        sumSq += ssq[i];                                                                                                                                                            
        sum2pq += s2pq[i];                                                                                                                                                          
        numNonZeros += nnz[i];                                                                                                                                                      
        nnzPerChr[i] = nnz[i];                                                                                                                                                      
        for (unsigned k=0; k<ndist; ++k) {
            for (unsigned j=0; j<snpsetChr[i][k].size(); ++j) {
                snpset[k].push_back(snpsetChr[i][k][j]);
                snpStore[k]++;
            }
        }
    }
    ++iter;
    
    values = VectorXd::Map(valuesPtr, size); 
}

void ApproxBayesR::SnpEffects::sampleFromFC(const VectorXd &ZPy, const VectorXd &ZPZdiag, const MatrixXd &Z, const double n_ref, const double n_gwas,
                                            const double sigmaSq, const VectorXd &pis, const VectorXd &gamma, const double vare,
                                            VectorXd &snpStore, VectorXd &ghat, const double varg, const bool originalModel) {
        sumSq = 0.0;
        numNonZeros = 0;
            
        ghat.setZero(n_ref);
        double oldSample;
        double my_rhs, rhs;
        // -----------------------------------------
        // Initialise the parameters in MCMC sampler
        // -----------------------------------------
        // ----------------
        // Bayes R specific
        // ----------------
        int ndist, indistflag;
        double v1,  b_ls, ssculm, r;
        VectorXd gp, ll, ll2, pll, snpindist, var_b_ls;
        ndist = pis.size();
        snpStore.setZero(pis.size());
        pll.setZero(pis.size());
        // --------------------------------------------------------------------------------
        // Scale the variances in each of the normal distributions by the genetic variance
        // and initialise the class membership probabilities
        // --------------------------------------------------------------------------------
        if (originalModel)
            gp = gamma * 0.01 * varg;
        else
            gp = gamma * sigmaSq;
    //    cout << varg << " " << gp.transpose() << endl;
        snpset.resize(ndist);
        for (unsigned k=0; k<ndist; ++k) {
            snpset[k].resize(0);
        }
        
        for (unsigned i=0; i<size; ++i) {
            // ------------------------------
            // Derived Bayes R implementation
            // ------------------------------
            // ----------------------------------------------------
            // Add back the content for the corrected rhs for SNP k
            // ----------------------------------------------------
            //my_rhs = Z.col(i).dot(ycorr);
            oldSample = values[i];
            rhs = ZPy[i] - n_gwas/n_ref*Z.col(i).dot(ghat) + ZPZdiag[i]*oldSample;
            // ------------------------------------------------------
            // Calculate the beta least squares updates and variances
            // ------------------------------------------------------
            b_ls = rhs / ZPZdiag[i];
            var_b_ls = gp.array() + vare / ZPZdiag[i];
            // ------------------------------------------------------
            // Calculate the likelihoods for each distribution
            // ------------------------------------------------------
            // ll  = (-1.0 / 2.0) * var_b_ls.array().log()  - (b_ls * b_ls)  / (2 * var_b_ls.array());
            ll = (-1.0 / 2.0) * var_b_ls.array().log()  - (b_ls * b_ls)  / (2 * var_b_ls.array()) + pis.array().log();
            // --------------------------------------------------------------
            // Calculate probability that snp is in each of the distributions
            // in this iteration
            // --------------------------------------------------------------
            // pll = (ll.array().exp().cwiseProduct(pis.array())) / ((ll.array().exp()).cwiseProduct(pis.array())).sum();
            for (unsigned k=0; k<pis.size(); ++k) {
                pll[k] = 1.0 / (exp(ll.array() - ll[k])).sum();
            }
            // --------------------------------------------------------------
            // Sample the group based on the calculated probabilities
            // --------------------------------------------------------------
            ssculm = 0.0;
            r = Stat::ranf();
            indistflag = 1;
            for (int kk = 0; kk < ndist; kk++)
            {
                ssculm += pll(kk);
                if (r < ssculm)
                {
                    indistflag = kk + 1;
                    snpStore(kk) = snpStore(kk) + 1;
                    break;
                }
            }
            snpset[indistflag-1].push_back(i);
            // --------------------------------------------------------------
            // Sample the effect given the group and adjust the rhs
            // --------------------------------------------------------------
            if (indistflag != 1)
            {
                v1 = ZPZdiag[i] + vare / gp((indistflag - 1));
                values[i] = normal.sample(rhs / v1, vare / v1);
                ghat  += Z.col(i) * (values[i] - oldSample);
                sumSq += (values[i] * values[i]) / gamma[indistflag - 1];
                ++numNonZeros;
            } else {
                if (oldSample) ghat -= Z.col(i) * oldSample;
                values[i] = 0.0;
            }
        }
}



void ApproxBayesR::SnpEffects::adjustByCG(const VectorXd &ZPy, const vector<SparseVector<double> > &ZPZsp, VectorXd &rcorr) {
    // construct mixed model equations for those SNPs with nonzero effects and solve the equations using conjugate gradient method
    // then adjust the Gibbs samples with the CG solutions
    
    VectorXd ZPyNZ(numNonZeros);

    vector<Triplet<double> > tripletList;
    tripletList.reserve(numNonZeros);
    
    for (unsigned i=0; i<numNonZeros; ++i) {
        unsigned row = deltaNzIdx[i];
        VectorXd val;
        val.setZero(size);
        for (SparseVector<double>::InnerIterator it(ZPZsp[row]); it; ++it) {
            val[it.index()] = it.value();
        }
        val[row] += lambdaVec[row];
//        cout << "val " << val.transpose() << endl;
        for (unsigned j=0; j<numNonZeros; ++j) {
            unsigned col = deltaNzIdx[j];
//            cout << i << " " << j << " " << row << " " << col << endl;
            tripletList.push_back(Triplet<double>(i, j, val[col]));
//            cout << i << " " << j << " " << val[col] << endl;
        }
        ZPyNZ[i] = ZPy[row];
    }
    
    SpMat C(numNonZeros, numNonZeros);
    C.setFromTriplets(tripletList.begin(), tripletList.end());
    C.makeCompressed();
    tripletList.clear();

//    cout << "C \n" << C.block(0,0,10,10) << endl;
    
    SimplicialLLT<SpMat> solverC;
    solverC.compute(C);
    
    if(solverC.info()!=Success) {
        cout << "Oh: Very bad" << endl;
    }
    
    SpMat eye(numNonZeros, numNonZeros);
    eye.setIdentity();
    
    SpMat Cinv = solverC.solve(eye);
    
    LLT<MatrixXd> llt;
    llt.compute(Cinv); // cholesky decomposition
    VectorXd nrnd(numNonZeros);
    for (unsigned i=0; i<numNonZeros; ++i) {
        nrnd[i] = Stat::snorm();
    }

//    ConjugateGradient<SpMat, Lower|Upper> cg;
//    cg.compute(C);
    VectorXd sol(numNonZeros);
//    sol = cg.solve(ZPyNZ + llt.matrixL()*nrnd);
    
    sol = Cinv * ZPyNZ + llt.matrixL()*nrnd;

//    cout << "numNonZeros " << numNonZeros << endl;
//    cout << "size C " << C.size() << endl;
////
//    cout << "ZPyNZ " << ZPyNZ << endl;
//        cout << "sol " << sol << endl;
//    cout << "uhatVec " << uhatVec << endl;
//
//    cout << "#nonZero:        " << numNonZeros << endl;
//    cout << "#iterations:     " << cg.iterations() << endl;
//    cout << "estimated error: " << cg.error()      << endl;

    double oldSample;
    for (unsigned i=0, j=0; i<size; ++i) {
        if (deltaNZ[i]) {
            oldSample = values[i];
            values[i] = sol[j];
            //values[i] *= sol[j]/uhatVec[i];
//            cout << i << " " << sol[j]/uhatVec[i] << endl;
            for (SparseVector<double>::InnerIterator it(ZPZsp[i]); it; ++it) {
                rcorr[it.index()] += it.value() * (oldSample - values[i]);
            }
            ++j;
        }
    }
    
//    cout << "old sumsq " << sumSq << endl;
    sumSq = values.cwiseProduct(invGammaVec).dot(values);
//    cout << "new sumsq " << sumSq << endl;

}

void ApproxBayesR::SnpEffects::sampleFromFC(const VectorXd &ZPy, const SpMat &ZPZsp, const VectorXd &ZPZdiag,
                                            VectorXd &rcorr, const VectorXd &LDsamplVar,
                                            const double sigmaSq, const VectorXd &pis, const VectorXd &gamma, VectorXd &snpStore,
                                            const double varg, const double vare, const double ps, const double overdispersion, const bool originalModel) {
    // CG-accelerated Gibbs sampling algorithm
    // first sample delta conditional on beta for all SNPs
    // then construct mixed model equations for which the solutions are samples from the Gibbs sampling
    // and solve the equations by conjugate gradient method
    
    VectorXd lambdaVec(size);
    VectorXd invGammaVec(size);
    
    unsigned ndist = gamma.size();
    snpStore.setZero(ndist);
    
    double varei;
    double rhs;
    
    ArrayXd wtdSigmaSq(ndist);
    ArrayXd invWtdSigmaSq(ndist);
    ArrayXd logWtdSigmaSq(ndist);
    ArrayXd logPis = pis.array().log();
    ArrayXd invLhs(ndist);
    ArrayXd uhat(ndist);
    ArrayXd logDelta(ndist);
    ArrayXd probDelta(ndist);
    
    unsigned delta;
    
    if (originalModel) {
        wtdSigmaSq = gamma * 0.01 * varg;
    } else {
        wtdSigmaSq = gamma * sigmaSq;
    }
    
    invWtdSigmaSq = wtdSigmaSq.inverse();
    logWtdSigmaSq = wtdSigmaSq.log();
    
    VectorXd invGamma = gamma.inverse();
    invGamma[0] = 0;


    for (unsigned i=0; i<size; ++i) {
        
        varei = LDsamplVar[i]*varg + vare + ps + overdispersion;
        
        rhs  = rcorr[i] + ZPZdiag[i]*values[i];
        
        invLhs = (ZPZdiag[i] + varei*invWtdSigmaSq).inverse();
        uhat = invLhs*rhs;
        
        logDelta = 0.5*(invLhs.log() - logWtdSigmaSq + uhat*rhs) + logPis;
        logDelta[0] = logPis[0];
        
        for (unsigned k=0; k<ndist; ++k) {
            probDelta[k] = 1.0f/(logDelta-logDelta[k]).exp().sum();
        }
        
        delta = bernoulli.sample(probDelta);
        
        deltaNZ[i] = delta ? 1:0;
        
        snpset[delta].push_back(i);
        snpStore[delta]++;

        lambdaVec[i] = varei*invWtdSigmaSq[delta];
        invGammaVec[i] = invGamma[delta];
    }
    
    numNonZeros = deltaNZ.sum();
    
    VectorXd lambdaNZ(numNonZeros);
    VectorXd RHS(numNonZeros);
    SpMat eye(numNonZeros, numNonZeros);
    vector<Triplet<double> > tripletList;
    tripletList.reserve(numNonZeros);
    for (unsigned i=0, j=0; i<size; ++i) {
        if (deltaNZ[i]) {
            tripletList.push_back(Triplet<double>(i,i,1));
            lambdaNZ[j] = lambdaVec[i];
            RHS[j] = ZPy[i] + normal.sample(0.0, ZPZdiag[i] + lambdaVec[i]);
            ++j;
        }
    }
    eye.setFromTriplets(tripletList.begin(), tripletList.end());
    eye.makeCompressed();
    tripletList.clear();

    SpMat LHS = eye * ZPZsp * eye;
    LHS.diagonal() += lambdaNZ;
    
    ConjugateGradient<SpMat, Lower|Upper> cg;
    cg.compute(LHS);
    values = cg.solve(RHS);
    
    sumSq = values.cwiseProduct(invGamma).dot(values);
    
    rcorr = ZPy - ZPZsp * values;
}

void ApproxBayesR::SnpEffects::sampleFromFC(vector<VectorXd> &wcorrBlocks, const vector<MatrixXd> &Qblocks, vector<VectorXd> &whatBlocks,
                                            const vector<LDBlockInfo*> keptLdBlockInfoVec, const VectorXd &nGWASblocks, const VectorXd &vareBlocks,
                                            const double sigmaSq, const VectorXd &pis, const VectorXd &gamma, VectorXd &snpStore, const double varg,
                                            const bool originalModel) {
    // -----------------------------------------
    // This method uses low-rank model with eigen-decomposition of LD matrices
    // -----------------------------------------
    long nBlocks = keptLdBlockInfoVec.size();
    
    whatBlocks.resize(nBlocks);
    for (unsigned i=0; i<nBlocks; ++i) {
        whatBlocks[i].resize(wcorrBlocks[i].size());
    }

    double ssq[nBlocks], s2pq[nBlocks], nnz[nBlocks];
    memset(ssq,0, sizeof(double)*nBlocks);
    memset(s2pq,0,sizeof(double)*nBlocks);
    memset(nnz,0, sizeof(double)*nBlocks);
    
    double *valuesPtr = values.data(); // for openmp, otherwise when one thread writes to the vector, the vector locking prevents the writing from other threads

    vector<double> urnd(size), nrnd(size);
    for (unsigned i=0; i<size; ++i) { // need this for openmp to work
        urnd[i] = Stat::ranf();
        nrnd[i] = Stat::snorm();
    }
    
    // R specific parameters
    int ndist = pis.size();
    ArrayXd logPis = pis.array().log();
    ArrayXd wtdSigmaSq(ndist);
    ArrayXd invWtdSigmaSq(ndist);
    ArrayXd logWtdSigmaSq(ndist);

    if (originalModel) {
        wtdSigmaSq = gamma * 0.01 * varg;
    } else {
        wtdSigmaSq = gamma * sigmaSq;
    }
    
    invWtdSigmaSq = wtdSigmaSq.inverse();
    logWtdSigmaSq = wtdSigmaSq.log();

    vector<vector<vector<unsigned> > > snpsetBlocks(nBlocks);
    for (unsigned i=0; i<nBlocks; ++i) {
        snpsetBlocks[i].resize(ndist);
        for (unsigned k=0; k<ndist; ++k) {
            snpsetBlocks[i][k].resize(0);
        }
    }

    // --------------------------------------------------------------------------------
    // Cycle over all variants in the window and sample the genetics effects
    // --------------------------------------------------------------------------------

    #pragma omp parallel for schedule(dynamic)
    for(unsigned blk = 0; blk < nBlocks; blk++){
        Ref<const MatrixXd> Q = Qblocks[blk];
        Ref<VectorXd> wcorr = wcorrBlocks[blk];
        Ref<VectorXd> what = whatBlocks[blk];

        what.setZero();
        
        LDBlockInfo *blockInfo = keptLdBlockInfoVec[blk];
        
        unsigned blockStart = blockInfo->startSnpIdx;
        unsigned blockEnd   = blockInfo->endSnpIdx;
        
        double invVareDn = nGWASblocks[blk] / vareBlocks[blk];

        ArrayXd invLhs = 1.0/(invVareDn + invWtdSigmaSq);
        ArrayXd logInvLhsMsigma = invLhs.log() - logWtdSigmaSq;

        for(unsigned i = blockStart; i <= blockEnd; i++){
            double oldSample = valuesPtr[i];
            Ref<const VectorXd> Qi = Q.col(i - blockStart);
            double rhs = (Qi.dot(wcorr) + oldSample)*invVareDn;
            ArrayXd uhat = invLhs * rhs;
            ArrayXd logDelta = 0.5*(logInvLhsMsigma + uhat*rhs) + logPis;
            logDelta[0] = logPis[0];
            
            ArrayXd probDelta(ndist);
            for (unsigned k=0; k<ndist; ++k) {
                probDelta[k] = 1.0f/(logDelta-logDelta[k]).exp().sum();
            }
                        
//            #pragma omp critical
//            {
            unsigned delta = bernoulli.sample(probDelta, urnd[i]);
            snpsetBlocks[blk][delta].push_back(i);
//            }
            
            if (delta) {
                valuesPtr[i] = uhat[delta] + nrnd[i]*sqrtf(invLhs[delta]);
                wcorr += Qi*(oldSample - valuesPtr[i]);
                what  += Qi* valuesPtr[i];
                ssq[blk] += (valuesPtr[i] * valuesPtr[i]) / gamma[delta];
                ++nnz[blk];
            }
            else {
                if (oldSample) wcorr += Qi * oldSample;
                valuesPtr[i] = 0.0;
            }
        }

    }

    // ---------------------------------------------------------------------
    // Tally up the effect sum of squares and the number of non-zero effects
    // ---------------------------------------------------------------------
    sumSq = 0.0;
    numNonZeros = 0;
    nnzPerBlk.setZero(nBlocks);
    ssqBlocks.setZero(nBlocks);
    snpStore.setZero(ndist);
    snpset.resize(ndist);
    for (unsigned k=0; k<ndist; ++k) {
        snpset[k].resize(0);
    }
    for (unsigned blk=0; blk<nBlocks; ++blk) {
        sumSq += ssq[blk];
        numNonZeros += nnz[blk];
        nnzPerBlk[blk] = nnz[blk];
        ssqBlocks[blk] = ssq[blk];
        for (unsigned k=0; k<ndist; ++k) {
            for (unsigned j=0; j<snpsetBlocks[blk][k].size(); ++j) {
                snpset[k].push_back(snpsetBlocks[blk][k][j]);
                snpStore[k]++;
            }
        }
    }
    values = VectorXd::Map(valuesPtr, size);

}

