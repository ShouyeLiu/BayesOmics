// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef BAYESOMICS_OMICS_MIXTURE_HPP
#define BAYESOMICS_OMICS_MIXTURE_HPP

#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>
#include <vector>

class Data;

// Conditional integration for the CO-EIEO spike-and-slab prior.
// EIEO's latent effect is drawn even for the zero component; only sqrt(gamma)*u
// enters the phenotype likelihood. Residual covariance is zero.
namespace OmicsMixture {
struct Draws {
    std::vector<std::size_t> offset;
    Eigen::VectorXd uniform, normal;
    explicit Draws(const Data &data);
};
struct Conditional {
    Eigen::VectorXd mean, variance, probability;
};
inline Conditional eieo(double precision, double score, double priorPrecision,
                        double crossPrecision, double otherLatent,
                        const Eigen::VectorXd &gamma, const Eigen::VectorXd &pi) {
    if (gamma.size()!=pi.size() || gamma.size()<2 || gamma[0]!=0 ||
        (gamma.array()<0).any() || (pi.array()<=0).any() ||
        std::abs(pi.sum()-1)>1e-8 || precision<0 || priorPrecision<=0)
        throw std::invalid_argument("Invalid EIEO conditional mixture parameters");
    Conditional out;
    out.variance=(gamma.array()*precision+priorPrecision).inverse();
    Eigen::VectorXd rhs=gamma.array().sqrt()*score-crossPrecision*otherLatent;
    out.mean=out.variance.array()*rhs.array();
    Eigen::VectorXd logWeight=pi.array().log()+0.5*(out.variance.array().log()+out.mean.array()*rhs.array());
    out.probability=(logWeight.array()-logWeight.maxCoeff()).exp();
    out.probability/=out.probability.sum();
    return out;
}
inline unsigned choose(const Eigen::VectorXd &probability, double uniform) {
    double cumulative=0;
    for (int k=0;k<probability.size()-1;++k) {
        cumulative+=probability[k];
        if (uniform<cumulative) return k;
    }
    return probability.size()-1;
}
}
#endif
