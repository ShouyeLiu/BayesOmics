// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef BAYESOMICS_PARALLEL_GIBBS_HPP
#define BAYESOMICS_PARALLEL_GIBBS_HPP

#include <Eigen/Core>
#include "GenotypeView.hpp"
#include <algorithm>
#include <omp.h>
#include <vector>
#include <stdexcept>

// Parallelize observations, retaining the sequential Gibbs order over SNPs.
// Fixed chunks and ordered sums give the same arithmetic at every thread count.
namespace ParallelGibbs {
constexpr Eigen::Index chunkSize=8192;
constexpr Eigen::Index threshold=32768;
inline int workers(Eigen::Index tasks) {return int(std::max<Eigen::Index>(1,std::min<Eigen::Index>(tasks,omp_get_max_threads())));}

// Large teams are inefficient for short per-SNP operations. Keep the validated
// four-worker scale at ~100K observations; expand teams for much longer vectors.
inline int observationWorkers(Eigen::Index chunks) { return workers(std::min(chunks,std::max<Eigen::Index>(4,chunks/4))); }

inline double dot(Eigen::Ref<const Eigen::VectorXd> x,Eigen::Ref<const Eigen::VectorXd> y) {
    if(x.size()!=y.size())throw std::invalid_argument("Gibbs dot dimension mismatch");
    if(x.size()<threshold)return x.dot(y);
    const Eigen::Index chunks=(x.size()+chunkSize-1)/chunkSize;
    std::vector<double> sums(chunks);
    #pragma omp parallel for schedule(static) num_threads(observationWorkers(chunks)) if(!omp_in_parallel() && omp_get_max_threads()>1)
    for(Eigen::Index c=0;c<chunks;++c) {
        const Eigen::Index start=c*chunkSize,count=std::min(chunkSize,x.size()-start);
        sums[c]=x.segment(start,count).dot(y.segment(start,count));
    }
    double total=0;for(double part:sums)total+=part;return total;
}

inline void axpy(Eigen::Ref<const Eigen::VectorXd> x,double coefficient,Eigen::Ref<Eigen::VectorXd> y) {
    if(x.size()!=y.size())throw std::invalid_argument("Gibbs residual dimension mismatch");
    if(coefficient==0)return;
    if(x.size()<threshold){y.noalias()+=x*coefficient;return;}
    const Eigen::Index chunks=(x.size()+chunkSize-1)/chunkSize;
    #pragma omp parallel for schedule(static) num_threads(observationWorkers(chunks)) if(!omp_in_parallel() && omp_get_max_threads()>1)
    for(Eigen::Index c=0;c<chunks;++c) {
        const Eigen::Index start=c*chunkSize,count=std::min(chunkSize,x.size()-start);
        y.segment(start,count).noalias()+=x.segment(start,count)*coefficient;
    }
}

inline void columnNorms(const Eigen::MatrixXd &matrix,Eigen::VectorXd &norms) {
    norms.resize(matrix.cols());
    #pragma omp parallel for schedule(static) num_threads(workers(matrix.cols())) if(matrix.size()>=threshold && !omp_in_parallel())
    for(Eigen::Index j=0;j<matrix.cols();++j)norms[j]=dot(matrix.col(j),matrix.col(j));
}

inline void sparseProduct(const Eigen::MatrixXd &matrix,const Eigen::VectorXd &coefficients,Eigen::VectorXd &prediction) {
    if(matrix.cols()!=coefficients.size())throw std::invalid_argument("Gibbs prediction dimension mismatch");
    std::vector<Eigen::Index> active;
    for(Eigen::Index j=0;j<coefficients.size();++j)if(coefficients[j]!=0)active.push_back(j);
    prediction.setZero(matrix.rows());
    const Eigen::Index chunks=(matrix.rows()+chunkSize-1)/chunkSize;
    #pragma omp parallel for schedule(static) num_threads(observationWorkers(chunks)) if(matrix.rows()>=threshold && !omp_in_parallel() && omp_get_max_threads()>1)
    for(Eigen::Index c=0;c<chunks;++c) {
        const Eigen::Index start=c*chunkSize,count=std::min(chunkSize,matrix.rows()-start);
        auto out=prediction.segment(start,count);
        for(auto j:active)out.noalias()+=matrix.col(j).segment(start,count)*coefficients[j];
    }
}
inline void columnNorms(const GenotypeView &matrix,Eigen::VectorXd &norms) {
    if(!matrix.streamed()){columnNorms(matrix.dense(),norms);return;}
    norms.resize(matrix.cols());
    for(Eigen::Index j=0;j<matrix.cols();++j){const auto x=matrix.col(j);norms[j]=dot(x,x);}
}
inline void sparseProduct(const GenotypeView &matrix,const Eigen::VectorXd &coefficients,Eigen::VectorXd &prediction) {
    if(!matrix.streamed()){sparseProduct(matrix.dense(),coefficients,prediction);return;}
    if(matrix.cols()!=coefficients.size())throw std::invalid_argument("Gibbs prediction dimension mismatch");
    prediction.setZero(matrix.rows());
    for(Eigen::Index j=0;j<coefficients.size();++j)if(coefficients[j]!=0)axpy(matrix.col(j),coefficients[j],prediction);
}

}
#endif
