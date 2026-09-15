// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "ParallelGibbs.hpp"
#include "MemoryBudget.hpp"
#include <exception>

// Exact arithmetic cache for sequential individual-level Gibbs sampling.
// At each panel boundary s = X_panel' r. After r += x_j*d, s += X_panel'x_j*d.
// Full individual residuals are still updated after every draw; correlations
// across panels are included when the next panel starts. No LD approximation,
// mixture distribution, random-number generator or model state lives here.
class PanelScores {
    bool ready=false,active=false;
    Eigen::Index width=64,start=0,count=0;
    std::vector<Eigen::MatrixXd> gram;
    Eigen::VectorXd scores;
public:
    bool enabled() const{return active;}
    Eigen::Index panelWidth() const{return active?width:0;}
    void begin(Eigen::Index j,const GenotypeView &x,const Eigen::VectorXd &residual,
               const Eigen::VectorXd &norms) {
        if(!ready) {
            ready=true;
            // Streamed BED keeps its existing bounded, single-column path.
            if(x.streamed() || x.cols()<2 || x.rows()<4096)return;
            width=std::min<Eigen::Index>(width,x.cols());
            if(MemoryBudget::enabled()) {
                const auto perColumn=MemoryBudget::bytes(x.cols(),1);
                width=std::min<Eigen::Index>(width,MemoryBudget::status().available/16/std::max<uint64_t>(1,perColumn));
            }
            if(width<2)return;
            MemoryBudget::require(MemoryBudget::bytes(x.cols(),width),"Bounded individual cross-product panels");
            const Eigen::Index panels=(x.cols()+width-1)/width;
            gram.resize(panels);std::vector<std::exception_ptr> errors(panels);
            const auto &matrix=x.dense();
            #pragma omp parallel for schedule(dynamic) num_threads(ParallelGibbs::workers(panels)) if(panels>1)
            for(Eigen::Index p=0;p<panels;++p) {
                try {
                    const Eigen::Index begin=p*width,n=std::min(width,x.cols()-begin);
                    const auto panel=matrix.middleCols(begin,n);
                    gram[p].noalias()=panel.transpose()*panel;
                    // Match the established fixed-chunk diagonal arithmetic.
                    gram[p].diagonal()=norms.segment(begin,n);
                } catch(...) {errors[p]=std::current_exception();}
            }
            for(const auto &error:errors)if(error)std::rethrow_exception(error);
            scores.resize(width);active=true;
        }
        if(!active || j%width)return;
        start=j;count=std::min(width,x.cols()-start);
        const auto &matrix=x.dense();
        #pragma omp parallel for schedule(static) num_threads(ParallelGibbs::workers(count)) if(count>1)
        for(Eigen::Index k=0;k<count;++k)scores[k]=ParallelGibbs::dot(matrix.col(start+k),residual);
    }
    double score(Eigen::Index j) const{return scores[j-start];}
    void update(Eigen::Index j,double change) {
        if(!active || change==0)return;
        const Eigen::Index local=j-start,remaining=count-local;
        scores.segment(local,remaining).noalias()+=gram[start/width].col(local).segment(local,remaining)*change;
    }
};
