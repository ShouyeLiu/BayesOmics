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


#ifndef stat_hpp
#define stat_hpp

#include <stdio.h>
//#include <random>
#include <boost/random.hpp>
#include <boost/math/distributions.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/inverse_chi_squared.hpp>
#include <Eigen/Eigen>

#include "Logger.hpp"

using namespace Eigen;
//using namespace std;


namespace Stat {
  
    typedef boost::mt19937 random_engine;
    typedef boost::uniform_01<> uniform_01;
    typedef boost::random::uniform_real_distribution<> uniform_real_distribution;
    typedef boost::normal_distribution<> normal_distribution;
    typedef boost::gamma_distribution<> gamma_distribution;
    typedef boost::math::inverse_chi_squared_distribution<> inverse_chi_squared_distribution;
    typedef boost::math::beta_distribution<> beta_distribution;
    
    typedef boost::variate_generator<random_engine&, uniform_01> uniform01_generator;
    typedef boost::variate_generator<random_engine&, normal_distribution> normal_generator;
    typedef boost::variate_generator<random_engine&, gamma_distribution> gamma_generator;
    
    static random_engine engine;
    static uniform01_generator ranf(engine, uniform_01());
    static normal_generator snorm(engine, normal_distribution(0,1));  // standard normal
    
    void seedEngine(const int seed);
    // various distribution
    // Uniform class definition
    class Uniform {
    public:
        // Static method to generate a uniform random number in a specified range [lower, upper]
        static double uniform(double lower, double upper) {
            boost::random::uniform_real_distribution<> dist(lower, upper);
            boost::variate_generator<random_engine&, boost::random::uniform_real_distribution<>> gen(engine, dist);
            return gen();
        }
    };

    class Normal {
    public:
        boost::math::normal_distribution <> d;
        double SQRT_2;
        double Inv_SQRT_2;
        
        Normal(){
            d = boost::math::normal_distribution <> (0 ,1);
            SQRT_2 = sqrt(2);
            Inv_SQRT_2 = 1.0/SQRT_2;
        }
        
        double sample(const double mean, const double variance);
        double cdf_01(const double value);
        double quantile_01(const double value);
    };
    
    class Flat : public Normal {
    public:
        // flat prior is a normal with infinit variance
    };

    class ChiSq {
        public:
        double sample(const double df,const double scale);
        double sample(double df);
        static double qchisq(double percentile, double df);
        static double pchisq(double chiSquStat, double df);
    };
    
    class InvChiSq {
    public:
        double sample(const double df, const double scale);
    };
    
    class Gamma {
    public:
        double sample(const double shape, const double scale);
    };
    
    class Beta {
    public:
        double sample(const double a, const double b);
    };
    
    class Bernoulli {
    public:
        unsigned sample(const double p);
        unsigned sample(const VectorXd &p); // multivariate sampling, return the index of component.
        unsigned sample(const VectorXd &p, const double rnd); // sampling with a given random number.
    };

    class Dirichlet {
    public:
        Gamma gamma;
        VectorXd sample(const int n, const VectorXd &irx);
    };
    
    class NormalZeroMixture {
    public:
        Normal normal;
        Bernoulli bernoulli;
        
        double sample(const double mean, const double variance, const double p);
    };
    
    class MixtureNormals {
    public:
        
    };

    class TruncatedNormal : public Normal {
    public:
        
        double sample_tail_01_rejection(const double a);
        double sample_lower_truncated(const double mean, const double sd, const double a);  // a < x < inf
        double sample_upper_truncated(const double mean, const double sd, const double b);  // -inf < x < b
    };

    class MultiNormal : public Normal {
        public:
        
        Normal normal;

        template <typename T1> void triMultLX(Ref<MatrixXd> Y, const Eigen::MatrixBase<T1>& L,const Ref<const MatrixXd>& X);
        Eigen::MatrixXd sample(int N, Eigen::MatrixXd mu,Eigen::MatrixXd Sigma);
    };

    class WishartBase {
    public:
        // storage
        int q_;
        LLT<MatrixXd> cholX_;
        LLT<MatrixXd> cholPsi_;
        MatrixXd Z_;
        MatrixXd XL_;
        Normal normal;
        ChiSq  chisq;
    public:
        /// Constructor
        void assignInitValue(int q);
        /// Random draw with scale matrix input
        void GenerateLowerTri(Ref<MatrixXd> XL,const Ref<const MatrixXd>& PsiL, double nu); 
        /// Random draw with precision matrix input
        void GenerateLowerTriXi(Ref<MatrixXd> XL,const Ref<const MatrixXd>& XiL, double nu); 
        template <typename T1> void triMultLiX(Ref<MatrixXd> X, const Eigen::MatrixBase<T1>& L); 
        template <typename T1, typename T2> void triMultLL(const Eigen::MatrixBase<T1>& X, const Ref<const MatrixXd>& L1, const Eigen::MatrixBase<T2>& L2); 
  };  //  eend of wshart

    class Wishart {
    public:
        WishartBase wishBase;
        Eigen::MatrixXd sample(int N, Eigen::MatrixXd Psi, Eigen::VectorXd nu,bool inverse = false);
        void ReverseCholesky(Ref<MatrixXd> L,const Ref<const MatrixXd>& V,LLT<MatrixXd>& llt);
        template <typename T1> void InverseLLt(Ref<MatrixXd> X, const Eigen::MatrixBase<T1>& L,Ref<MatrixXd> L2,Ref<MatrixXd> U,const Ref<const MatrixXd>& I);
        template <typename T1, typename T2> void CrossProdLLt(const Eigen::MatrixBase<T1>& X, const Ref<const MatrixXd>& L, const Eigen::MatrixBase<T2>& U);
        template <typename T1, typename T2> void CrossProdLtL(const Eigen::MatrixBase<T1>& X,const Ref<const MatrixXd>& L,const Eigen::MatrixBase<T2>& U); 

        // arma::mat riwish_C(int v, arma::mat S);
    };

} // end of Stat

#endif /* stat_hpp */
