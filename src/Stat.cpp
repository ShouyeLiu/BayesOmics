//
//  stat.cpp
//  gctb
//
//  Created by Jian Zeng on 14/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#include "Stat.hpp"
#include "Logger.hpp"

void Stat::seedEngine(const int seed){
    if (seed) {
        srand(seed);
        engine.seed(seed);
    } else {
        srand((int)time(NULL));
        engine.seed((int)time(NULL));
    }
    LOGGER << "Seed is set as " << std::to_string(seed) << "." << std::endl;
}

double Stat::Normal::sample(const double mean, const double variance){
    return mean + snorm()*sqrtf(variance);
}

double Stat::Normal::cdf_01(const double value){
    return cdf(d, value);
    //return 0.5 * boost::math::erfc(-value * Inv_SQRT_2);
}

double Stat::Normal::quantile_01(const double value){
    return quantile(d, value);
    //return - SQRT_2 * boost::math::erfc_inv(2.0 * value);
}

double Stat::ChiSq::sample(const double df, const double scale){
	gamma_generator sgamma(engine, gamma_distribution(0.5f*df, 1));
	return scale * (2.0f*sgamma());

}
double Stat::ChiSq::qchisq(double percentile, double df){
	boost::math::chi_squared dist(df);
	double quantileValue = boost::math::quantile(dist, percentile);
	return quantileValue;
}
double Stat::ChiSq::sample(double df){
	boost::random::chi_squared_distribution<> chiSquared(df);
	return chiSquared(engine);
}

double Stat::ChiSq::pchisq(double chiSquStat, double df){
	boost::math::chi_squared dist(df);
	double pValue = 1 - boost::math::cdf(dist, chiSquStat);
	return pValue;
}



double Stat::InvChiSq::sample(const double df, const double scale){
    // inverse_chi_squared_distribution invchisq(df, scale);
    // return boost::math::quantile(invchisq, ranf());   // don't know why this is not correct
    if(df <= 0) LOGGER.e(0, "df should be greater than 0 in InvChiSq::sample");
	if(scale <= 0) LOGGER.e(0, "scale should be greater than 0 in InvChiSq::sample");
    gamma_generator sgamma(engine, gamma_distribution(0.5f*df, 1));
    return scale/(2.0f*sgamma());
}

double Stat::Gamma::sample(const double shape, const double scale){
    gamma_generator sgamma(engine, gamma_distribution(shape, scale));
    return sgamma();
}

double Stat::Beta::sample(const double a, const double b){
    beta_distribution beta(a,b);
    return boost::math::quantile(beta,ranf());
	// return boost::math::ibeta_inv(a, b, ranf()); 
}

unsigned Stat::Bernoulli::sample(const double p){
    if (std::isnan(p)) return ranf() < 0.5 ? 1:0;
    else return ranf() < p ? 1:0;
}

unsigned Stat::Bernoulli::sample(const VectorXd &p){
    double cum = 0;
    double rnd = ranf();
    long size = p.size();
    unsigned ret = 0;
    for (unsigned i=0; i<size; ++i) {
        if (!std::isnan(p[i])) cum += p[i];
        if (rnd < cum) {
            ret = i;
            break;
        }
    }
    return ret;
}

unsigned Stat::Bernoulli::sample(const VectorXd &p, const double rnd){
    // sampler with a given random number
    double cum = 0;
    //double rnd = ranf();
    long size = p.size();
    unsigned ret = 0;
    for (unsigned i=0; i<size; ++i) {
        if (!std::isnan(p[i])) cum += p[i];
        if (rnd < cum) {
            ret = i;
            break;
        }
    }
    return ret;
}

double Stat::NormalZeroMixture::sample(const double mean, const double variance, const double p){
    return bernoulli.sample(p) ? normal.sample(mean, variance) : 0;
}

// Sample Dirichlet

VectorXd Stat::Dirichlet::sample(const int n, const VectorXd &irx){
    VectorXd ps(n);
    double sx = 0.0;
    for (int i = 0; i < n; i++)
    {
        ps[i] = gamma.sample(irx(i), 1.0);
        sx += ps[i];
    }
    ps = ps / sx;
    return ps;
}

double Stat::TruncatedNormal::sample_tail_01_rejection(const double a){  // a < x < inf
    double b = 0.5 * a*a;
    double w, v;
    do {
        double u = ranf();
        w = b - log(u);
        v = ranf();
    } while (v > w/b);
    return sqrt(2.0 * w);
}

double Stat::TruncatedNormal::sample_lower_truncated(const double mean, const double sd, const double a){  // a < x < inf
//    int seed = rand();
//    double x = truncated_normal_a_sample (mean, sd, a, seed);
//    return x;
    
    if (a - mean > 5*sd) {
        return mean + sd * sample_tail_01_rejection((a - mean)/sd);
    } else {
        double alpha = (a - mean)/sd;
        double alpha_cdf = cdf_01(alpha);
        double u = ranf();
        double x = alpha_cdf + (1.0 - alpha_cdf) * u;  // x ~ Uniform(alpha_cdf, 1)
        if (x <= 0 || x >= 1) LOGGER << "alpha " << alpha << " alpha_cdf " << alpha_cdf << " u " << u << " x " << x << std::endl;
        return mean + sd * quantile_01(x);
    }
}

double Stat::TruncatedNormal::sample_upper_truncated(const double mean, const double sd, const double b){  // -inf < x < b
//    int seed = rand();
//    double x = truncated_normal_b_sample (mean, sd, b, seed);
//    return x;

    if (mean - b > 5*sd) {
        return mean - sd * sample_tail_01_rejection((mean - b)/sd);
    } else{
        double beta = (b - mean)/sd;
        double beta_cdf = cdf_01(beta);
        double u;
        do {
            u = ranf();
        } while (!u);
        double x = beta_cdf * u;  // x ~ Uniform(0, beta_cdf);
        if (x <= 0 || x >= 1) LOGGER << "beta " << beta << " beta_cdf " << beta_cdf << " u " << u << " x " << x << std::endl;
        return mean + sd * quantile_01(x);
    }
}


////////// MVN
/// Generate a random sample from the Multivariate Normal distribution.
///
/// Generate `N` independent draws from a `q` dimensional Multivariate Normal distribution.  Each argument can be vectorized, meaning that it can have length either `N` or `1`, denoted here as `n`.
///
/// @param [in] N Integer number of random draws.
/// @param [in] mu Matrix of size `q x n` of mean vectors.
/// @param [in] Sigma Matrix of size `q x nq` of variance matrices.
///
/// @return Matrix of size `q x N` of random draws.
//[[Rcpp::export]]
Eigen::MatrixXd Stat::MultiNormal::sample(int N, Eigen::MatrixXd mu,
					   Eigen::MatrixXd Sigma) {
  int q = mu.rows();
  bool singleMu = (mu.cols() == 1);
  bool singleSigma = (Sigma.cols() == q);
  int ii,jj;
  // output variables
  MatrixXd X(q,N);
  // internal variables
  LLT<MatrixXd> lltq(q);
  MatrixXd SigmaL(q,q);
  VectorXd lambda(q);
  MatrixXd Z(q,N);
  // generate iid normals
  for(ii=0; ii<N; ii++) {
    for(jj=0; jj<q; jj++) {
    //   Z(jj,ii) = norm_rand();
    Z(jj,ii) = normal.sample(0,1);
    }
  }
  // multiply by Sigma^{1/2}
  if(singleSigma) {
    lltq.compute(Sigma);
    SigmaL = lltq.matrixL();
    triMultLX(X, SigmaL, Z);
  } else {
    for(ii=0; ii<N; ii++) {
      lltq.compute(Sigma.block(0,ii*q,q,q));
      SigmaL = lltq.matrixL();
      triMultLX(X.col(ii), SigmaL, Z.col(ii));
    }
  }
  // add mu
  if(singleMu) {
    X.colwise() += mu.col(0);
  } else {
    X += mu;
  }
  return X;
}


  /// @brief Left-multiplication by a lower triangular matrix
  ///
  /// Performs the matrix multiplication `Y = L * X`, where `L` is a lower triangular matrix.
  ///
  /// @param [in] X Matrix of size `n x p` to be multiplied.
  /// @param [in] L Lower triangular matrix of size `n x n` left-multiplying `X`.
  /// @param [out] Y Matrix of size `n x p` containing the product of `L` and `X`.
  ///
  /// @note Does not work properly if any of the inputs arguments are also outputs.
  /// @note For safest results, all "off-triangular" elements of triangular arguments should be set to zero.
  template <typename T1>
  void Stat::MultiNormal::triMultLX(Ref<MatrixXd> Y, const Eigen::MatrixBase<T1>& L,
		 const Ref<const MatrixXd>& X) {
    Y.noalias() = L.template triangularView<Eigen::Lower>() * X;
    return;
  }

/////////////////////////// IW //////////////////////////

inline void Stat::WishartBase::assignInitValue(int q) {
	q_ = q;
	cholX_.compute(MatrixXd::Identity(q_,q_));
	cholPsi_.compute(MatrixXd::Identity(q_,q_));
	Z_ = MatrixXd::Zero(q_,q_);
	XL_ = MatrixXd::Zero(q_,q_);
}
/// Simulate a random draw from the lower Cholesky factor of a WishartBase distribution with given scale matrix and shape parameters.
///
/// @param [out] XL Matrix of size `q x q` containing the random draw, which lives only on the lower triangular half of the matrix.
/// @param [in] PsiL Matrix of size `q x q` containing the lower Cholesky factor of the scale matrix parameter `Psi`.
/// @param [in] nu Shape parameter.
inline void Stat::WishartBase::GenerateLowerTri(Ref<MatrixXd> XL,
												const Ref<const MatrixXd>& PsiL,
												double nu) {
	int ii, jj;
	for(ii=0; ii<q_; ii++) {
		// diagonal
		// XL_(ii,ii) = sqrt(chisq_rand(nu-ii));
		XL(ii,ii) = sqrt(chisq.sample(nu-ii));
		// off-diagonals
		for(jj=0; jj<ii; jj++) {
			// XL_(ii,jj) = norm_rand();
			XL_(ii,jj) = normal.sample(0,1);
		}
	}
	// multiply by precision matrix XL = PsiL * XL
	triMultLL(XL, PsiL, XL_);
	return;
}

/// Simulate a random draw from the lower Cholesky factor of a WishartBase distribution with given precision matrix and shape parameters. 
///
/// @param [out] XL Matrix of size `q x q` containing the random draw, which lives only on the lower triangular half of the matrix.
/// @param [in] XiL Matrix containing the inverse of the lower Cholesky factor of the scale matrix parameter `Psi`, namely `Psi^{-1} = XiL' * XiL`.
/// @param [in] nu Shape parameter.
inline void Stat::WishartBase::GenerateLowerTriXi(Ref<MatrixXd> XL,
												  const Ref<const MatrixXd>& XiL,
												  double nu) {
	int ii, jj;
	for(ii=0; ii<q_; ii++) {
		// diagonal
		// XL(ii,ii) = sqrt(chisq_rand(nu-ii));
		XL(ii,ii) = sqrt(chisq.sample(nu-ii));
		// off-diagonals
		for(jj=0; jj<ii; jj++) {
			// XL(ii,jj) = norm_rand();
			XL(ii,jj) =normal.sample(0,1);
		}
	}
	// multiply by precision matrix
	triMultLiX(XL, XiL);
	return;
}

/// Generate a random sample from the WishartBase or Inverse-WishartBase distribution.
///
/// Generate `N` independent draws from a `q x q` WishartBase or Inverse-WishartBase.  Each argument except `inverse` can be vectorized, meaning that it can have length either `N` or `1`, denoted here as `n`.
///
/// @param [in] N Integer number of random draws to produce
/// @param [in] Psi Matrix of size `q x nq` of scale matrices.
/// @param [in] nu Vector of `n` degrees-of-freedom parameters.
/// @param [in] inverse Whether to evaluate the log-density of the WishartBase or Inverse-WishartBase distribution.
///
/// @return A `q x Nq` matrix of `N` random draws.
//[[Rcpp::export]]


/// @brief Multiplication of two lower triangular matrices
///
/// Performs the matrix multiplication `X = L1 * L2`, where `L1` and `L2` are lower triangular matrices.
///
/// @param [in] L1 First `n x n` lower triangular matrix.
/// @param [in] L2 Second `n x n` lower triangular matrix.
/// @param [out] X Matrix of size `n x n` containing the product of `L1` and `L2`.
///
/// @note Does not work properly if any of the inputs arguments are also outputs.
/// @note For safest results, all "off-triangular" elements of triangular arguments should be set to zero.
template <typename T1, typename T2>
void Stat::WishartBase::triMultLL(const Eigen::MatrixBase<T1>& X,
								  const Ref<const MatrixXd>& L1,
								  const Eigen::MatrixBase<T2>& L2) {
	// hack to use triangularView: see Eigen documentation
	//"Writing Functions Taking Eigen Types as Parameters"
	// #define _VL const_cast<Eigen::MatrixBase<T1>& >(VL)
	// #define _XL const_cast<Eigen::MatrixBase<T2>& >(XL)
#define _X const_cast<Eigen::MatrixBase<T1>& >(X)
#define _L2 const_cast<Eigen::MatrixBase<T2>& >(L2)
	_X.template triangularView<Eigen::Lower>() = L1 * _L2.template triangularView<Eigen::Lower>();
#undef _X
#undef _L2
	return;
}

template <typename T1>
void Stat::WishartBase::triMultLiX(Ref<MatrixXd> X, const Eigen::MatrixBase<T1>& L) {
#ifdef _MSC_VER
	L.template triangularView<Eigen::Lower>().solveInPlace(X);
#else
	L.template triangularView<Eigen::Lower>().template solveInPlace(X);
#endif
	return;
}


Eigen::MatrixXd Stat::Wishart::sample(int N, Eigen::MatrixXd Psi, Eigen::VectorXd nu, bool inverse ) {
	int q = Psi.rows();
	bool singlePsi = (Psi.cols() == q);
	bool singleNu = (nu.size() == 1);
	// output variables
	MatrixXd X(q,N*q);
	// internal variables
	// TempPQ *tmp = new TempPQ(1,q);
	LLT<MatrixXd> lltq(q);
	MatrixXd Lq = MatrixXd::Zero(q,q);
	MatrixXd Uq = MatrixXd::Zero(q,q);
	MatrixXd Iq = MatrixXd::Identity(q,q);
	MatrixXd PsiL = MatrixXd::Zero(q,q);
	MatrixXd VL = MatrixXd::Zero(q,q);
	//Wishart wish(q);
	wishBase.assignInitValue(q);
	// wishart lower triangular factor
	if(singlePsi) {
		if(!inverse) {
			// tmp->lltq.compute(Psi);
			// PsiL = tmp->lltq.matrixL();
			lltq.compute(Psi);
			PsiL = lltq.matrixL();
		}
		else {
			// ReverseCholesky(PsiL, Psi, tmp->lltq);
			ReverseCholesky(PsiL, Psi, lltq);
		}
	}
	// for-loop
	for(int ii=0; ii<N; ii++) {
		if(!inverse) {
			if(!singlePsi) {
				// tmp->lltq.compute(Psi.block(0,ii*q,q,q));
				// PsiL = tmp->lltq.matrixL();
				lltq.compute(Psi.block(0,ii*q,q,q));
				PsiL = lltq.matrixL();
			}
			// GenerateWishartLowerTri(VL, PsiL, nu(ii*(!singleNu)), tmp->Lq);
			// CrossProdLLt(X.block(0,ii*q,q,q), VL, tmp->Uq);
			wishBase.GenerateLowerTri(VL, PsiL, nu(ii*(!singleNu)));
			CrossProdLLt(X.block(0,ii*q,q,q), VL, Uq);
		}
		else {
			if(!singlePsi) {
				// ReverseCholesky(PsiL, Psi.block(0,ii*q,q,q), tmp->lltq);
				ReverseCholesky(PsiL, Psi.block(0,ii*q,q,q), lltq);
			}
			// GenerateWishartLowerTriXi(VL, PsiL, nu(ii*(!singleNu)));
			wishBase.GenerateLowerTriXi(VL, PsiL, nu(ii*(!singleNu)));
			// InverseLLt(X.block(0,ii*q,q,q), VL, tmp->Lq, tmp->Uq, tmp->Iq);
			// wish.GenerateLowerTriXi(VL, PsiL, nu(ii*(!singleNu)));
			InverseLLt(X.block(0,ii*q,q,q), VL, Lq, Uq, Iq);
		}
	}
	// delete tmp;
	return X;
}

/// @brief Reverse-Cholesky decomposition of a positive-definite matrix
///
/// Calculates the lower triangular matrix `L` satisfying `V = L&apos;L`, where `V` is a positive-definite matrix.
///
/// @param [out] L Lower triangular matrix of size `n x n`.
/// @param [in] V Positive-definite matrix of size `n x n`.
/// @param [in] llt Cholesky solver via reference to `Eigen::LLT` object used for intermediate calculations.
///
/// @note These calculations leave the upper triangular half of `L` unchanged, such that the output is truly triangular only if the upper triangular half of `L` was initialized to zero.
///
inline void Stat::Wishart::ReverseCholesky(Ref<MatrixXd> L,
										   const Ref<const MatrixXd>& V,
										   LLT<MatrixXd>& llt) {
	int q = L.cols();
	int ii, jj;
	// L = anti_transpose(lower_tri(V))
	for(ii=0; ii<q; ii++) {
		for(jj=0; jj<=ii; jj++) {
			L(q-1-jj,q-1-ii) = V(ii,jj);
		}
	}
	// llt "=" chol(L)
	llt.compute(L.triangularView<Eigen::Lower>());
	// L "=" anti_transpose(llq)
	for(ii=0; ii<q; ii++) {
		for(jj=0; jj<=ii; jj++) {
			L(q-1-jj,q-1-ii) = llt.matrixL()(ii,jj);
		}
	}
	// // upper_tri(L) = 0
	// if(upperZero) {
	L = L.triangularView<Eigen::Lower>();
	// }
	return;
}

template <typename T1, typename T2>
void Stat::Wishart::CrossProdLLt(const Eigen::MatrixBase<T1>& X,
								 const Ref<const MatrixXd>& L,
								 const Eigen::MatrixBase<T2>& U) {
	// hack
#define _X const_cast<Eigen::MatrixBase<T1>& >(X)
#define _U const_cast<Eigen::MatrixBase<T2>& >(U)
	_U.template triangularView<Eigen::Upper>() = L.adjoint();
	_X.template triangularView<Eigen::Upper>() = L * U.template triangularView<Eigen::Upper>();
	_X.template triangularView<Eigen::Lower>() = X.adjoint();
#undef _X
#undef _U
	return;
}

template <typename T1>
void Stat::Wishart::InverseLLt(Ref<MatrixXd> X,
							   const Eigen::MatrixBase<T1>& L,
							   Ref<MatrixXd> L2,
							   Ref<MatrixXd> U,
							   const Ref<const MatrixXd>& I) {
	// L2 = L^{-1}
	L2 = L.template triangularView<Eigen::Lower>().solve(I);
	// X = L2' * L2 = L'{-1} * L^{-1} = (L * L')^{-1}
	CrossProdLtL(X, L2, U);
	return;
}

/// @brief Transpose-product of lower triangular matrices
///
/// Performs the multiplication `X = L&apos; * L`, where `L` is a lower triangular matrix.
///
/// @param [out] X Matrix of size `n x n` containing the transpose-product of `L`.
/// @param [in] L Lower triangular matrix of size `n x n`.
/// @param [in] U Upper triangular matrix of size `n x n` used for intermediate calculations.
///
/// @note Does not work properly if any of the inputs arguments are also outputs.
/// @note For safest results, all "off-triangular" elements of triangular arguments should be set to zero.
template <typename T1, typename T2>
void Stat::Wishart::CrossProdLtL(const Eigen::MatrixBase<T1>& X,
								 const Ref<const MatrixXd>& L,
								 const Eigen::MatrixBase<T2>& U) {
	// hack
#define _X const_cast<Eigen::MatrixBase<T1>& >(X)
#define _U const_cast<Eigen::MatrixBase<T2>& >(U)
	_U.template triangularView<Eigen::Upper>() = L.adjoint();
	_X.template triangularView<Eigen::Upper>() = U.template triangularView<Eigen::Upper>() * L;
	_X.template triangularView<Eigen::Lower>() = X.adjoint();
#undef _X
#undef _U
	return;
}

// using namespace arma;
// arma::mat Stat::Wishart::riwish_C(int v, arma::mat S){
//   // Generates a random draw from a inverse-Wishart distribution.
//   // arma::mat CC = chol( inv(S) );
  
//   arma::mat CC = chol( inv(S) );
//   int p = S.n_cols;
//   // Make a diagonal matrix with the elements:
//   // R::rchisq( df ) // with df in a sequence v:(v - p + 1)
//   arma::vec chi_sample(p);
//   // regspace will make a sequence.
//   arma::vec df = arma::regspace(v, v-p+1); // The degrees of freedom when sampling.
//   for( int i=0; i < p; i++ ) {
// 	//  boost::random::chi_squared_distribution<> chiSquared(df);
//     // chi_sample[i] = sqrt( R::rchisq( df[i] ) );
//   }
//   arma::mat Z = diagmat( chi_sample );
//   // Need to fill the upper triangular elements with samples from a standard normal.
//   for( int i=0; i < p-1; i++) {
//     // randn uses a normal distribution to sample.
//     Z(i,arma::span((i+1), (p-1))) = trans(arma::randn(p-(i+1)));
//   }
//   arma::mat out = Z * CC;
//   return inv( trans(out) * out );
// }
