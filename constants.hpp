#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include "Eigen/Core"
#include "Eigen/Sparse"
#include "Eigen/SPQRSupport"
#include "Eigen/QR"

/******************************************************************************/
/***** Define the type used to store the coefficients of the monomials    *****/
/***** and thus also the entries of the matrix whose kernel we want.      *****/
/***** WARNING: SparseQR (at least) seems to always fail with integers!   *****/
/***** I'm not sure why this is; my code never divides coeff_class.       *****/
/******************************************************************************/

//typedef mpq_class coeff_class;				// arbitrary precision rational
//typedef mpfr::mpreal coeff_class;				// arbitrary precision float
typedef double coeff_class;						// ordinary double-width float
//typedef long coeff_class;						// ordinary double-width integer

typedef Eigen::SparseMatrix<coeff_class> Matrix;
typedef Eigen::Matrix<coeff_class, Eigen::Dynamic, Eigen::Dynamic> DMatrix;
typedef Eigen::SparseVector<coeff_class> Vector;
typedef Eigen::Matrix<coeff_class, Eigen::Dynamic, 1> DVector;
typedef Eigen::Triplet<coeff_class> Triplet;

/******************************************************************************/
/***** Define which solver to use to find the kernel of sparse matrices.  *****/
/***** SparseQR and SPQR are both direct (i.e. exact, non-iterative)      *****/
/***** solvers. SPQR requires an external dependency, suitesparse;        *****/
/***** SparseQR is built into Eigen, but appears to be a lot slower.      *****/
/******************************************************************************/

//typedef Eigen::SparseQR<Matrix, Eigen::COLAMDOrdering<int>> QRSolver;
typedef Eigen::SPQR<Matrix> QRSolver; // direct interface to SPQR also possible

/******************************************************************************/
/***** Define which solver to use for dense matrix orthogonalization.     *****/
/******************************************************************************/

//typedef Eigen::ColPivHouseholderQR<DMatrix> DQRSolver; // faster, less info
typedef Eigen::FullPivHouseholderQR<DMatrix> DQRSolver; // slower, full info

/******************************************************************************/
/***** Constants and struct definitions which should be widely accessible *****/
/******************************************************************************/

//constexpr coeff_class delta = 0.25;
constexpr coeff_class EPSILON = 1e-10;
constexpr unsigned int MAX_THREADS = 8u;

struct arguments{
	int numP;
	int degree;
	coeff_class delta;
	int options;
};

struct particle{
	int pm;	// P_-
	int pt; // P_\perp
	int pp; // P_+

	bool operator==(const particle& other) const;
	bool operator!=(const particle& other) const { return !(*this == other); }
};

enum options { OPT_BRUTE = 1 << 0, OPT_VERSION = 1 << 1, OPT_DEBUG = 1 << 2,
				OPT_PARITYONLY = 1 << 3, OPT_EQNMOTION = 1 << 4,
				OPT_MSORTING = 1 << 5, OPT_DIRICHLET = 1 << 6,
				OPT_OUTPUT = 1 << 7, OPT_IPTEST = 1 << 8, OPT_ALLMINUS = 1 << 9};

/******************************************************************************/
/***** Compile-time constant math functions                               *****/
/******************************************************************************/

constexpr coeff_class Pochhammer(const coeff_class A, const int n){
	coeff_class ret = 1;
	for(int i = 0; i < n; ++i){
		ret *= A + i;
	}
	return ret;
}

constexpr coeff_class Factorial(const int n){
	return Pochhammer(1, n);
}

// this would be constexpr instead of inline if the cmath functions were 
// properly tagged. It likely compiles as constexpr in GCC.
inline coeff_class Binomial(const int n, const int k){
	coeff_class logRet = 0;
	logRet += std::lgamma(n+1);
	logRet -= std::lgamma(k+1);
	logRet -= std::lgamma(n-k+1);
	return std::exp(logRet);

	//return Factorial(n)/(Factorial(k) * Factorial(n-k));
}

#endif
