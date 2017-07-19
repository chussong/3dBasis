#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include "Eigen/Core"
#include "Eigen/Sparse"
#include "Eigen/SPQRSupport"

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

//constexpr coeff_class delta = 0.25;
constexpr coeff_class EPSILON = 1e-5;

typedef Eigen::SparseMatrix<coeff_class> Matrix;
typedef Eigen::Matrix<coeff_class, Eigen::Dynamic, Eigen::Dynamic> DMatrix;
typedef Eigen::SparseVector<coeff_class> Vector;
typedef Eigen::Matrix<coeff_class, Eigen::Dynamic, 1> DVector;
typedef Eigen::Triplet<coeff_class> Triplet;

/******************************************************************************/
/***** Define which solver to use to find the kernel of the matrix.       *****/
/***** SparseQR and SPQR are both direct (i.e. exact, non-iterative)      *****/
/***** solvers. SPQR requires an external dependency, suitesparse;        *****/
/***** SparseQR is built into Eigen, but appears to be a lot slower.      *****/
/******************************************************************************/
//typedef Eigen::SparseQR<Matrix, Eigen::COLAMDOrdering<int>> QRSolver;
typedef Eigen::SPQR<Matrix> QRSolver; // direct interface to SPQR also possible

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
				OPT_MSORTING = 1 << 5, OPT_DIRICHLET = 1 << 6 };

#endif
