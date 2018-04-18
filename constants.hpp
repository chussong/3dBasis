#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <iostream> // for the output stream in the argument struct
#ifndef NO_GUI
#include <QtCore/QTextStream>
#include <sstream>
#endif

#include "Eigen/Core"
// #include "Eigen/Sparse"
// #include "Eigen/SPQRSupport"
#include "Eigen/QR"
#include "Eigen/Eigenvalues"

/******************************************************************************/
/***** Define the type used to store the coefficients of the monomials    *****/
/***** and thus also the entries of the matrix whose kernel we want.      *****/
/******************************************************************************/

//typedef mpq_class coeff_class;				// arbitrary precision rational
//typedef mpfr::mpreal coeff_class;				// arbitrary precision float
//typedef double coeff_class;					// ordinary double-width float
//typedef long coeff_class;						// ordinary double-width integer
//typedef long double coeff_class;
typedef __float128 coeff_class;

typedef Eigen::Matrix<coeff_class, Eigen::Dynamic, Eigen::Dynamic> DMatrix;
typedef Eigen::Matrix<coeff_class, Eigen::Dynamic, 1> DVector;
// typedef Eigen::SparseMatrix<coeff_class> SMatrix;
// typedef Eigen::SparseVector<coeff_class> SVector;
// typedef Eigen::Triplet<coeff_class> Triplet;

// also define a built-in class to which coeff_class is implicitly convertible
// to help with template resolution in stdlib functions. If coeff_class is
// built-in, then this can just be the same thing; otherwise, it should be
// something with specializations for std::sqrt, etc. You could also just define
// new specializations for these functions.
typedef double builtin_class;

// if builtin_class and coeff_class are different, this ostream operator also
// needs to be defined

#ifdef NO_GUI
typedef std::ostream OStream;
using endl = std::endl;
#else
typedef QTextStream OStream;
#endif

inline std::ostream& operator<<(std::ostream& os, const coeff_class& out){
    // std::cout << "sending out a coeff_class" << std::endl;
    return os << static_cast<builtin_class>(out);
}

inline std::ostream& operator<<(std::ostream& os, const DMatrix& out) {
    // std::cout << "sending out a matrix" << std::endl;
    return os << out.cast<builtin_class>();
}

inline std::ostream& operator<<(std::ostream& os, const DVector& out) {
    // std::cout << "sending out a vector" << std::endl;
    return os << out.cast<builtin_class>();
}

/******************************************************************************/
/***** Define which solver to use to find the kernel of sparse matrices.  *****/
/***** SparseQR and SPQR are both direct (i.e. exact, non-iterative)      *****/
/***** solvers. SPQR requires an external dependency, suitesparse;        *****/
/***** SparseQR is built into Eigen, but appears to be a lot slower.      *****/
/******************************************************************************/

// typedef Eigen::SparseQR<SMatrix, Eigen::COLAMDOrdering<int>> QRSolver;
// typedef Eigen::SPQR<SMatrix> SQRSolver; // direct interface to SPQR also possible

/******************************************************************************/
/***** Define which solver to use for dense matrix orthogonalization.     *****/
/******************************************************************************/

//typedef Eigen::ColPivHouseholderQR<DMatrix> DQRSolver; // faster, less info
// typedef Eigen::FullPivHouseholderQR<DMatrix> DQRSolver; // slower, full info
typedef Eigen::FullPivHouseholderQR<
		Eigen::Matrix<builtin_class, Eigen::Dynamic, Eigen::Dynamic> > DQRSolver;

/******************************************************************************/
/***** Define which solver to use for dense eigenvalue decomposition.     *****/
/******************************************************************************/

// the below only works for symmetric matrices
typedef Eigen::SelfAdjointEigenSolver<
		Eigen::Matrix<builtin_class, Eigen::Dynamic, Eigen::Dynamic> > EigenSolver;

/******************************************************************************/
/***** Constants and struct definitions which should be widely accessible *****/
/******************************************************************************/

//constexpr coeff_class delta = 0.25;
constexpr coeff_class EPSILON = 1e-8;
constexpr unsigned int MAX_THREADS = 8u;

struct Arguments {
    int numP = -1;
    int degree = -1;
    std::size_t partitions = 4; // \mu partitions per operator pair
    coeff_class msq = 1; // the coefficient of the mass term
    int options = 0;
    OStream* outStream = nullptr;
    OStream* console = nullptr;
};

struct particle {
    char pm;    // P_-
    char pt;    // P_\perp

    bool operator==(const particle& other) const;
    bool operator!=(const particle& other) const { return !(*this == other); }
};

inline bool particle::operator==(const particle& other) const {
	return (pm == other.pm) && (pt == other.pt);
}

enum options { OPT_BRUTE = 1 << 0, OPT_VERSION = 1 << 1, OPT_DEBUG = 1 << 2,
                OPT_PARITYONLY = 1 << 3, OPT_EQNMOTION = 1 << 4,
                OPT_MSORTING = 1 << 5, OPT_DIRICHLET = 1 << 6,
                OPT_OUTPUT = 1 << 7, OPT_IPTEST = 1 << 8, OPT_ALLMINUS = 1 << 9,
                OPT_MULTINOMTEST = 1 << 10, OPT_TEST = 1 << 11,
                OPT_STATESONLY = 1 << 12, OPT_MATHEMATICA = 1 << 13 };

enum MATRIX_TYPE { MAT_KINETIC, MAT_INNER, MAT_MASS, MAT_INTER_SAME_N, 
    MAT_INTER_N_PLUS_2 };

/******************************************************************************/
/***** Compile-time constant math functions                               *****/
/******************************************************************************/

constexpr coeff_class Pochhammer(const coeff_class A, const int n) {
    coeff_class ret = 1;
    for(int i = 0; i < n; ++i){
        ret *= A + i;
    }
    return ret;
}

constexpr coeff_class Factorial(const int n) {
    return Pochhammer(1, n);
}

// this would be constexpr instead of inline if the cmath functions were 
// properly tagged. It likely compiles as constexpr in GCC.
/*inline coeff_class Binomial(const int n, const int k){
	coeff_class logRet = 0;
	logRet += std::lgamma(n+1);
	logRet -= std::lgamma(k+1);
	logRet -= std::lgamma(n-k+1);
	return std::exp(logRet);

	//return Factorial(n)/(Factorial(k) * Factorial(n-k));
}*/

// return total size of a vector of containers
template<class T>
inline std::size_t TotalSize(const std::vector<T>& vectorOfContainers) {
    std::size_t totalSize = 0;
    for(const T& element : vectorOfContainers){
        totalSize += element.size();
    }
    return totalSize;
}

// return the vector sum of A and B
template<typename T>
std::vector<T> AddVectors(const std::vector<T>& A, const std::vector<T>& B) {
    const std::vector<T>* aP = &A;
    const std::vector<T>* bP = &B;
    if (A.size() < B.size()) std::swap(aP, bP);
    std::vector<T> output(*aP);
    for (typename std::vector<T>::size_type i = 0; i < bP->size(); ++i) {
        output[i] += (*bP)[i];
    }
    return output;
}

#endif
