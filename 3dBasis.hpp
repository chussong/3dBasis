#ifndef MAIN_HPP
#define MAIN_HPP

#include <iostream>
#include <exception>
#include <string>
#include <vector>
#include <array>
#include <list>
#include <utility>		// std::pair
#include <algorithm>	// std::remove_if
#include "Eigen/Core"
#include "Eigen/Sparse"
#include "Eigen/SPQRSupport"

constexpr char VERSION[] = "0.4.2";
constexpr char RELEASE_DATE[] = __DATE__;

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

#include "mono.hpp"
#include "poly.hpp"
#include "basis.hpp"

bool ParticlePrecedence(const particle& a, const particle& b);

std::ostream& operator<<(std::ostream& os, const Triplet& out);

// startup and input parsing --------------------------------------------------

enum options { OPT_BRUTE = 1 << 0, OPT_VERSION = 1 << 1, OPT_DEBUG = 1 << 2,
				OPT_PARITYONLY = 1 << 3, OPT_EQNMOTION = 1 << 4,
				OPT_MSORTING = 1 << 5 };
arguments ParseArguments(int argc, char* argv[]);
int ParseOptions(std::vector<std::string> options);

int FindPrimaries(const arguments& args);
int FindPrimariesParityOnly(const arguments& args);
int FindPrimariesBruteForce(const arguments& args);

// functions interfacing with Eigen ------------------------------------------

Matrix KMatrix(const basis& startingBasis, const basis& targetBasis,
		const coeff_class delta);
Matrix K13Matrix(const basis& startingBasis, const basis& targetBasis,
		const coeff_class delta);
Matrix K2Matrix(const basis& startingBasis, const basis& targetBasis,
		const coeff_class delta);
std::array<Matrix,4> KMatrices(const splitBasis& startingBasis,
		const splitBasis& targetBasis, const coeff_class delta);
std::list<Triplet> ConvertToRows(const std::vector<poly>& polyForms, 
		const basis& targetBasis, const Eigen::Index rowOffset);
std::vector<poly> Kernel(const Matrix& KActions, const basis& startBasis,
		const bool outputKernel);
std::vector<poly> CombineKernels(const std::vector<poly>& kernel1,
		const std::vector<poly>& kernel2);
poly VectorToPoly(const Vector& kernelVector, const basis& startBasis);
poly VectorToPoly(const DVector& kernelVector, const basis& startBasis);
poly ColumnToPoly(const Matrix& kernelMatrix, const Eigen::Index col, 
		const basis& startBasis);
poly ColumnToPoly(const DMatrix& kernelMatrix, const Eigen::Index col, 
		const basis& startBasis);


// templates -----------------------------------------------------------------

template<typename ParseTo>
ParseTo ReadArg(const std::string& arg){
	std::cerr << "Error: attempted to parse the argument " << arg << " to a "
		<< "type with no known parsing function. Please specialize the ReadArg "
		<< "template to your coeff_class." << std::endl;
	return ParseTo();
}

template<>
int ReadArg<int>(const std::string& arg){
	int ret;
	try{ret = std::stoi(arg);}
	catch(const std::invalid_argument &e){
		std::cerr << "Error: this non-option argument could not be "
			<< "converted to an integer: " << arg << std::endl;
		throw;
	}
	catch(const std::out_of_range &e){
		std::cerr << "Error: specification of N or degree is too "
			<< "large to store. This computation would never finish"
			<< " anyway..." << std::endl;
		throw;
	}
	return ret;
}

template<>
double ReadArg<double>(const std::string& arg){
	double ret;
	try{ret = std::stod(arg);}
	catch(const std::invalid_argument &e){
		std::cerr << "Error: this non-option argument could not be "
			<< "converted to an integer: " << arg << std::endl;
		throw;
	}
	catch(const std::out_of_range &e){
		std::cerr << "Error: specification of Delta is too large to store."
			<< std::endl;
		throw;
	}
	return ret;
}

// T can be any type or class with an == operator; value indexes T by uint
template<typename Accessor>
inline std::vector<int> IdentifyNodes(Accessor A, const size_t size){
	std::vector<int> nodes;
	int newNode;
	newNode = 1;
	for(unsigned int i = 1; i < size; ++i){
		if(A(i) != A(i-1)){
			nodes.push_back(newNode);
			newNode = 1;
		} else {
			++newNode;
		}
	}
	nodes.push_back(newNode);
	return nodes;
}

// calls the generic IdentifyNodes using the class's particles and (*value)
template<typename T>
inline std::vector<int> mono::IdentifyNodes(T (*value)(particle)) const{
	return IdentifyNodes([this, value](unsigned int i){return value(this->particles[i]);},
			NParticles());
}

// should work for any containerish thing T with operator[]
template<typename T>
inline std::vector<int> IdentifyNodes(const T& container){
	return IdentifyNodes([container](unsigned int i){return container[i];},
			container.size());
}

// default behavior for a vector of particles; used by mono::IdentifyNodes()
template<>
inline std::vector<int> IdentifyNodes(const std::vector<particle>& particles){
	return IdentifyNodes([particles](unsigned int i){return 
			std::array<int, 3>({{particles[i].pm, particles[i].pt, particles[i].pp}}); },
			particles.size());
}

// stream output operator template for vectors
template<typename T>
inline std::ostream& operator<<(std::ostream& os, const std::vector<T>& out){
	os << "{";
	for(auto& element : out){
		if(element >= 0) os << " ";
		os << element << ",";
	}
	os << "\b }";
	return os;
}

// specialization of above template for vectors of polynomials
template<>
inline std::ostream& operator<<(std::ostream& os, const std::vector<poly>& out){
	if(out.size() == 0) return os << "{ }";
	os << "{ ";
	for(auto& element : out){
		os << element.HumanReadable() << " | ";
	}
	os << "\b\b \b}";
	return os;
}

// specialization of above template which "transposes" particle vectors
template<>
inline std::ostream& operator<<(std::ostream& os, const std::vector<particle>& out){
	os << "{";
	for(auto& p : out){
		if(p.pm >= 0) os << " ";
		os << p.pm << ",";
	}
	os << "\b }{";
	for(auto& p : out){
		if(p.pt >= 0) os << " ";
		os << p.pt << ",";
	}
	os << "\b }{";
	for(auto& p : out){
		if(p.pp >= 0) os << " ";
		os << p.pp << ",";
	}
	os << "\b }";
	return os;
}

#endif
