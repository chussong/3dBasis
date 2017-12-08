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
#include <type_traits>	// std::is_same

constexpr char VERSION[] = "0.8.2";
constexpr char RELEASE_DATE[] = __DATE__;

#include "constants.hpp"
#include "construction.hpp"
#include "mono.hpp"
#include "poly.hpp"
#include "basis.hpp"
#include "io.hpp"
#include "timer.hpp"
#include "gram-schmidt.hpp"
#include "matrix.hpp"
#include "multinomial.hpp" // the binomial coefficients are initialized here

// startup and input parsing --------------------------------------------------

arguments ParseArguments(int argc, char* argv[]);
int ParseOptions(std::vector<std::string> options);

// int FindPrimaries(const arguments& args);
// int FindPrimariesParityOnly(const arguments& args);
// int FindPrimariesBruteForce(const arguments& args);
// int FindPrimariesByM(const arguments& args);
int InnerProductTest(const arguments& args);

// functions interfacing with Eigen ------------------------------------------

// std::list<Triplet> ConvertToRows(const std::vector<Poly>& PolyForms, 
		// const Basis<Mono>& targetBasis, const Eigen::Index rowOffset);
std::vector<Poly> CombineKernels(const std::vector<Poly>& kernel1,
		const std::vector<Poly>& kernel2);
// Poly VectorToPoly(const SVector& kernelVector, const Basis<Mono>& startBasis);
Poly VectorToPoly(const DVector& kernelVector, const Basis<Mono>& startBasis);
// Poly ColumnToPoly(const SMatrix& kernelMatrix, const Eigen::Index col, 
		// const Basis<Mono>& startBasis);
Poly ColumnToPoly(const DMatrix& kernelMatrix, const Eigen::Index col, 
		const Basis<Mono>& startBasis);

void ClearZeros(DMatrix* toClear);
//DMatrix GramMatrix(const Basis<Mono>& basis);

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

/*template<class T>
inline std::vector<Poly> Kernel(const SMatrix& KActions, const Basis<T>& startBasis,
		const int options, const bool outputKernel){
	if(KActions.rows() == 0 || KActions.cols() == 0) return std::vector<Poly>();
	if(options & OPT_DEBUG){
		std::cout << "Computing kernel from below matrix..." << std::endl;
		std::cout << KActions << std::endl;
	}
	SQRSolver solver;
	solver.setPivotThreshold(EPSILON); // norms smaller than this are zero
	solver.compute(KActions.transpose());
	//std::cout << "Solved. Found rank " << solver.rank() << ", i.e. "
		//<< startBasis.size() - solver.rank() << " kernel elements." << std::endl;

	//std::cout << "Converting the kernel to polynomials..." << std::endl;

	DVector projector = Eigen::VectorXd::Zero(startBasis.size());
	DVector kernelVector(startBasis.size());
	std::vector<Poly> ret;
	ret.resize(startBasis.size() - solver.rank());

	if(!outputKernel) return ret;

	for(auto col = 0u; col < startBasis.size() - solver.rank(); ++col){
		projector(solver.rank() + col-1) = 0;
		projector(solver.rank() + col) = 1;
		kernelVector = solver.matrixQ()*projector;
		if(options & OPT_DEBUG){
			std::cout << "Projecting out with this: " << projector << std::endl;
			std::cout << kernelVector << "\n----------" << std::endl;
		}
		ret[col] = VectorToPoly(kernelVector, startBasis);
	}

	return ret;
}*/

#endif
