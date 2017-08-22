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

constexpr char VERSION[] = "0.5.3";
constexpr char RELEASE_DATE[] = __DATE__;

#include "constants.hpp"
#include "construction.hpp"
#include "mono.hpp"
#include "poly.hpp"
#include "basis.hpp"
#include "io.hpp"
#include "timer.hpp"

std::ostream& operator<<(std::ostream& os, const Triplet& out);

// startup and input parsing --------------------------------------------------

arguments ParseArguments(int argc, char* argv[]);
int ParseOptions(std::vector<std::string> options);

int FindPrimaries(const arguments& args);
int FindPrimariesParityOnly(const arguments& args);
int FindPrimariesBruteForce(const arguments& args);
int FindPrimariesByM(const arguments& args);
unsigned int AddPrimariesAtL(const mBasis& startBasis, const mBasis& targetBasis,
		const unsigned int L, std::vector<poly>& primaries, 
		const coeff_class delta, const int options);
int InnerProductTest(const arguments& args);

int Orthogonalize(const std::vector<Basis<mono>>& inputBases);

// functions interfacing with Eigen ------------------------------------------

Matrix KMatrix(const Basis<mono>& startingBasis, const Basis<mono>& targetBasis,
		const coeff_class delta, const int options);
Matrix K13Matrix(const Basis<mono>& startingBasis, const Basis<mono>& targetBasis,
		const coeff_class delta);
Matrix K2Matrix(const Basis<mono>& startingBasis, const Basis<mono>& targetBasis,
		const coeff_class delta);
std::array<Matrix,4> KMatrices(const splitBasis<mono>& startingBasis,
		const splitBasis<mono>& targetBasis, const coeff_class delta);
std::list<Triplet> ConvertToRows(const std::vector<poly>& polyForms, 
		const Basis<mono>& targetBasis, const Eigen::Index rowOffset);
/*std::vector<poly> Kernel(const Matrix& KActions, const Basis<mono>& startBasis,
		const bool outputKernel);*/
std::vector<poly> CombineKernels(const std::vector<poly>& kernel1,
		const std::vector<poly>& kernel2);
poly VectorToPoly(const Vector& kernelVector, const Basis<mono>& startBasis);
poly VectorToPoly(const DVector& kernelVector, const Basis<mono>& startBasis);
poly ColumnToPoly(const Matrix& kernelMatrix, const Eigen::Index col, 
		const Basis<mono>& startBasis);
poly ColumnToPoly(const DMatrix& kernelMatrix, const Eigen::Index col, 
		const Basis<mono>& startBasis);

DMatrix ExtractQMatrix(const Eigen::FullPivHouseholderQR<DMatrix>& solver, 
		               const int dimension);
DMatrix ExtractQMatrix(const Eigen::ColPivHouseholderQR<DMatrix>& solver, 
		               const int dimension);
//DMatrix GramMatrix(const Basis<mono>& basis);

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

template<class T>
inline poly VectorToPoly(const DVector& kernelVector, const Basis<T>& startBasis){
	poly ret;
	if(static_cast<size_t>(kernelVector.rows()) != startBasis.size()){
		std::cerr << "Error: the given Q column has " << kernelVector.rows()
			<< " rows, " << "but the given basis has " << startBasis.size() 
			<< " monomials. These must be the same." << std::endl;
		return ret;
	}
	for(auto row = 0; row < kernelVector.rows(); ++row){
		if(std::abs(kernelVector.coeff(row)) < EPSILON) continue;
		ret += kernelVector.coeff(row)*startBasis[row];
	}

	if(ret.size() == 0) return ret;
	coeff_class smallestCoeff = std::abs(ret[0].Coeff());
	for(auto& term : ret) smallestCoeff = std::min(std::abs(term.Coeff()), smallestCoeff);
	for(auto& term : ret) term /= smallestCoeff;
	return ret;
}

template<class T>
inline std::vector<poly> Kernel(const Matrix& KActions, const Basis<T>& startBasis,
		const int options, const bool outputKernel){
	if(KActions.rows() == 0 || KActions.cols() == 0) return std::vector<poly>();
	if(options & OPT_DEBUG){
		std::cout << "Computing kernel from below matrix..." << std::endl;
		std::cout << KActions << std::endl;
	}
	QRSolver solver;
	solver.setPivotThreshold(EPSILON); // norms smaller than this are zero
	solver.compute(KActions.transpose());
	//std::cout << "Solved. Found rank " << solver.rank() << ", i.e. "
		//<< startBasis.size() - solver.rank() << " kernel elements." << std::endl;

	//std::cout << "Converting the kernel to polynomials..." << std::endl;

	DVector projector = Eigen::VectorXd::Zero(startBasis.size());
	DVector kernelVector(startBasis.size());
	std::vector<poly> ret;
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
}

template<class T>
DMatrix GramMatrix(const Basis<T>& basis){
	std::vector<Triplet> entries;
	for(auto row = 0u; row < basis.size(); ++row){
		for(auto col = row; col < basis.size(); ++col){
			entries.emplace_back(row, col, 
					T::InnerProduct(basis[row], basis[col]));
			if(row != col){
				entries.emplace_back(col, row, 
						T::InnerProduct(basis[row], basis[col]));
			}
		}
	}
	Matrix gram(basis.size(), basis.size());
	gram.setFromTriplets(entries.begin(), entries.end());
	return gram;
}

template<class T>
DMatrix GramMatrix(const std::vector<Basis<T>>& allBases){
	// could also do this by defining some kind of super iterator that traverses 
	// a vector of bases?
	std::vector<Triplet> entries;
	auto row = 0u;
	size_t totalSize = 0;
	for(auto& basisA : allBases){
		for(auto& elementA : basisA){
			auto col = 0u;
			for(auto& basisB : allBases){
				for(auto& elementB : basisB){
					coeff_class product = T::InnerProduct(elementA, elementB);
					entries.emplace_back(row, col, product);
					if(row != col) entries.emplace_back(col, row, product);
					++col;
					if(col > row) break;
				}
				if(col > row) break;
			}
			++row;
		}
		totalSize += basisA.size();
	}
	Matrix gram(totalSize, totalSize);
	gram.setFromTriplets(entries.begin(), entries.end());
	return gram;
}

/*template<typename T> struct ValidQSolver { static constexpr bool value = false; };
template<> struct ValidQSolver<Eigen::FullPivHouseholderQR<DMatrix>> { 
	static constexpr bool value = true;
};
template<> struct ValidQSolver<Eigen::ColPivHouseholderQR<DMatrix>> { 
	static constexpr bool value = true;
};

template<class T>
inline DMatrix ExtractQMatrix(const T&, const int){
	static_assert(ValidQSolver<T>::value, "Solver must be set to one of the "
			"properly handled types enumerated in 3dBasis.hpp.");
	return DMatrix();
}

template<>
inline DMatrix ExtractQMatrix(const Eigen::FullPivHouseholderQR<DMatrix>& solver, 
		               const int){
	return solver.matrixQ();
}

template<>
inline DMatrix ExtractQMatrix(const Eigen::ColPivHouseholderQR<DMatrix>& solver, 
		               const int dimension){
	return solver.householderQ()*DMatrix::Identity(dimension, dimension);
}*/

template<class T>
inline size_t TotalSize(const std::vector<T>& vectorOfContainers){
	size_t totalSize = 0;
	for(const T& element : vectorOfContainers){
		totalSize += element.size();
	}
	return totalSize;
}

#endif
