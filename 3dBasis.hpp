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
#include <thread>

constexpr char VERSION[] = "0.7.2";
constexpr char RELEASE_DATE[] = __DATE__;

#include "constants.hpp"
#include "construction.hpp"
#include "mono.hpp"
#include "poly.hpp"
#include "basis.hpp"
#include "io.hpp"
#include "timer.hpp"
#include "cache.hpp"
#include "matrix.hpp"
//#include "multinomial.hpp" // only needed for testing, really

std::ostream& operator<<(std::ostream& os, const Triplet& out);

// startup and input parsing --------------------------------------------------

arguments ParseArguments(int argc, char* argv[]);
int ParseOptions(std::vector<std::string> options);

int FindPrimaries(const arguments& args);
int FindPrimariesParityOnly(const arguments& args);
int FindPrimariesBruteForce(const arguments& args);
int FindPrimariesByM(const arguments& args);
int InnerProductTest(const arguments& args);

int Orthogonalize(const std::vector<Basis<Mono>>& inputBases,
					const GammaCache& cache, const KVectorCache& kCache);
std::vector<Poly> GramSchmidt(const std::vector<Basis<Mono>> input,
		const GammaCache& cache, const KVectorCache& kCache);
std::vector<Poly> GramSchmidt(const Basis<Mono> input, const GammaCache& cache, 
		const KVectorCache& kCache);
std::vector<Poly> GramSchmidt_WithMatrix(const std::vector<Basis<Mono>> input,
		const DMatrix& gramMatrix);
std::vector<Poly> GramSchmidt_WithMatrix(const Basis<Mono> input, 
		const DMatrix& gramMatrix);
DVector GSProjection(const DVector& toProject, const DVector& projectOnto,
		const DMatrix& gramMatrix);
coeff_class GSNorm(const DVector& vector, const DMatrix& gramMatrix);
std::vector<Poly> GramSchmidt_MatrixOnly(const DMatrix& input, 
		const std::vector<Basis<Mono>>& inputBases);

// functions interfacing with Eigen ------------------------------------------

std::list<Triplet> ConvertToRows(const std::vector<Poly>& PolyForms, 
		const Basis<Mono>& targetBasis, const Eigen::Index rowOffset);
std::vector<Poly> CombineKernels(const std::vector<Poly>& kernel1,
		const std::vector<Poly>& kernel2);
Poly VectorToPoly(const Vector& kernelVector, const Basis<Mono>& startBasis);
Poly VectorToPoly(const DVector& kernelVector, const Basis<Mono>& startBasis);
Poly ColumnToPoly(const Matrix& kernelMatrix, const Eigen::Index col, 
		const Basis<Mono>& startBasis);
Poly ColumnToPoly(const DMatrix& kernelMatrix, const Eigen::Index col, 
		const Basis<Mono>& startBasis);

DMatrix ExtractQMatrix(const Eigen::FullPivHouseholderQR<DMatrix>& solver, 
		               const int dimension);
DMatrix ExtractQMatrix(const Eigen::ColPivHouseholderQR<DMatrix>& solver, 
		               const int dimension);
void ClearZeros(DMatrix* toClear);
//DMatrix GramMatrix(const Basis<Mono>& basis);

// miscellaneous -------------------------------------------------------------

Basis<Mono> MinimalBasis(const std::vector<Poly>& polynomials);

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
inline Poly VectorToPoly(const DVector& kernelVector, const Basis<T>& startBasis){
	Poly ret;
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
inline std::vector<Poly> Kernel(const Matrix& KActions, const Basis<T>& startBasis,
		const int options, const bool outputKernel){
	if(KActions.rows() == 0 || KActions.cols() == 0) return std::vector<Poly>();
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
}

// both GramMatrix functions construct a dense matrix containing the inner 
// product of every basis element with every other basis element. They split the
// inner products up among (up to) MAX_THREADS different threads, as defined in
// constants.hpp, with each thread doing approximately the same number of 
// elements.
template<class T>
void ThreadFillGram(const Basis<T>& basis, const GammaCache& cache,
		const KVectorCache& kCache, std::vector<Triplet>* output, 
		const size_t startRow, const size_t endRow) {
	for(size_t row = startRow; row <= endRow; ++row){
		output->at(basis.size()*row + row) = Triplet(row, row,
				T::InnerProduct(basis[row], basis[row], cache, kCache) );
		for(size_t col = row + 1; col < basis.size(); ++col){
			coeff_class prod = T::InnerProduct(basis[row], basis[col], cache, 
					kCache);
			output->at(basis.size()*row + col) = Triplet(row, col, prod);
			output->at(basis.size()*col + row) = Triplet(col, row, prod);
			// we could just be filling these in order instead of jumping around
			// if the final order really matters you can sort them at the end
		}
	}
}

template<class T>
DMatrix GramMatrix(const Basis<T>& basis, const GammaCache& cache,
		const KVectorCache& kCache){
	if(basis.size() == 0) return DMatrix();
	unsigned int numThreads = std::min(static_cast<unsigned long>(MAX_THREADS), 
			basis.size()/2);
	std::vector<std::thread> threads;
	std::vector<Triplet> entries(basis.size()*basis.size());
	size_t entryCount = 0;
	size_t totalCount = (basis.size() * (basis.size() + 1) )/2;
	size_t startRow = 0;

	numThreads = 1;

	for(size_t row = 0; row < basis.size(); ++row){
		entryCount += basis.size() - row;
		if(entryCount > totalCount/numThreads/* || row == basis.size()-1*/){
			threads.emplace_back(ThreadFillGram<T>, std::cref(basis), 
					std::cref(cache), std::cref(kCache), &entries, startRow, row);
			entryCount = 0;
			startRow = row + 1;
		}
	}
	// the main thread fills in the last batch on its own
	ThreadFillGram<T>(basis, cache, kCache, &entries, startRow, 
			basis.size()-1);
	for(auto& thread : threads) thread.join();
	//for(auto& entry : entries) std::cout << entry << std::endl;
	Matrix gram(basis.size(), basis.size());
	gram.setFromTriplets(entries.begin(), entries.end());
	return gram;
}

template<class T>
const T& Get(const std::vector<Basis<T>>& multipleBases, size_t index){
	for(const auto& basis : multipleBases){
		if(index < basis.size()){
			return basis[index];
		} else {
			index -= basis.size();
		}
	}
	throw std::runtime_error("out of range error in Get(vector<Basis>)");
}

template<class T>
void ThreadFillGramFromVec(const std::vector<Basis<T>>& bases, 
		const GammaCache& cache, const KVectorCache& kCache, 
		std::vector<Triplet>* output, const size_t startRow, 
		const size_t endRow) {
	size_t totalSize = 0;
	for(const auto& basis : bases) totalSize += basis.size();
	for(size_t row = startRow; row <= endRow; ++row){
		for(size_t col = row; col < totalSize; ++col){
			output->at(totalSize*row + col) = Triplet(row, col, 
					T::InnerProduct(Get(bases, row), Get(bases, col), 
						cache, kCache) );
			if(row != col){
				output->at(totalSize*col + row) = Triplet(col, row, 
						output->at(totalSize*row + col).value());
			}
		}
	}
}

// note: this function spawns a bunch of threads, but it should do part of the
// computation on its own as well instead of just waiting around.
template<class T>
DMatrix GramMatrix(const std::vector<Basis<T>>& allBases,
					const GammaCache& cache, const KVectorCache& kCache){
	// could also do this by defining some kind of super iterator that traverses 
	// a vector of bases?
	size_t totalSize = 0;
	for(auto& basis : allBases) totalSize += basis.size();
	if(totalSize == 0) return DMatrix();
	unsigned int numThreads = std::min(static_cast<unsigned long>(MAX_THREADS), 
			totalSize/2);
	std::vector<std::thread> threads;
	std::vector<Triplet> entries(totalSize * totalSize);
	size_t entryCount = 0;
	size_t totalCount = (totalSize * (totalSize + 1) )/2;
	size_t startRow = 0;

	numThreads = 1;

	for(size_t row = 0; row < totalSize; ++row){
		entryCount += totalSize - row;
		if(entryCount > totalCount/numThreads/* || row == totalSize-1*/){
			threads.emplace_back(ThreadFillGramFromVec<T>, std::cref(allBases), 
					std::cref(cache), std::cref(kCache), &entries, startRow, row);
			entryCount = 0;
			startRow = row + 1;
		}
	}
	// the main thread fills in the last batch on its own
	ThreadFillGramFromVec<T>(allBases, cache, kCache, &entries, startRow, 
			totalSize-1);
	for(auto& thread : threads) thread.join();
	//for(auto& entry : entries) std::cout << entry << std::endl;
	Matrix gram(totalSize, totalSize);
	gram.setFromTriplets(entries.begin(), entries.end());
	return gram;
}

template<class T>
inline size_t TotalSize(const std::vector<T>& vectorOfContainers){
	size_t totalSize = 0;
	for(const T& element : vectorOfContainers){
		totalSize += element.size();
	}
	return totalSize;
}

#endif
