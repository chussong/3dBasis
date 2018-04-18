#ifndef GRAM_SCHMIDT_HPP
#define GRAM_SCHMIDT_HPP

// functions related to orthogonalization go here, including potentially non-GS 
// methods which would be more stable

#include <vector>
#include <iostream>

#include "constants.hpp"
#include "timer.hpp"
#include "mono.hpp"
#include "poly.hpp"
#include "basis.hpp"
#include "matrix.hpp"

// this should be the only function called from outside of this file ----------

std::vector<Poly> Orthogonalize(const std::vector<Basis<Mono>>& inputBases, 
                OStream& console, const bool odd);

// custom gram-schmidt --------------------------------------------------------

std::vector<Poly> GramSchmidt(const std::vector<Basis<Mono>> input);
std::vector<Poly> GramSchmidt(const Basis<Mono> input);
std::vector<Poly> GramSchmidt_WithMatrix(const std::vector<Basis<Mono>> input,
		const DMatrix& gramMatrix);
std::vector<Poly> GramSchmidt_WithMatrix(const Basis<Mono> input, 
		const DMatrix& gramMatrix);
std::vector<Poly> GramSchmidt_WithMatrix_A(const Basis<Mono> inputBasis, 
		const DMatrix& gramMatrix);
std::vector<Poly> GramSchmidt_WithMatrix_B(const Basis<Mono> inputBasis, 
		const DMatrix& gramMatrix);
DVector GSProjection(const DVector& toProject, const DVector& projectOnto,
		const DMatrix& gramMatrix);
coeff_class GSNorm(const DVector& vector, const DMatrix& gramMatrix);
std::vector<Poly> GramSchmidt_MatrixOnly(const DMatrix& input, 
		const Basis<Mono>& inputBases);

// interface with matrix QR decompositions (as alternative to GS) ------------

std::vector<Poly> PolysFromQMatrix(const DMatrix& QMatrix, 
		const Basis<Mono>& basis, const DMatrix& gramMatrix, 
		const Eigen::Index rank);

// miscellaneous -------------------------------------------------------------

DMatrix ExtractQMatrix(const Eigen::FullPivHouseholderQR<DMatrix>& solver, 
		               const int dimension);
DMatrix ExtractQMatrix(const Eigen::ColPivHouseholderQR<DMatrix>& solver, 
		               const int dimension);

// deprecated pre-fock templates

// both GramMatrix functions construct a dense matrix containing the inner 
// product of every basis element with every other basis element. They split the
// inner products up among (up to) MAX_THREADS different threads, as defined in
// constants.hpp, with each thread doing approximately the same number of 
// elements.
/*template<class T>
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
		if(entryCount > totalCount/numThreads) {
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
	SMatrix gram(basis.size(), basis.size());
	gram.setFromTriplets(entries.begin(), entries.end());
	return gram;
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
		if(entryCount > totalCount/numThreads){
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
	SMatrix gram(totalSize, totalSize);
	gram.setFromTriplets(entries.begin(), entries.end());
	return gram;
}*/

#endif
