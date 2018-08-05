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

// these represent alternative ways to do the computation. Implementation 
// details are in gram-schmidt.cpp.
std::vector<Poly> GramSchmidt(const std::vector<Basis<Mono>> input);
std::vector<Poly> GramSchmidt(const Basis<Mono> input);
std::vector<Poly> GramSchmidt_WithMatrix(const std::vector<Basis<Mono>> input,
		const HPMatrix& gramMatrix);
std::vector<Poly> GramSchmidt_WithMatrix(const Basis<Mono> input, 
		const HPMatrix& gramMatrix);
std::vector<Poly> GramSchmidt_WithMatrix_A(const Basis<Mono> inputBasis, 
		const HPMatrix& gramMatrix);
std::vector<Poly> GramSchmidt_WithMatrix_B(const Basis<Mono> inputBasis, 
		const HPMatrix& gramMatrix);
HPVector GSProjection(const HPVector& toProject, const HPVector& projectOnto,
		const HPMatrix& gramMatrix);
hp_class GSNorm(const HPVector& vector, const HPMatrix& gramMatrix);
std::vector<Poly> GramSchmidt_MatrixOnly(const DMatrix& input, 
		const Basis<Mono>& inputBases);

// interface with matrix QR decompositions (as alternative to GS) ------------

std::vector<Poly> PolysFromQMatrix(const DMatrix& QMatrix, 
		const Basis<Mono>& basis, const DMatrix& gramMatrix, 
		const Eigen::Index rank);

// the two solvers need to have their Q matrices extracted in different ways
DMatrix ExtractQMatrix(const Eigen::FullPivHouseholderQR<DMatrix>& solver, 
		               const int dimension);
DMatrix ExtractQMatrix(const Eigen::ColPivHouseholderQR<DMatrix>& solver, 
		               const int dimension);

#endif
