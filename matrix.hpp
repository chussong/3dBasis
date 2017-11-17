#ifndef MATRIX_HPP
#define MATRIX_HPP

// This header is for functions related to computing matrix elements from 
// polynomials; perhaps it should be split into one file for the mass matrix
// and one for the phi^4 term, depending on how different they end up being.

#include <vector>
#include <cmath>
#include <unordered_map> // for caching integral results
#include <gsl/gsl_sf_hyperg.h>
#include <boost/functional/hash.hpp>

#include "constants.hpp"
#include "binomial.hpp"
#include "multinomial.hpp"
#include "mono.hpp"
#include "basis.hpp"
#include "io.hpp"

coeff_class InnerFock(const Mono& A, const Mono& B);
DMatrix GramFock(const Basis<Mono>& basis);
DMatrix MassMatrix(const Basis<Mono>& basis);

namespace MatrixInternal {

// the main point of this header
coeff_class MassMatrixTerm(const Mono& A, const Mono& B, 
		const bool isInnerProduct);

// two structs used in the coordinate transformations for MassMatrixTerm
struct MatrixTerm_Intermediate {
	coeff_class coefficient = 1;
	std::vector<char> uPlus;
	std::vector<char> uMinus;
	std::vector<char> yTilde;

	MatrixTerm_Intermediate() = default;
	explicit MatrixTerm_Intermediate(const size_t n);
	void Resize(const size_t n);
};
struct MatrixTerm_Final {
	coeff_class coefficient = 1;
	std::vector<char> uPlus;
	std::vector<char> uMinus;
	std::vector<char> sinTheta;
	std::vector<char> cosTheta;

	MatrixTerm_Final() = default;
	explicit MatrixTerm_Final(const size_t n);
	// maybe should be rvalue references instead? hopefully it compiles the same
	MatrixTerm_Final(const coeff_class coefficient, 
			const std::vector<char>& uPlus, const std::vector<char>& uMinus, 
			const std::vector<char>& sinTheta, const std::vector<char>& cosTheta);
	void Resize(const size_t n);
};

// coordinate transform functions, called from MassMatrixTerm
std::vector<MatrixTerm_Final> MatrixTermsFromMono(const Mono& input);
std::vector<MatrixTerm_Final> MatrixTermsFromMono_Permuted(const Mono& input,
		const std::vector<std::size_t>& permutationVector);
std::vector<MatrixTerm_Final> MatrixTermsFromXandY(
		const std::array<std::vector<char>,2>& xAndy, const int nParticles);
std::array<std::vector<char>,2> ExponentExtractXY(const Mono& extractFromThis);
std::array<std::vector<char>,2> PermuteXandY(
		const std::array<std::vector<char>,2>& xAndy,
		const std::vector<std::size_t>& permutationVector);
std::vector<char> ExponentUFromX(const std::vector<char>& x);
std::vector<MatrixTerm_Intermediate> ExponentYTildeFromY(const std::vector<char>& y);
std::vector<MatrixTerm_Final> ExponentThetaFromYTilde(
		std::vector<MatrixTerm_Intermediate>& intermediateTerms);

// coordinate transform helper functions, called from transforms
MatrixTerm_Intermediate YTildeTerm(const unsigned int i, const char a, 
		const char l, const std::vector<char>& mVector);
MatrixTerm_Intermediate YTildeLastTerm(const unsigned int n, const char a, 
		const char l, const std::vector<char>& mVector);
coeff_class YTildeCoefficient(const char a, const char l, 
		const std::vector<char>& mVector);
coeff_class YTildeLastCoefficient(const char a, const char l, 
		const std::vector<char>& mVector);

// functions representing individual steps of the MassMatrixTerm computation
std::vector<char> AddVectors(const std::vector<char>& A, 
		const std::vector<char>& B);
std::vector<MatrixTerm_Final> CombineTwoFs(const std::vector<MatrixTerm_Final>& F1,
		const std::vector<MatrixTerm_Final>& F2);
/*std::vector<char> AddVectors(const std::vector<char>& A, 
		std::vector<char> B, const std::vector<std::size_t>& permVector);
std::vector<MatrixTerm_Final> CombineTwoFs(const std::vector<MatrixTerm_Final>& F1,
		const std::vector<MatrixTerm_Final>& F2, std::vector<std::size_t>& perm);*/
coeff_class FinalResult(const std::vector<MatrixTerm_Final>& exponents,
		const bool isInnerProduct);

// prefactor and integrals used in FinalResult
coeff_class InnerProductPrefactor(const char n);
coeff_class MassMatrixPrefactor(const char n);
coeff_class UIntegral(const coeff_class a, const coeff_class b);
coeff_class ThetaIntegral_Short(const coeff_class a, const coeff_class b);
coeff_class ThetaIntegral_Long(const coeff_class a, const coeff_class b);

} // namespace MassMatrix
#endif
