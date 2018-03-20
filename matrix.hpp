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
#include "multinomial.hpp"
#include "mono.hpp"
#include "basis.hpp"
#include "io.hpp"
#include "discretization.hpp"

// these should be the only functions you have to call from other files -------

coeff_class InnerFock(const Mono& A, const Mono& B);
DMatrix GramFock(const Basis<Mono>& basis, const std::size_t partitions, 
        const coeff_class partWidth);
DMatrix MassMatrix(const Basis<Mono>& basis, const std::size_t partitions, 
        const coeff_class partWidth);
DMatrix InteractionMatrix(const Basis<Mono>& basis, const std::size_t partitions, 
        const coeff_class partWidth);

// internal stuff -------------------------------------------------------------

namespace MatrixInternal {

// this seems to be necessary if I've declared any operator<< inside of
// MatrixInternal, which I've done for YTerm and MatrixTerm_Intermediate
using ::operator<<;

// the main point of this header
DMatrix Matrix(const Basis<Mono>& basis, const std::size_t partitions, 
        const coeff_class partWidth, const MATRIX_TYPE type);
coeff_class MatrixTerm(const Mono& A, const Mono& B, const MATRIX_TYPE type);

// five structs used in the coordinate transformations for MatrixTerm

struct YTerm {
	coeff_class coeff;
	std::string y;

	YTerm() = delete;
	YTerm(const coeff_class coeff, const std::string& y, 
			const std::string& nAndm);
	
	char operator[](std::size_t i) const { return y[i]; }
	std::size_t size() const { return y.size(); }
};
std::ostream& operator<<(std::ostream& os, const YTerm& out);

struct MatrixTerm_Intermediate {
	coeff_class coefficient = 1;
	std::vector<char> uPlus;
	std::vector<char> uMinus;
	std::vector<char> yTilde;

	MatrixTerm_Intermediate() = default;
	explicit MatrixTerm_Intermediate(const size_t n);
	void Resize(const size_t n);
};
MatrixTerm_Intermediate operator*(const MatrixTerm_Intermediate& A,
			MatrixTerm_Intermediate B);
std::ostream& operator<<(std::ostream& os, const MatrixTerm_Intermediate& out);

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

struct InteractionTerm_Step2 {
    // coeff(F1) * coeff(F2), not including degeneracies or coeffs of A&B (yet?)
    coeff_class coefficient;
    // u1+, u1-, u2+, u2-, ... , u(n-1)+, u(n-1)-, u'(n-1)+, u'(n-1)-
    std::vector<char> u;
    // sin(theta_1), cos(theta_1), ... , sin(theta_(n-2)), cos(theta_(n-2))
    std::vector<char> theta;
    // r, sqrt(1 - r^2), sqrt(1 - alpha^2 r^2)
    std::array<char,3> r;
};
std::ostream& operator<<(std::ostream& os, const InteractionTerm_Step2& out);

struct InteractionTerm_Output {
    // entire constant part of term, including prefactors, degeneracies, etc.
    coeff_class coefficient;
    // r, sqrt(1 - r^2), sqrt(1 - alpha^2 r^2)
    std::array<char,3> r;

    InteractionTerm_Output(coeff_class coefficient, std::array<char,3> r):
        coefficient(coefficient), r(r) {}
};

// the two functions for actually computing the two types of MatrixTerms:
// direct for inner product and mass, inter for interactions
coeff_class MatrixTerm_Direct(
        const Mono& A, const Mono& B, const MATRIX_TYPE type);
std::vector<InteractionTerm_Output> MatrixTerm_Inter(
        const Mono& A, const Mono& B, const MATRIX_TYPE type);

// coordinate transform functions, called from MatrixTerm
std::vector<MatrixTerm_Final> MatrixTermsFromMono(const Mono& input);
std::vector<MatrixTerm_Final> MatrixTermsFromMono_Permuted(const Mono& input,
		const std::vector<std::size_t>& permutationVector);
std::vector<MatrixTerm_Final> MatrixTermsFromXandY(
		const std::array<std::string,2>& xAndy, const int nParticles);
std::array<std::string,2> ExtractXY(const Mono& extractFromThis);
std::array<std::string,2> PermuteXandY(
		const std::array<std::string,2>& xAndy,
		const std::vector<std::size_t>& permutationVector);
bool PermuteXY(std::string& xAndy);
std::array<std::string,2> CombineXandY(const std::array<std::string,2>& xAndy_A,
		std::array<std::string,2> xAndy_B);

std::vector<char> UFromX(const std::string& x);
std::vector<MatrixTerm_Final> ThetaFromY(const std::string y);
std::vector<MatrixTerm_Intermediate> YTildeFromY(const std::string& y);
std::vector<MatrixTerm_Final> ThetaFromYTilde(
		std::vector<MatrixTerm_Intermediate>& intermediateTerms);

// coordinate transform helper functions, called from transforms
std::vector<YTerm> EliminateYn(const std::string& y);
std::vector<MatrixTerm_Intermediate> YTildeTerms(
		const unsigned int i, const char a, const char l, std::string nAndm);
std::vector<MatrixTerm_Intermediate> MultiplyIntermediateTerms(
		const std::vector<MatrixTerm_Intermediate>& termsA, 
		const std::vector<MatrixTerm_Intermediate>& termsB);
//MatrixTerm_Intermediate YTildeLastTerm(const unsigned int n, const char a, 
		//const char l, const std::vector<char>& mVector);
coeff_class YTildeCoefficient(const char a, const char l, 
		const std::string& nAndm);
//coeff_class YTildeLastCoefficient(const char a, const char l, 
		//const std::vector<char>& mVector);

// functions specific to DIRECT computations
std::vector<char> AddVectors(const std::vector<char>& A, 
		const std::vector<char>& B);
std::vector<MatrixTerm_Final> CombineTwoFs(const std::vector<MatrixTerm_Final>& F1,
		const std::vector<MatrixTerm_Final>& F2);
/*std::vector<char> AddVectors(const std::vector<char>& A, 
		std::vector<char> B, const std::vector<std::size_t>& permVector);
std::vector<MatrixTerm_Final> CombineTwoFs(const std::vector<MatrixTerm_Final>& F1,
		const std::vector<MatrixTerm_Final>& F2, std::vector<std::size_t>& perm);*/
coeff_class FinalResult(std::vector<MatrixTerm_Final>& exponents,
		const MATRIX_TYPE type);

// functions specific to INTERACTION computations
const std::vector<MatrixTerm_Intermediate>& InteractionTermsFromXY(
        const std::string& xAndy);
std::vector<InteractionTerm_Step2> CombineInteractionFs(
        const std::vector<MatrixTerm_Intermediate>& F1, 
        const std::vector<MatrixTerm_Intermediate>& F2 );
InteractionTerm_Step2 CombineInteractionFs_OneTerm(
        const MatrixTerm_Intermediate& F1, const MatrixTerm_Intermediate& F2 );
std::vector<InteractionTerm_Output> InteractionOutput(
        std::vector<InteractionTerm_Step2>& combinedFs, 
        const MATRIX_TYPE type, const coeff_class prefactor);

// numerical prefactors used in the various computations
coeff_class Prefactor(const Mono& A, const Mono& B, const MATRIX_TYPE type);
coeff_class PrefactorN(const char n);
coeff_class InnerProductPrefactor(const char n);
coeff_class MassMatrixPrefactor(const char n);
coeff_class InteractionMatrixPrefactor(const char n);

// integrals used in FinalResult
coeff_class DoAllIntegrals(const MatrixTerm_Final& term);
coeff_class DoAllIntegrals(InteractionTerm_Step2& term, const MATRIX_TYPE type);
builtin_class UIntegral(const builtin_class a, const builtin_class b);
builtin_class ThetaIntegral_Short(const builtin_class a, const builtin_class b);
builtin_class ThetaIntegral_Long(const builtin_class a, const builtin_class b);
builtin_class RIntegral(const builtin_class a, const builtin_class alpha);

} // namespace MassMatrix
#endif
