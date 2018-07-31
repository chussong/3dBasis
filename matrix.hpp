#ifndef MATRIX_HPP
#define MATRIX_HPP

// This header is for functions related to computing matrix elements from 
// polynomials; perhaps it should be split into one file for the mass matrix
// and one for the phi^4 term, depending on how different they end up being.

#include <vector>
#include <cmath>
#include <unordered_map> // for caching integral results
#include <algorithm> // std::remove_if
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_gamma.h> // beta function
#include <boost/functional/hash.hpp>

#include "constants.hpp"
#include "multinomial.hpp"
#include "mono.hpp"
#include "basis.hpp"
#include "io.hpp"
#include "discretization.hpp"

// these should be the only functions you have to call from other files -------

hp_class InnerFock(const Mono& A, const Mono& B);
HPMatrix GramFock(const Basis<Mono>& basis);
hp_class InnerProduct(const Mono& A, const Mono& B);
DMatrix GramMatrix(const Basis<Mono>& basis, const std::size_t partitions);
DMatrix MassMatrix(const Basis<Mono>& basis, const std::size_t partitions);
DMatrix KineticMatrix(const Basis<Mono>& basis, const std::size_t partitions);
DMatrix InteractionMatrix(const Basis<Mono>& basis, const std::size_t partitions);
DMatrix NPlus2Matrix(const Basis<Mono>& basisA, const Basis<Mono>& basisB,
                     const std::size_t partitions);

// internal stuff -------------------------------------------------------------

namespace MatrixInternal {

// this seems to be necessary if I've declared any operator<< inside of
// MatrixInternal, which I've done for YTerm and MatrixTerm_Intermediate
using ::operator<<;

// the main point of this header
DMatrix Matrix(const Basis<Mono>& basis, const std::size_t partitions, 
        const MATRIX_TYPE type);
coeff_class MatrixTerm(const Mono& A, const Mono& B, const MATRIX_TYPE type);
DMatrix MatrixBlock(const Mono& A, const Mono& B, const MATRIX_TYPE type,
        const std::size_t partitions);

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
    coeff_class coeff = 1;
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
    hp_class coeff = 1;
    std::vector<char> uPlus;
    std::vector<char> uMinus;
    std::vector<char> sinTheta;
    std::vector<char> cosTheta;

    MatrixTerm_Final() = default;
    explicit MatrixTerm_Final(const size_t n);
    // maybe should be rvalue references instead? hopefully it compiles the same
    MatrixTerm_Final(const coeff_class coeff, 
                    const std::vector<char>& uPlus, const std::vector<char>& uMinus, 
                    const std::vector<char>& sinTheta, const std::vector<char>& cosTheta);
    void Resize(const size_t n);
};

struct InteractionTerm_Step2 {
    // coeff(F1) * coeff(F2), not including degeneracies or coeffs of A&B (yet?)
    coeff_class coeff;
    // u1+, u1-, u2+, u2-, ... , u(n-1)+, u(n-1)-, u'(n-1)+, u'(n-1)-
    std::vector<char> u;
    // sin(theta_1), cos(theta_1), ... , sin(theta_(n-2)), cos(theta_(n-2))
    std::vector<char> theta;
    // r, sqrt(1 - r^2), sqrt(1 - alpha^2 r^2)
    std::array<char,3> r;
    // the overall power of alpha derived from all yTilde' except the last one
    char alpha;
};
std::ostream& operator<<(std::ostream& os, const InteractionTerm_Step2& out);

struct InteractionTerm_Output {
    // entire constant part of term, including prefactors, degeneracies, etc.
    coeff_class coeff;
    // r, sqrt(1 - r^2), sqrt(1 - alpha^2 r^2)
    std::array<char,3> r;

    InteractionTerm_Output(coeff_class coeff, std::array<char,3> r):
        coeff(coeff), r(r) {}
};

typedef std::unordered_map<std::array<char,2>, coeff_class,
                           boost::hash<std::array<char,2>> >NtoN_Final;

struct NPlus2Term_Step2 {
    // coeff(F1) * coeff(F2), not including degeneracies or coeffs of A&B (yet?)
    coeff_class coeff;
    // u1+, u1-, u2+, u2-, ..., u(n-1)+, u(n-1)-, u'(n-1)+, u'(n-1)-, u'n+, u'n-
    std::vector<char> u;
    // sin(theta_1), cos(theta_1), ... , sin(theta_(n-2)), cos(theta_(n-2)),
    // sin(theta'), cos(theta')
    std::vector<char> theta;
    // r, sqrt(1 - r^2), sqrt(1 - alpha^2 r^2)
    char r;
    char alpha;
};
std::ostream& operator<<(std::ostream& os, const NPlus2Term_Step2& out);

struct NPlus2Term_Output {
    // entire constant part of term, including prefactors, degeneracies, etc.
    coeff_class coeff;
    // exponent of r
    char r;
    // exponent of alpha
    char alpha;

    NPlus2Term_Output(const coeff_class coeff, const char r, const char alpha):
        coeff(coeff), r(r), alpha(alpha) {}
};

// the two functions for actually computing the two types of MatrixTerms:
// direct for inner product and mass, inter for interactions
hp_class MatrixTerm_Direct(
        const Mono& A, const Mono& B, const MATRIX_TYPE type);
NtoN_Final MatrixTerm_NtoN(
        const Mono& A, const Mono& B);
std::vector<NPlus2Term_Output> MatrixTerm_NPlus2(const Mono& A, const Mono& B);

// coordinate transform functions, called from MatrixTerm
std::string ExtractXY(const Mono& extractFromThis);
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
const std::vector<MatrixTerm_Final>& DirectTermsFromXY(const std::string& xAndy);
std::vector<MatrixTerm_Final> CombineTwoFs(const std::vector<MatrixTerm_Final>& F1,
		const std::vector<MatrixTerm_Final>& F2);
hp_class FinalResult(std::vector<MatrixTerm_Final>& exponents,
		const MATRIX_TYPE type);

// functions specific to INTERACTION computations
const std::vector<MatrixTerm_Intermediate>& InteractionTermsFromXY(
        const std::string& xAndy);
std::vector<InteractionTerm_Step2> CombineInteractionFs(
        const std::vector<MatrixTerm_Intermediate>& F1, 
        const std::vector<MatrixTerm_Intermediate>& F2 );
InteractionTerm_Step2 CombineInteractionFs_OneTerm(
        const MatrixTerm_Intermediate& F1, const MatrixTerm_Intermediate& F2 );
NtoN_Final InteractionOutput(const std::vector<InteractionTerm_Step2>& combinedFs, 
                             const coeff_class prefactor);
const NtoN_Final& Expand(const std::array<char,3>& r, const char alpha);

std::vector<NPlus2Term_Step2> CombineNPlus2Fs(
        const std::vector<MatrixTerm_Intermediate>& F1, 
        const std::vector<MatrixTerm_Intermediate>& F2);
NPlus2Term_Step2 CombineNPlus2Fs_OneTerm(
        const MatrixTerm_Intermediate& f1, const MatrixTerm_Intermediate& f2);
std::vector<NPlus2Term_Output> NPlus2Output(
        std::vector<NPlus2Term_Step2>& combinedFs, const coeff_class prefactor);

// numerical prefactors used in the various computations
hp_class Prefactor(const Mono& A, const Mono& B, const MATRIX_TYPE type);
hp_class InnerProductPrefactor(const char n);
hp_class MassMatrixPrefactor(const char n);
hp_class InteractionMatrixPrefactor(const char n);
hp_class NPlus2MatrixPrefactor(const char n);

// integrals used in FinalResult
hp_class DoAllIntegrals(const MatrixTerm_Final& term);
hp_class DoAllIntegrals(const InteractionTerm_Step2& term);
hp_class DoAllIntegrals_NPlus2(const NPlus2Term_Step2& term);
hp_class UPlusIntegral(const builtin_class a, const builtin_class b);
hp_class ThetaIntegral_Short(const builtin_class a, const builtin_class b);
hp_class ThetaIntegral_Long(const builtin_class a, const builtin_class b);

} // namespace MassMatrix
#endif
