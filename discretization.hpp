#ifndef DISCRETIZATION_HPP
#define DISCRETIZATION_HPP

#include <array>
#include <vector>
#include <cmath>
#include <unordered_map>
#include <iostream>

#include <boost/functional/hash.hpp>
#include <gsl/gsl_sf_gamma.h> // for beta function

#include "constants.hpp"
#include "hypergeo.hpp"
#include "multinomial.hpp"

SMatrix DiscretizePolys(const DMatrix& polysOnMinBasis, 
                        std::size_t partitions);

// direct matrices ------------------------------------------------------------

DMatrix MuPart(const std::size_t partitions, const MATRIX_TYPE type);
DMatrix MuPart_Kinetic(const std::size_t partitions);

// same-n interactions --------------------------------------------------------

const DMatrix& MuPart_NtoN(const unsigned int n, 
                           std::array<char,2> exponents, 
                           const std::size_t partitions);
const DMatrix& MuPart_2to2(const std::size_t partitions);

coeff_class NtoNWindow_Less(const std::array<char,2>& exponents,
                       const std::array<builtin_class,2>& mu1sq_ab,
                       const std::array<builtin_class,2>& mu2sq_ab);
coeff_class NtoNWindow_Less_Special(const builtin_class r, 
                                const std::array<builtin_class,2>& mu1sq_ab, 
                                const std::array<builtin_class,2>& mu2sq_ab);
coeff_class NtoNWindow_Greater(const std::array<char,2>& exponents,
                       const std::array<builtin_class,2>& mu1sq_ab,
                       const std::array<builtin_class,2>& mu2sq_ab);
coeff_class NtoNWindow_Equal(const std::array<char,2>& exponents,
                             const std::array<builtin_class,2>& musq_ab);
coeff_class NtoNWindow_Equal_Term(const std::array<builtin_class,2>& musq_ab,
                                  const builtin_class arg, 
                                  const builtin_class r,
                                  const bool useMuB);
builtin_class NtoNWindow_Equal_Hypergeometric(const builtin_class arg, 
                                              const builtin_class r,
                                              const builtin_class x);

// n+2 interactions -----------------------------------------------------------

const DMatrix& MuPart_NPlus2(const std::array<char,2>& nr, 
                             const std::size_t partitions);

coeff_class NPlus2Window_Less(const char n, const char r,
        const std::array<builtin_class,2>& mu1_ab,
        const std::array<builtin_class,2>& mu2_ab);
coeff_class NPlus2Window_Equal(const char n, const char r, 
        const builtin_class mu_a, const builtin_class mu_b);
coeff_class NPlus2Window_15_Less(const char n, const char r, 
                                 const std::array<builtin_class,2>& mu1_ab,
                                 const std::array<builtin_class,2>& mu2_ab);
// coeff_class NPlus2Window_1_Less_Term(const builtin_class a,
                                     // const std::array<builtin_class,2>& mu1_ab,
                                     // const std::array<builtin_class,2>& mu2_ab);
coeff_class NPlus2Window_15_Equal(const char n, const char r, 
                                  const builtin_class mu_a, 
                                  const builtin_class mu_b);
// coeff_class NPlus2Window_1_Equal_Term(const builtin_class a,
                                      // const builtin_class mu_a,
                                      // const builtin_class mu_b);

// trapezoid rule approximations for cases that aren't analytically solvable --
coeff_class Trapezoid_Rectangular(
        const std::function<coeff_class(builtin_class,builtin_class)> integrand,
        const std::array<builtin_class,2>& mu1_ab, 
        const std::array<builtin_class,2>& mu2_ab,
        const std::size_t samples);
coeff_class Trapezoid_Triangular(
        const std::function<coeff_class(builtin_class,builtin_class)> integrand,
        const builtin_class mu_a, const builtin_class mu_b,
        const std::size_t samples);
#endif
