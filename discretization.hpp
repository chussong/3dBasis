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

DMatrix DiscretizePolys(const DMatrix& polysOnMinBasis, 
        const std::size_t partitions);

// direct matrices ------------------------------------------------------------

DMatrix MuPart(const std::size_t partitions, const MATRIX_TYPE type);
DMatrix MuPart_Kinetic(const std::size_t partitions);

// same-n interactions --------------------------------------------------------

const DMatrix& MuPart_NtoN(const unsigned int n, 
                           const std::array<char,2>& exponents, 
                           const std::size_t partitions);
coeff_class NtoNWindow(const unsigned int n, 
                       const std::array<char,2>& exponents,
                       const std::array<builtin_class,2>& mu1_ab,
                       const std::array<builtin_class,2>& mu2_ab);

// memoized interfaces to hypergeometric functions
coeff_class Hypergeometric3F2(const std::array<builtin_class,3>& a, 
        const std::array<builtin_class,2>& b, const builtin_class x);
coeff_class Hypergeometric3F2_Reg(const std::array<builtin_class,3>& a, 
        const std::array<builtin_class,2>& b, const builtin_class x);

// n+2 interactions -----------------------------------------------------------

const DMatrix& MuPart_NPlus2(const std::array<char,2>& nr, 
                             const std::size_t partitions);

coeff_class NPlus2Window(const char n, const char r,
        const std::array<builtin_class,2>& mu1_ab,
        const std::array<builtin_class,2>& mu2_ab);
coeff_class NPlus2Window_Equal(const char n, const char r, 
        const builtin_class mu_a, const builtin_class mu_b);

coeff_class Hypergeometric2F1(const builtin_class a, const builtin_class b,
        const builtin_class c, const builtin_class x);
coeff_class Beta(const std::array<builtin_class,3>& args);

#endif
