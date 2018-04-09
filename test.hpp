#ifndef TEST_HPP
#define TEST_HPP

#include <iostream>
#include <algorithm> // random_shuffle

#include "matrix.hpp"
#include "io.hpp"
#include "discretization.hpp"
#include "gram-schmidt.hpp"
#include "hypergeo.hpp"

// This file contains unit tests for various functions; for a function named
// Namespace::Function, the test will be Test::Namespace::Function, and will be
// a bool function with no arguments; each test should inspect a variety of 
// potential arguments to the actual function itself, writing any results to
// cout, and return true if (and only if) all tests are passed.

namespace Test {

bool RunAllTests();

namespace MatrixInternal {

bool PermuteXY();
bool InteractionTermsFromXY();
bool CombineInteractionFs();

} // namespace MatrixInternal

// bool RIntegral();
bool Hypergeometric();
bool InteractionMatrix(const Basis<Mono>& basis, const std::size_t partitions,
        const coeff_class partWidth);

// templates for testing templates --------------------------------------------

template<int P, int Q>
bool HypergeometricPFQ_Case(std::array<builtin_class,P> a, 
        std::array<builtin_class,Q> b, builtin_class x, coeff_class expected) {
    constexpr coeff_class tol = 1e-4;
    coeff_class answer = ::HypergeometricPFQ<P,Q>(a, b, x);
    std::cerr << "HypergeometricPFQ<" << P << "," << Q << ">(" << a 
        << ", " << b << ", " << x << ") == " << answer;
    if (std::abs(static_cast<builtin_class>(answer - expected)) <= tol*answer) {
        std::cerr << " == " << expected << " (PASS)" << std::endl;
        return true;
    } else {
        std::cerr << " != " << expected << " (FAIL)" << std::endl;
        return false;
    }
}

template<int P, int Q>
bool HypergeometricPFQ_Reg_Case(std::array<builtin_class,P> a, 
        std::array<builtin_class,Q> b, builtin_class x, coeff_class expected) {
    constexpr coeff_class tol = 1e-4;
    coeff_class answer = ::HypergeometricPFQ_Reg<P,Q>(a, b, x);
    std::cerr << "HypergeometricPFQ_Reg<" << P << "," << Q << ">(" << a 
        << ", " << b << ", " << x << ") == " << answer;
    if (std::abs(static_cast<builtin_class>(answer - expected)) <= tol*answer) {
        std::cerr << " == " << expected << " (PASS)" << std::endl;
        return true;
    } else {
        std::cerr << " != " << expected << " (FAIL)" << std::endl;
        return false;
    }
}

} // namespace Test

#endif
