#ifndef TEST_HPP
#define TEST_HPP

#include <vector>
#include <algorithm> // random_shuffle

#include "constants.hpp"
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

bool RunAllTests(const Arguments& args);

namespace MatrixInternal {

bool PermuteXY(OStream& console);
bool InteractionTermsFromXY(OStream& console);
bool CombineInteractionFs(OStream& console);
bool UPlusIntegral(OStream& console);
bool UPlusIntegral_Case(const builtin_class a, const builtin_class b, 
        const builtin_class expected, OStream& console);

} // namespace MatrixInternal

bool Hypergeometric(OStream& console);
bool InteractionMatrix(const Basis<Mono>& basis, const Arguments& args);

// templates for testing templates --------------------------------------------

template<int P, int Q>
bool HypergeometricPFQ_Case(std::array<builtin_class,P> a, 
        std::array<builtin_class,Q> b, builtin_class x, coeff_class expected,
        OStream& console) {
    constexpr coeff_class tol = 1e-4;
    coeff_class answer = ::HypergeometricPFQ<P,Q>(a, b, x);
    console << "HypergeometricPFQ<" << P << "," << Q << ">(" << a 
        << ", " << b << ", " << x << ") == " << answer;
    if (std::abs(static_cast<builtin_class>(answer - expected)) <= tol*answer) {
        console << " == " << expected << " (PASS)" << endl;
        return true;
    } else {
        console << " != " << expected << " (FAIL)" << endl;
        return false;
    }
}

template<int P, int Q>
bool HypergeometricPFQ_Reg_Case(std::array<builtin_class,P> a, 
        std::array<builtin_class,Q> b, builtin_class x, coeff_class expected,
        OStream& console) {
    constexpr coeff_class tol = 1e-4;
    coeff_class answer = ::HypergeometricPFQ_Reg<P,Q>(a, b, x);
    console << "HypergeometricPFQ_Reg<" << P << "," << Q << ">(" << a 
        << ", " << b << ", " << x << ") == " << answer;
    if (std::abs(static_cast<builtin_class>(answer - expected)) <= tol*answer) {
        console << " == " << expected << " (PASS)" << endl;
        return true;
    } else {
        console << " != " << expected << " (FAIL)" << endl;
        return false;
    }
}

} // namespace Test

#endif
