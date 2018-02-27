#ifndef TEST_HPP
#define TEST_HPP

#include <iostream>
#include <algorithm> // random_shuffle

#include "matrix.hpp"
#include "io.hpp"

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

} // namespace Test

#endif
