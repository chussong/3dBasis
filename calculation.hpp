#ifndef CALCULATION_HPP
#define CALCULATION_HPP

#include <iostream>
#include <string>
#include <vector>
#include <gsl/gsl_errno.h>  // handling for GSL errors

#include "constants.hpp"
#include "construction.hpp"
#include "test.hpp"
#include "mono.hpp"
#include "poly.hpp"
#include "basis.hpp"
#include "io.hpp"
#include "timer.hpp"
#include "gram-schmidt.hpp"
#include "matrix.hpp"
#include "multinomial.hpp" // the coefficients are initialized in Calculate()
#include "discretization.hpp"

// startup and input parsing --------------------------------------------------

void GSLErrorHandler(const char* reason, const char* file, int line, int err);

// actual computations --------------------------------------------------------

int Calculate(const Arguments& args);
std::vector<Poly> ComputeBasisStates(const Arguments& args);
std::vector<Poly> ComputeBasisStates_SameParity(
        const std::vector<Basis<Mono>>& inputBases, const Arguments& args,
        const bool odd);
DMatrix PolysOnMinBasis(const Basis<Mono>& minimalBasis,
        const std::vector<Poly> orthogonalized, OStream& outStream);
DMatrix ComputeHamiltonian(const Arguments& args);
DMatrix ComputeHamiltonian_SameParity(const std::vector<Basis<Mono>>& inputBases,
                                      const Arguments& args, const bool odd);

#endif
