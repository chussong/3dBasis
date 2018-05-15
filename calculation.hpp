#ifndef CALCULATION_HPP
#define CALCULATION_HPP

#include <cctype>
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

// actual computations --------------------------------------------------------

struct Hamiltonian {
    int maxN;
    std::vector<DMatrix> diagonal;
    std::vector<DMatrix> nPlus2;
};

int Calculate(const Arguments& args);
std::vector<Poly> ComputeBasisStates(const Arguments& args);
std::vector<Poly> ComputeBasisStates_SameParity(
        const std::vector<Basis<Mono>>& inputBases, const Arguments& args,
        const bool odd);
DMatrix PolysOnMinBasis(const Basis<Mono>& minimalBasis,
        const std::vector<Poly> orthogonalized, OStream& outStream);
DMatrix ComputeHamiltonian(const Arguments& args);
Hamiltonian FullHamiltonian(Arguments args, const bool odd);
DMatrix DiagonalBlock(const Basis<Mono>& minimalBasis, 
                      const SMatrix& discPolys, 
                      const Arguments& args, const bool odd);
DMatrix NPlus2Block(const Basis<Mono>& basisA, const SMatrix& discPolysA,
                    const Basis<Mono>& basisB, const SMatrix& discPolysB,
                    const Arguments& args, const bool odd);

void AnalyzeHamiltonian(const Hamiltonian& hamiltonian, const Arguments& args);
void AnalyzeHamiltonian_Dense(const Hamiltonian& hamiltonian, 
                               const Arguments& args);
void AnalyzeHamiltonian_Sparse(const Hamiltonian& hamiltonian, 
                               const Arguments& args);

// stuff for printing results -------------------------------------------------

void OutputMatrix(const DMatrix& monoMatrix, const DMatrix& polyMatrix,
                  std::string name, const std::string& suffix, Timer& timer, 
                  const Arguments& args);
std::string MathematicaName(std::string name);

// bookkeeping stuff ----------------------------------------------------------

void GSLErrorHandler(const char* reason, const char* file, int line, int err);

#endif
