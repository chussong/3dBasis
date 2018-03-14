#ifndef DISCRETIZATION_HPP
#define DISCRETIZATION_HPP

#include <vector>
#include <cmath>

#include "constants.hpp"
#include "mono.hpp"
#include "basis.hpp"

// takes a vector of polynomial states and expands each one into many states 
// with \mu integrated over a small window

// steps:
// 1. determine spectral density \rho(\mu), which is the square inner product of
// O with all \mu_i for which \mu_i^2 == \mu^2. It should be a polynomial in \mu
// normalized so that it gives 1 when integrated over any single \mu window
// 2. each operator given expands into a block of \int d\mu^2 \rho(\mu) f(\mu) O
// where f(\mu) is a function which differs for different Hamiltonian terms. For
// the kinetic term, it's the identity, whereas for the mass term it's \mu^2

char MuExponent(const Mono& A, const Mono& B);
DVector MuIntegral_Mass(const Mono& A, const Mono& B, 
		const std::size_t partitions, const coeff_class partitionWidth);
DVector MuIntegral_InnerProduct(const Mono& A, const Mono& B, 
		const std::size_t partitions, const coeff_class partitionWidth);
DMatrix PartitionMu_Mass(const Basis<Mono>& minimalBasis, const DMatrix& mass,
		const std::size_t partitions, const coeff_class partWidth);

DMatrix DiscretizeMonos(const Basis<Mono>& minBasis, 
        const std::size_t partitions, const coeff_class partWidth);
DMatrix DiscretizePolys(const DMatrix& polysOnMinBasis, 
        const Basis<Mono>& minBasis, const std::size_t partitions);

#endif
