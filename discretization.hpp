#ifndef DISCRETIZATION_HPP
#define DISCRETIZATION_HPP

#include <vector>
#include <cmath>
#include <unordered_map>

#include <boost/functional/hash.hpp>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_errno.h>

#include "constants.hpp"
#include "mono.hpp"
#include "basis.hpp"

// struct DiscMono {
    // const Mono& mono;
    // const std::size_t slice;
// 
    // DiscMono(const Mono& mono, const std::size_t k): mono(mono), slice(k) {}
// };
// 
// class DiscBasis {
    // Basis<Mono> basis;
    // std::size_t partitions;
// 
    // public:
        // DiscMono operator()(const std::size_t m, const std::size_t k) {
            // return DiscMono(basis[m], k);
        // }
// };

class InteractionCache {
    std::size_t partitions;
    coeff_class partWidth;
    std::unordered_map< std::array<char,3>, DMatrix, 
        boost::hash<std::array<char,3>> > cache;

    public:
        void SetPartitions(const std::size_t partitions, 
                const coeff_class partWidth);
        bool HasPartitions(const std::size_t partitions,
                const coeff_class partWidth);

        void Emplace(const std::array<char,3>& key, DMatrix value);
        bool Contains(const std::array<char,3>& key);
        const DMatrix& operator[](const std::array<char,3>& key) const;
        void Clear();
};

char MuExponent(const Mono& A, const Mono& B);
DVector MuIntegral(const Mono& A, const Mono& B, const std::size_t partitions,
        const coeff_class partitionWidth, const MATRIX_TYPE calculationType);
DVector MuIntegral_Body(const char muExp, const std::size_t partitions,
                const coeff_class partitionWidth);
coeff_class MuNorm(const Mono& A, const std::size_t k, 
        const coeff_class partWidth);
DMatrix PartitionMu_Mass(const Basis<Mono>& minimalBasis, const DMatrix& mass,
		const std::size_t partitions, const coeff_class partWidth);

DMatrix MuPart(const Mono& A, const Mono& B, const std::size_t partitions, 
        const coeff_class partWidth, const MATRIX_TYPE type);
DMatrix MuPart_Kinetic(const std::size_t partitions, 
        const coeff_class partWidth);
DMatrix MuTotal(const Basis<Mono>& minBasis, const std::size_t partitions,
        const coeff_class partWidth, const MATRIX_TYPE calculationType);

DMatrix InteractionMu(const std::array<char,3> r, 
        const std::size_t partitions, const coeff_class partWidth);
coeff_class RIntegral(const char a, const char b, const char c, 
        const coeff_class alpha);
coeff_class Hypergeometric2F1(const coeff_class a, const coeff_class b,
        const coeff_class c, const coeff_class x);

DMatrix DiscretizeMonos(const Basis<Mono>& minBasis, 
        const std::size_t partitions, const coeff_class partWidth);
DMatrix DiscretizePolys(const DMatrix& polysOnMinBasis, 
        const std::size_t partitions);

#endif
