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

DMatrix MuPart(const Mono& A, const Mono& B, const std::size_t partitions, 
        const coeff_class partWidth, const MATRIX_TYPE type);
DMatrix MuPart_Kinetic(const std::size_t partitions, 
        const coeff_class partWidth);

const DMatrix& MuPart(const std::array<char,3>& r, 
        const std::size_t partitions, const coeff_class partWidth);
coeff_class RIntegral(const char a, const char b, const char c, 
        const coeff_class alpha);
coeff_class Hypergeometric2F1(const coeff_class a, const coeff_class b,
        const coeff_class c, const coeff_class x);

coeff_class InteractionWindow(const std::array<char,3>& r,
        const std::array<coeff_class,2>& mu1,
        const std::array<coeff_class,2>& mu2);
coeff_class Hypergeometric3F2(const std::array<coeff_class,3>& a, 
        const std::array<coeff_class,2>& b, const coeff_class x);
coeff_class Hypergeometric3F2_Reg(const std::array<coeff_class,3>& a, 
        const std::array<coeff_class,2>& b, const coeff_class x);

DMatrix DiscretizePolys(const DMatrix& polysOnMinBasis, 
        const std::size_t partitions);

#endif
