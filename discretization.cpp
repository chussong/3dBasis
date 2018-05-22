#include "discretization.hpp"

// Take a non-discretized polysOnMinBasis matrix and return one that expresses
// the k'th slice of each polynomial in terms of the k'th slices of its 
// constituent monomials
SMatrix DiscretizePolys(const DMatrix& polysOnMinBasis, 
                        std::size_t partitions) {
    // if (partitions == 0) return polysOnMinBasis;
    if (partitions == 0) partitions = 1;

    SMatrix output(polysOnMinBasis.rows()*partitions,
                   polysOnMinBasis.cols()*partitions);
    std::vector<Triplet> triplets;
    for (Eigen::Index i = 0; i < polysOnMinBasis.rows(); ++i) {
        for (Eigen::Index j = 0; j < polysOnMinBasis.cols(); ++j) {
            for (std::size_t p = 0; p < partitions; ++p) {
                triplets.emplace_back(i*partitions + p, 
                                      j*partitions + p, 
                                      polysOnMinBasis(i, j));
            }
            // output.block(i*partitions, j*partitions, partitions, partitions) = 
                // polysOnMinBasis(i, j)*DMatrix::Identity(partitions, partitions);
        }
    }

    output.setFromTriplets(triplets.begin(), triplets.end());
    return output;
}

// direct matrix computations -------------------------------------------------

// given monomials A and B, give the PxP tensor representing block A,B of the
// specified matrix
//
// this works for every matrix computation except the interactions, because 
// they pass different arguments
DMatrix MuPart(const std::size_t partitions, const MATRIX_TYPE type) {
    if (type == MAT_INNER || type == MAT_MASS) {
        return DMatrix::Identity(partitions, partitions);
    } else if (type == MAT_KINETIC) {
        return MuPart_Kinetic(partitions);
    } else {
        std::cerr << "Error: the requested MuPart type has not yet been "
            << "implemented." << std::endl;
        throw std::logic_error("MuPart: type not implemented");
    }
}

DMatrix MuPart_Kinetic(const std::size_t partitions) {
    coeff_class partWidth = coeff_class(1) / partitions;
    DMatrix output = DMatrix::Zero(partitions, partitions);
    for (std::size_t k = 0; k < partitions; ++k) {
        output(k, k) = partWidth/2;
        output(k, k) *= k + (k+1);
    }
    // std::cout << "KINETIC BLOCK:\n" << output << std::endl;
    return output;
}

// interaction (same n) matrix computations -----------------------------------

namespace {
    std::unordered_map<std::array<char,2>, DMatrix, 
        boost::hash<std::array<char,2>> > nPlus2Cache;

    std::unordered_map<std::size_t,DMatrix> zeroMatrix;
    std::unordered_map<std::size_t,DMatrix> nEquals2Matrix;

    // caches for expensive functions
    std::unordered_map<std::array<builtin_class,3>,coeff_class,
        boost::hash<std::array<builtin_class,3>> > betaCache;
} // anonymous namespace

// before transformation, first exponent is that of alpha, and the second is 
// that of r; afterward, the first is the exponent of sqrt(alpha), and the
// second is the exponent of r
const DMatrix& MuPart_NtoN(const unsigned int n,
                           std::array<char,2> exponents, 
                           const std::size_t partitions) {
    static std::unordered_map<std::array<char,2>, DMatrix, 
                              boost::hash<std::array<char,2>> > intCache;

    if (n == 2) return MuPart_2to2(partitions);

    exponents[0] = 2*exponents[0] + n - 3;
    exponents[1] = exponents[1] + n - 3;

    builtin_class partWidth = builtin_class(1) / partitions;
    if (intCache.count(exponents) == 0) {
        DMatrix block(partitions, partitions);
        for (std::size_t winA = 0; winA < partitions; ++winA) {
            block(winA, winA) = NtoNWindow_Equal(exponents,
                                                 {{winA*partWidth, 
                                                   (winA+1)*partWidth}} );
            for (std::size_t winB = winA+1; winB < partitions; ++winB) {
                std::array<builtin_class,2> mu1sq_ab{{winA*partWidth, 
                                                      (winA+1)*partWidth}};
                std::array<builtin_class,2> mu2sq_ab{{winB*partWidth, 
                                                      (winB+1)*partWidth}};
                block(winA, winB) = NtoNWindow_Less(exponents, mu1sq_ab,
                                                    mu2sq_ab);
                block(winB, winA) = NtoNWindow_Greater(exponents, mu2sq_ab,
                                                       mu1sq_ab);
            }
        }
        // this is from the normalization of the g_k
        block /= 2*M_PI*partWidth;
        intCache.emplace(exponents, std::move(block));
    }

    return intCache[exponents];
}

const DMatrix& MuPart_2to2(const std::size_t partitions) {
    static std::unique_ptr<DMatrix> output;

    if (output == nullptr) {
        output = std::make_unique<DMatrix>(DMatrix::Zero(partitions, partitions));

        builtin_class width = builtin_class(1)/partitions;
        for (std::size_t i = 0; i < partitions; ++i) {
            (*output)(i, i) = 2*(std::sqrt((i+1)*width) - std::sqrt(i*width));
        }
    }

    return *output;
}

coeff_class NtoNWindow_Less(const std::array<char,2>& exponents,
                            const std::array<builtin_class,2>& mu1sq_ab,
                            const std::array<builtin_class,2>& mu2sq_ab) {
    // if exponents[0] == 2, there's a 0/0 limit that must be done separately
    if (exponents[0] == 2) {
        return NtoNWindow_Less_Special(exponents[1], mu1sq_ab, mu2sq_ab);
    }

    const builtin_class a = exponents[0]/2.0; // exponent of alpha (not alpha^2)
    const builtin_class r = exponents[1];     // exponent of r     (not r^2)
    const coeff_class overall = std::sqrt(M_PI)*std::tgamma(0.5 + r/2.0) / 3.0;

    coeff_class hypergeos = 0;
    for (std::size_t i = 0; i < 2; ++i) {
        builtin_class mu1 = mu1sq_ab[i];
        for (std::size_t j = 0; j < 2; ++j) {
            builtin_class mu2 = mu2sq_ab[j];
            int sign = (i+j)%2 == 0 ? 1 : -1;
            builtin_class x = mu1 / mu2;

            coeff_class common = sign * mu1 * std::sqrt(mu2) * std::pow(x, a/2.0);

            hypergeos += common * std::tgamma((a+2.0)/2.0) *
                Hypergeometric3F2_Reg(0.5, 0.5 + r/2.0, (a+2.0)/2.0,
                                      r/2.0 + 1.0, (a+2.0)/2.0 + 1.0, x);

            hypergeos -= common * std::tgamma((a-1.0)/2.0) *
                Hypergeometric3F2_Reg(0.5, 0.5 + r/2.0, (a-1.0)/2.0, 
                                      r/2.0 + 1.0, (a-1.0)/2.0 + 1.0, x);
        }
    }

    if (!std::isfinite(static_cast<builtin_class>(hypergeos))) {
        std::cerr << "Error: NtoNWindow_Less(" << exponents << ", " 
            << mu1sq_ab << ", " << mu2sq_ab << ") not finite.\n";
    }
    if (hypergeos < 0) {
        std::cerr << "Error: NtoNWindow_Less(" << exponents << ", " 
            << mu1sq_ab << ", " << mu2sq_ab << ") = " << hypergeos << " < 0.\n";
    }

    return overall * hypergeos;
}

// this is the special case of NtoNWindow_Less that happens when a == 1, i.e. 
// there is exactly 1 power of alpha (equivalently, 1/2 power of alpha^2)
coeff_class NtoNWindow_Less_Special(const builtin_class r, 
                                const std::array<builtin_class,2>& mu1sq_ab, 
                                const std::array<builtin_class,2>& mu2sq_ab) {
    // if the intervals are adjacent, there's a term that becomes indeterminate,
    // so we'll just use an approximation instead of the real answer; we'd like
    // to use the trapezoid rule, but that's indeterminate along the boundary as
    // well, so we're reduced to using the value in the middle
    if (mu1sq_ab[1] == mu2sq_ab[0]) {
        builtin_class mu1sq = (mu1sq_ab[1]+mu1sq_ab[0])/2.0;
        builtin_class mu2sq = (mu2sq_ab[1]+mu2sq_ab[0])/2.0;

        coeff_class midValue = std::sqrt(M_PI) * std::sqrt(mu1sq) / (2*mu2sq)
                             * std::tgamma((1.0+r)/2.0)
                             * Hypergeometric2F1_Reg(0.5, 0.5+r/2.0, 
                                                     1.0+r/2.0, mu1sq/mu2sq);

        return midValue * (mu1sq_ab[1]-mu1sq_ab[0]) * (mu2sq_ab[1]-mu2sq_ab[0]);
    }

    coeff_class overall = std::sqrt(M_PI) * std::tgamma((r + 3.0)/2.0)
                        / std::tgamma((r + 4.0)/2.0);

    coeff_class hypergeos = 0;
    for (std::size_t i = 0; i < 2; ++i) {
        builtin_class mu1 = mu1sq_ab[i];
        for (std::size_t j = 0; j < 2; ++j) {
            builtin_class mu2 = mu2sq_ab[j];
            int sign = (i+j)%2 == 0 ? 1 : -1;
            builtin_class x = mu1 / mu2;

            coeff_class common = sign * std::pow(mu1, 1.5);

            hypergeos += (common * 8.0 * (r+2.0)) / (9.0 * (r+1.0))
                       * Hypergeometric2F1(1.5, 0.5 + r/2.0, 1.0 + r/2.0, x);

            hypergeos -= (common * 2.0 * x) / 5.0
                       * Hypergeometric3F2(1.5, 2.5, 1.5 + r/2.0,
                                           3.5, 2.0 + r/2.0, x);

            hypergeos -= (common * 8.0 * x) / 15.0
                       * Hypergeometric3F2(2.5, 2.5, 1.5 + r/2.0,
                                           3.5, 2.0 + r/2.0, x);

            hypergeos -= (common * 0.5 * x)
                       * Hypergeometric4F3(1.0, 1.0, 2.5, 1.5 + r/2.0,
                                           2.0, 2.0, 2.0 + r/2.0, x);
        }
    }

    coeff_class output = std::sqrt(M_PI) * std::tgamma(0.5 + r/2.0)
                       / (3 * std::tgamma(1.0 + r/2.0));
    output *= std::pow(mu1sq_ab[1], 1.5) - std::pow(mu1sq_ab[0], 1.5);
    output *= std::log(mu2sq_ab[1] / mu2sq_ab[0]);
    output += overall * hypergeos;

    if (!std::isfinite(static_cast<builtin_class>(output))) {
        std::cerr << "Error: NtoNWindow_Less_Special(" << r << ", " 
            << mu1sq_ab << ", " << mu2sq_ab << ") not finite.\n";
    }
    if (output < 0) {
        std::cerr << "Error: NtoNWindow_Less_Special(" << r << ", " 
            << mu1sq_ab << ", " << mu2sq_ab << ") = " << hypergeos << " < 0.\n";
    }

    return output;
}

coeff_class NtoNWindow_Greater(const std::array<char,2>& exponents,
                       const std::array<builtin_class,2>& mu1sq_ab,
                       const std::array<builtin_class,2>& mu2sq_ab) {
    return NtoNWindow_Less({{static_cast<char>(2*exponents[1]-exponents[0]), 
                             exponents[1]}}, 
                             mu2sq_ab, mu1sq_ab);
}

coeff_class NtoNWindow_Equal(const std::array<char,2>& exponents,
                             const std::array<builtin_class,2>& musq_ab) {
    const builtin_class a = exponents[0]/2.0; // exponent of alpha (was sqrt(a))
    const builtin_class r = exponents[1];     // exponent of r
    const coeff_class overall = std::sqrt(M_PI)*std::tgamma(0.5 + r/2.0) / 3.0;

    coeff_class hypergeos = 0;
    hypergeos += NtoNWindow_Equal_Term(musq_ab, (a+2.0)/2.0,     r, true );
    hypergeos -= NtoNWindow_Equal_Term(musq_ab, (a-1.0)/2.0,     r, false);
    hypergeos += NtoNWindow_Equal_Term(musq_ab, ((r-a)+2.0)/2.0, r, true );
    hypergeos -= NtoNWindow_Equal_Term(musq_ab, ((r-a)-1.0)/2.0, r, false);

    if (!std::isfinite(static_cast<builtin_class>(hypergeos))) {
        std::cerr << "Error: NtoNWindow_Equal(" << exponents << ", " 
            << musq_ab << ") not finite.\n";
    }
    if (hypergeos < 0) {
        std::cerr << "Error: NtoNWindow_Equal(" << exponents << ", " 
            << musq_ab << ") is negative.\n";
    }

    return overall * hypergeos;
}

coeff_class NtoNWindow_Equal_Term(const std::array<builtin_class,2>& musq_ab,
                                  const builtin_class arg, 
                                  const builtin_class r,
                                  const bool useMuB) {
    const builtin_class& msA = musq_ab[0];
    const builtin_class& msB = musq_ab[1];
    auto HGR = &NtoNWindow_Equal_Hypergeometric;

    if (std::abs(arg) < EPSILON) {
        coeff_class output = (r + 1.0)/4.0;
        output *= Hypergeometric4F3_Reg(1.0, 1.0, 1.5, 1.5 + r/2.0,
                                        2.0, 2.0, 2.0 + r/2.0, 1)
                - msA/msB*Hypergeometric4F3_Reg(1.0, 1.0, 1.5, 1.5 + r/2.0,
                                                2.0, 2.0, 2.0 + r/2.0, msA/msB);
        if (msA != 0.0) output -= std::log(msA/msB) / std::tgamma(1.0 + r/2.0);
        output *= (useMuB ? std::pow(msB, 1.5) : std::pow(msA, 1.5));
        return output;
    }

    if (arg < 0 && std::abs(arg - std::round(arg)) < EPSILON) {
        std::cerr << "Error: NtoNWindow_Equal_Term(" << arg << ") diverges but "
            << "there's no special case for it.\n";
        return 0.0;
    }

    coeff_class output = useMuB ? std::pow(msB, 1.5) : std::pow(msA, 1.5);
    if (msA == 0) {
        // if msA == 0, the second HGR is just 1, so the second term is either
        // 0 or infinity. If it's infinity it's going to have to cancel, so we
        // drop it; if it's zero, it'll be zero either way.
        output *= HGR(arg, r, 1);
    } else {
        output *= HGR(arg, r, 1) - std::pow(msA/msB, arg)*HGR(arg, r, msA/msB);
    }
    return std::tgamma(arg) * output;
}

// this function called "HGR" is the hypergeometric form that appears many times
// in the NtoNWindow_Equal terms
builtin_class NtoNWindow_Equal_Hypergeometric(const builtin_class arg, 
                                              const builtin_class r,
                                              const builtin_class x) {
    return Hypergeometric3F2_Reg(0.5, (r+1.0)/2.0, arg,
                                 (r+2.0)/2.0, arg + 1, x);
}

const DMatrix& MuPart_NPlus2(const std::array<char,2>& nr, 
                             const std::size_t partitions) {
    // if (nr[1]%2 == 1) {
        // if (zeroMatrix.count(partitions) == 0) {
            // zeroMatrix.emplace(partitions, DMatrix::Zero(partitions, partitions));
        // }
        // return zeroMatrix[partitions];
    // }

    coeff_class partWidth = coeff_class(1) / partitions;
    if (nPlus2Cache.count(nr) == 0) {
        DMatrix block = DMatrix::Zero(partitions, partitions);
        for (std::size_t winA = 0; winA < partitions; ++winA) {
            // entry is 0 when alpha > 1, so winB >= winA; when winB == winA, we
            // need to use a special answer as well
            block(winA, winA) = NPlus2Window_Equal(nr[0], nr[1], 
                                                   winA*partWidth, 
                                                   (winA+1)*partWidth);
            for (std::size_t winB = winA+1; winB < partitions; ++winB) {
                block(winA, winB) = NPlus2Window_Less(nr[0], nr[1], 
                        {{static_cast<builtin_class>(winA*partWidth),
                        static_cast<builtin_class>((winA+1)*partWidth)}},
                        {{static_cast<builtin_class>(winB*partWidth), 
                        static_cast<builtin_class>((winB+1)*partWidth)}} );
            }
        }
        // this is from the normalization of the g_k
        block /= 2*M_PI*partWidth;
        nPlus2Cache.emplace(nr, std::move(block));
    }

    return nPlus2Cache[nr];
}

coeff_class NPlus2Window_Less(const char n, const char r, 
        const std::array<builtin_class,2>& mu1_ab,
        const std::array<builtin_class,2>& mu2_ab) {
    builtin_class a = 0.5 * r;
    coeff_class overall = 8.0 / 3.0;

    coeff_class hypergeos = 0;
    for (std::size_t i = 0; i < 2; ++i) {
        coeff_class mu1 = mu1_ab[i];
        for (std::size_t j = 0; j < 2; ++j) {
            coeff_class mu2 = mu2_ab[j];
            int sign = (i+j)%2 == 0 ? 1 : -1;

            builtin_class x = mu1 / mu2;

            coeff_class term = sign * std::pow(mu1, (n+1.0)/4.0) 
                             / std::pow(mu2, (n-5.0)/4.0);
            hypergeos += term * 
                Hypergeometric2F1(-a, (n+1.0)/4.0, (n+5.0)/4.0, x) / (n + 1.0);
            hypergeos -= term * 
                Hypergeometric2F1(-a, (n-5.0)/4.0, (n-1.0)/4.0, x) / (n - 5.0);
        }
    }

    if (!std::isfinite(static_cast<builtin_class>(hypergeos))) {
        std::cerr << "Error: NPlus2Window_Less(" << (int)n << ", " << a 
            << ", " << mu1_ab << ", " << mu2_ab << ") not finite." 
            << std::endl;
    }
    return overall * hypergeos;
}

// This is different from the generic case because it needs to stop when alpha
// is 1, i.e. when mu2 >= mu1; luckily it's still pretty simple
coeff_class NPlus2Window_Equal(const char n, const char r, 
        const builtin_class mu_a, const builtin_class mu_b) {
    builtin_class a = 0.5 * r;
    
    coeff_class gammaPart = std::pow(mu_b, 1.5) * std::tgamma((n+1.0)/4.0) 
                          / std::tgamma((n+5.0)/4.0 + a);
    gammaPart -= std::pow(mu_a, 1.5) * std::tgamma((n-5.0)/4.0) 
               / std::tgamma((n-1.0)/4.0 + a);
    gammaPart *= 2.0 * std::tgamma(a+1.0) / 3.0;

    coeff_class hyperPart = 
        Hypergeometric2F1(-a, (n-5.0)/4.0, (n-1.0)/4.0, mu_a/mu_b) / (n - 5.0);
    hyperPart -= 
        Hypergeometric2F1(-a, (n+1.0)/4.0, (n+5.0)/4.0, mu_a/mu_b) / (n + 1.0);
    hyperPart *= (8.0 * std::pow(mu_a, (n+1.0)/4.0))
               / (3.0 * std::pow(mu_b, (n-5.0)/4.0));

    return gammaPart + hyperPart;
}
