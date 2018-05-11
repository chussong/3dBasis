#include "discretization.hpp"

// Take a non-discretized polysOnMinBasis matrix and return one that expresses
// the k'th slice of each polynomial in terms of the k'th slices of its 
// constituent monomials
DMatrix DiscretizePolys(const DMatrix& polysOnMinBasis, 
        const std::size_t partitions) {
    if (partitions == 0) return polysOnMinBasis;

    DMatrix output(polysOnMinBasis.rows()*partitions,
            polysOnMinBasis.cols()*partitions);
    for (Eigen::Index i = 0; i < polysOnMinBasis.rows(); ++i) {
        for (Eigen::Index j = 0; j < polysOnMinBasis.cols(); ++j) {
            output.block(i*partitions, j*partitions, partitions, partitions) = 
                polysOnMinBasis(i, j)*DMatrix::Identity(partitions, partitions);
        }
    }

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

    if (n > 2) {
        exponents[0] = 2*exponents[0] + n - 3;
        exponents[1] = exponents[1] + n - 3;
    }

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
        intCache.emplace(exponents, std::move(block));
    }

    return intCache[exponents];
}

coeff_class NtoNWindow_Less(const std::array<char,2>& exponents,
                            const std::array<builtin_class,2>& mu1sq_ab,
                            const std::array<builtin_class,2>& mu2sq_ab) {
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
            // FIXME?? below assumes that all infinite terms must exactly cancel
            // each other
            if (exponents[0] != 2) {
                hypergeos -= common * std::tgamma((a-1.0)/2.0) *
                    Hypergeometric3F2_Reg(0.5, 0.5 + r/2.0, (a-1.0)/2.0, 
                                          r/2.0 + 1.0, (a-1.0)/2.0 + 1.0, x);
            }
        }
    }

    if (!std::isfinite(static_cast<builtin_class>(hypergeos))) {
        std::cerr << "Error: NtoNWindow_Less(" << exponents << ", " 
            << mu1sq_ab << ", " << mu2sq_ab << ") not finite." << std::endl;
    }

    return overall * hypergeos;
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
            << musq_ab << ") not finite." << std::endl;
    }

    return overall * hypergeos;
}

coeff_class NtoNWindow_Equal_Term(const std::array<builtin_class,2>& musq_ab,
                                  const builtin_class arg, 
                                  const builtin_class r,
                                  const bool useMuB) {
    // FIXME?? dubious; but we're assuming that if arg is a pole of the gamma
    // function it will just produce an infinity which has to cancel against
    // something
    if (arg <= 0 && arg - std::round(arg) < EPSILON) {
        return 0;
    }

    const builtin_class& msA = musq_ab[0];
    const builtin_class& msB = musq_ab[1];
    auto HGR = &NtoNWindow_Equal_Hypergeometric;
    // output *= HGR(arg, r, 1);
    // output -= std::pow(msA/msB, arg) * std::pow(msA, 1.5) * HGR(arg, r, msA/msB);
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

builtin_class NtoNWindow_Equal_Hypergeometric(const builtin_class arg, 
                                              const builtin_class r,
                                              const builtin_class x) {
    return Hypergeometric3F2_Reg(0.5, (r+1.0)/2.0, arg,
                                 (r+2.0)/2.0, arg + 1, x);
}

const DMatrix& MuPart_NPlus2(const std::array<char,2>& nr, 
                             const std::size_t partitions) {
    if (nr[1]%2 == 1) {
        if (zeroMatrix.count(partitions) == 0) {
            zeroMatrix.emplace(partitions, DMatrix::Zero(partitions, partitions));
        }
        return zeroMatrix[partitions];
    }

    coeff_class partWidth = coeff_class(1) / partitions;
    if (nPlus2Cache.count(nr) == 0) {
        DMatrix block(partitions, partitions);
        for (std::size_t winA = 0; winA < partitions; ++winA) {
            // entry is 0 when alpha > 1, so winB >= winA; when winB == winA, we
            // need to use a special answer as well
            block(winA, winA) = NPlus2Window_Equal(nr[0], nr[1], 
                                                   winA*partWidth, 
                                                   (winA+1)*partWidth);
            for (std::size_t winB = winA+1; winB < partitions; ++winB) {
                block(winA, winB) = NPlus2Window(nr[0], nr[1], 
                        {{static_cast<builtin_class>(winA*partWidth),
                        static_cast<builtin_class>((winA+1)*partWidth)}},
                        {{static_cast<builtin_class>(winB*partWidth), 
                        static_cast<builtin_class>((winB+1)*partWidth)}} );
            }
        }
        nPlus2Cache.emplace(nr, std::move(block));
    }

    return nPlus2Cache[nr];
}

coeff_class NPlus2Window(const char n, const char r, 
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
        std::cerr << "Error: NPlus2Window(" << (int)n << ", " << a 
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

coeff_class Hypergeometric2F1(const builtin_class a, const builtin_class b,
        const builtin_class c, const builtin_class x) {
    static std::unordered_map< std::array<builtin_class,4>,coeff_class,
                          boost::hash<std::array<builtin_class,4>> > hg2f1Cache;

    const std::array<builtin_class,4> params = {{a, b, c, x}};
    if (hg2f1Cache.count(params) == 0) {
        hg2f1Cache.emplace(params, HypergeometricPFQ<2,1>({{a,b}}, {{c}}, x));
    }

    return hg2f1Cache[params];
}

coeff_class Hypergeometric3F2_Reg(const builtin_class a1, 
                                  const builtin_class a2,
                                  const builtin_class a3,
                                  const builtin_class b1,
                                  const builtin_class b2,
                                  const builtin_class x) {
    return Hypergeometric3F2_Reg({{a1, a2, a3, b1, b2, x}});
}

coeff_class Hypergeometric3F2_Reg(const std::array<builtin_class,3>& a, 
        const std::array<builtin_class,2>& b, const builtin_class x) {
    return Hypergeometric3F2_Reg({{a[0], a[1], a[2], b[0], b[1], x}});
}

coeff_class Hypergeometric3F2_Reg(const std::array<builtin_class,6>& params) {
    static std::unordered_map<std::array<builtin_class,6>,coeff_class,
                          boost::hash<std::array<builtin_class,6>> > hgfrCache;

    if (hgfrCache.count(params) == 0) {
        // coeff_class reg = std::tgamma(b[0]) * std::tgamma(b[1]);
        // hgfrCache.emplace(params, Hypergeometric3F2(a, b, x) / reg);
        coeff_class value;
        try {
            value = HypergeometricPFQ_Reg<3,2>({{params[0], params[1], 
                                                 params[2]}},
                                               {{params[3], params[4]}}, 
                                               params[5]);
            if (!std::isfinite(static_cast<builtin_class>(value))) {
                std::cerr << "Error: 3F2(" << params << ") = " << value << '\n';
            }
        }
        catch (const std::runtime_error& err) {
            std::cerr << "Error: 3F2(" << params << ") did not converge.\n";
            value = 0.0/0.0;
        }
        // std::cout << "Hypergeometric3F2_Reg(" << params << ") = " <<  value 
            // << '\n';
        hgfrCache.emplace(params, value);
    }

    return hgfrCache[params];
}
