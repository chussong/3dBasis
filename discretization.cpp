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
        boost::hash<std::array<char,2>> > intCache;
    std::unordered_map<std::array<char,2>, DMatrix, 
        boost::hash<std::array<char,2>> > nPlus2Cache;

    std::unordered_map<std::size_t,DMatrix> zeroMatrix;
    std::unordered_map<std::size_t,DMatrix> nEquals2Matrix;

    // caches for expensive functions
    std::unordered_map<std::array<builtin_class,6>,coeff_class,
        boost::hash<std::array<builtin_class,6>> > hgfrCache;
    std::unordered_map<std::array<builtin_class,4>,coeff_class,
        boost::hash<std::array<builtin_class,4>> > hg2f1Cache;
    std::unordered_map<std::array<builtin_class,3>,coeff_class,
        boost::hash<std::array<builtin_class,3>> > betaCache;
} // anonymous namespace

const DMatrix& MuPart_NtoN(const unsigned int n,
                           const std::array<char,2>& exponents, 
                           const std::size_t partitions) {
    coeff_class partWidth = coeff_class(1) / partitions;
    if (intCache.count(exponents) == 0) {
        DMatrix block(partitions, partitions);
        for (std::size_t winA = 0; winA < partitions; ++winA) {
            for (std::size_t winB = 0; winB < partitions; ++winB) {
                block(winA, winB) = NtoNWindow(n, exponents, 
                        {{static_cast<builtin_class>(winA*partWidth),
                        static_cast<builtin_class>((winA+1)*partWidth)}},
                        {{static_cast<builtin_class>(winB*partWidth), 
                        static_cast<builtin_class>((winB+1)*partWidth)}} );
            }
        }
        intCache.emplace(exponents, std::move(block));
    }

    return intCache[exponents];
}

// first exponent is that of alpha^2, and the second is that of r (not squared)
// FIXME: this needs a special case for 4*exponents[0] + n == 5
coeff_class NtoNWindow(const unsigned int n,
                       const std::array<char,2>& exponents,
                       const std::array<builtin_class,2>& mu1sq_ab,
                       const std::array<builtin_class,2>& mu2sq_ab) {
    const builtin_class a = exponents[0] + builtin_class(n)/4.0;
    const builtin_class b = exponents[1];

    coeff_class output = 16;
    output *= std::pow(mu1sq_ab[1], a + 0.25) - std::pow(mu1sq_ab[0], a + 0.25);
    output /= (b + 1)*(4*a - 5)*(4*a + 1);

    if (mu2sq_ab[0] == 0) {
        output *= std::pow(mu2sq_ab[1], 1.25 - a);
    } else {
        output *= std::pow(mu2sq_ab[0], 1.25)*std::pow(mu2sq_ab[1], a)
                - std::pow(mu2sq_ab[1], 1.25)*std::pow(mu2sq_ab[0], a);
        output /= std::pow(mu2sq_ab[0]*mu2sq_ab[1], a);
    }

    if (std::isnan(static_cast<builtin_class>(output))) {
        std::cerr << "Error: NtoNWindow(" << n << ", " <<  exponents << ", " 
            << mu1sq_ab << ", " << mu2sq_ab << ") returns a NaN. (a,b) = (" 
            << a << ", " << b << ")" << std::endl;
    }
    return output;
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
            for (std::size_t winB = 0; winB < partitions; ++winB) {
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

            coeff_class term = sign * std::pow(mu1, 0.5) * std::pow(mu2, 1.0);
            hypergeos += term * 
                Hypergeometric2F1(-a, (n+1.0)/4.0, (n+5.0)/4.0, x) / (n + 1.0);
            hypergeos -= term * 
                Hypergeometric2F1(-a, (n-5.0)/4.0, (n-1.0)/4.0, x) / (n - 5.0);
        }
    }

    if (std::isnan(static_cast<builtin_class>(hypergeos))) {
        std::cerr << "Error: NPlus2Window(" << (int)n << ", " << a 
            << ", " << mu1_ab << ", " << mu2_ab << ") returns a NaN." 
            << std::endl;
    }
    return overall * hypergeos;
}

coeff_class NPlus2Window_Less(const char n, const builtin_class a, 
        const std::array<builtin_class,2>& mu1_ab,
        const std::array<builtin_class,2>& mu2_ab) {
    coeff_class overall = 6.0;
    if (static_cast<int>(std::round(a))%2 == 0) overall = -overall;

    coeff_class betas = 0;
    for (std::size_t i = 0; i < 2; ++i) {
        coeff_class mu1 = mu1_ab[i];
        for (std::size_t j = 0; j < 2; ++j) {
            coeff_class mu2 = mu2_ab[j];
            int sign = (i+j)%2 == 0 ? 1 : -1;

            builtin_class x = (mu2*mu2) / (mu1*mu1);

            betas += sign * (mu2*mu2*mu2*Beta({{x, -0.25 - a - 0.25*n, 1.0 + a}})
                    - mu1*mu1*mu1*Beta({{x, 1.25 - a - 0.25*n, 1.0 + a}}));
        }
    }

    if (std::isnan(static_cast<builtin_class>(betas))) {
        std::cerr << "Error: NPLus2Window_Less(" << (int)n << ", " << a 
            << ", " << mu1_ab << ", " << mu2_ab << ") returns a NaN." 
            << std::endl;
    }
    return betas / overall;
}

// all three window types seem to actually be the same
coeff_class NPlus2Window_Equal(const char n, const builtin_class a, 
        const std::array<builtin_class,2>& mu_ab) {
    return NPlus2Window_Less(n, a, mu_ab, mu_ab);
}

coeff_class NPlus2Window_Greater(const char n, const builtin_class a, 
        const std::array<builtin_class,2>& mu1_ab,
        const std::array<builtin_class,2>& mu2_ab) {
    return NPlus2Window_Less(n, a, mu1_ab, mu2_ab);
}

coeff_class Hypergeometric2F1(const builtin_class a, const builtin_class b,
        const builtin_class c, const builtin_class x) {
    const std::array<builtin_class,4> params = {{a, b, c, x}};
    if (hg2f1Cache.count(params) == 0) {
        hg2f1Cache.emplace(params, HypergeometricPFQ<2,1>({{a,b}}, {{c}}, x));
    }

    return hg2f1Cache[params];
}

coeff_class Hypergeometric3F2_Reg(const std::array<builtin_class,3>& a, 
        const std::array<builtin_class,2>& b, const builtin_class x) {
    std::array<builtin_class,6> params = {{a[0], a[1], a[2], b[0], b[1], x}};
    if (hgfrCache.count(params) == 0) {
        coeff_class reg = std::tgamma(b[0]) * std::tgamma(b[1]);
        // hgfrCache.emplace(params, Hypergeometric3F2(a, b, x) / reg);
        hgfrCache.emplace(params, HypergeometricPFQ<3,2>(a, b, x) / reg);
    }

    return hgfrCache[params];
}

// this is the INCOMPLETE BETA FUNCTION, B_args[0](args[1], args[2])
coeff_class Beta(const std::array<builtin_class,3>& args) {
    if (betaCache.count(args) == 0) {
        // need to multiply beta in to go from gsl's regularized IBF -> IBF
        coeff_class beta = std::exp(std::lgamma(args[1]) + std::lgamma(args[2]) 
                - std::lgamma(args[1] + args[2]) );
        betaCache.emplace(args, beta*gsl_sf_beta_inc(args[1], args[2], args[0]));
    }

    return betaCache[args];
}
