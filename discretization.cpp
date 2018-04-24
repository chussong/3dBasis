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
    std::unordered_map<std::array<char,3>, DMatrix, 
        boost::hash<std::array<char,3>> > intCache;
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

// this is the MuPart for the same-n interaction matrix
//
// it is the analytic answer, NOT a trapezoid-style approximation of any kind;
// furthermore, the 3F2s used are mainly (entirely?) just finite polynomials
// due to the presence of negative integers in the p lists
//
// this function returns part of the block of the interaction matrix 
// corresponding to one pair of monomials. The sum of all of the returned blocks
// will be the block's actual value.
const DMatrix& MuPart(const std::array<char,3>& r, 
        const std::size_t partitions) {
    if (r[0]%2 == 1) {
        if (zeroMatrix.count(partitions) == 0) {
            zeroMatrix.emplace(partitions, DMatrix::Zero(partitions, partitions));
        }
        return zeroMatrix[partitions];
    }

    coeff_class partWidth = coeff_class(1) / partitions;
    if (intCache.count(r) == 0) {
        DMatrix block(partitions, partitions);
        for (std::size_t winA = 0; winA < partitions; ++winA) {
            for (std::size_t winB = 0; winB < partitions; ++winB) {
                block(winA, winB) = InteractionWindow(r, 
                        {{static_cast<builtin_class>(winA*partWidth),
                        static_cast<builtin_class>((winA+1)*partWidth)}},
                        {{static_cast<builtin_class>(winB*partWidth), 
                        static_cast<builtin_class>((winB+1)*partWidth)}} );
            }
        }
        intCache.emplace(r, std::move(block));
    }

    return intCache[r];
}

// this is the special MuPart block for n == 2, which must be treated separately
// because there is no r integral in this case
const DMatrix& MuPart_NEquals2(const std::size_t partitions) {
    if (nEquals2Matrix.count(partitions) == 0) {
        DMatrix newMatrix = DMatrix::Zero(partitions, partitions);
        for (Eigen::Index row = 0; row < newMatrix.rows(); ++row) {
            for (Eigen::Index col = 0; col < newMatrix.cols(); ++col) {
                // FIXME: this is probably bogus. I'm just doing the same thing
                // as n > 2 but with an integrand of 1 so it's (\Delta \mu)^2
                newMatrix(row, col) = (coeff_class(1) / partitions);
                newMatrix(row, col) *= newMatrix(row, col);
            }
        }
        nEquals2Matrix.emplace(partitions, newMatrix);
    }

    return nEquals2Matrix[partitions];
}

coeff_class InteractionWindow(const std::array<char,3>& r,
        const std::array<builtin_class,2>& mu1_ab,
        const std::array<builtin_class,2>& mu2_ab) {
    builtin_class a = r[0]/2.0;
    builtin_class b = r[1]/2.0;
    builtin_class c = r[2]/2.0;

    if (mu1_ab[0] < mu2_ab[0]) {
        return InteractionWindow_Less(a, b, c, mu1_ab, mu2_ab);
    } else if (mu1_ab[0] == mu2_ab[0]) {
        return InteractionWindow_Equal(a, b, c, mu1_ab);
    } else {
        return InteractionWindow_Greater(a, b, c, mu1_ab, mu2_ab);
    }
}

coeff_class InteractionWindow_Less(const builtin_class a, const builtin_class b, 
        const builtin_class c, const std::array<builtin_class,2>& mu1_ab,
        const std::array<builtin_class,2>& mu2_ab) {
    coeff_class overall = std::tgamma(static_cast<builtin_class>(1.0 + b));
    overall /= std::tgamma(static_cast<builtin_class>(0.5 - a));
    overall *= std::pow(M_PI, 1.5) / 4;
    if (static_cast<int>(2*a)%4 == 2) overall = -overall;

    coeff_class hypergeos = 0;
    for (coeff_class mu1 : mu1_ab) {
        for (std::size_t i = 0; i < 2; ++i) {
            coeff_class mu2 = mu2_ab[i];
            int sign = i%2 == 0 ? -1 : 1;
            hypergeos += sign * mu1*mu2 * Hypergeometric3F2_Reg(
                    {{0.5+a, -c, -0.5}}, {{1.5+a+b, 1.5}}, (mu1*mu1)/(mu2*mu2));
        }
    }

    if (std::isnan(static_cast<builtin_class>(hypergeos))) {
        std::cerr << "Error: InteractionWindow_Less(" << a << ", " << b 
            << ", " << c << ", " << mu1_ab << ", " << mu2_ab 
            << ") returns a NaN." << std::endl;
    }
    return overall * hypergeos;
}

coeff_class InteractionWindow_Equal(const builtin_class a, const builtin_class b, 
        const builtin_class c, const std::array<builtin_class,2>& mu_ab) {
    coeff_class overall = M_PI / 8.0;
    overall /= std::tgamma(static_cast<builtin_class>(0.5 - a));
    if (static_cast<int>(2*a)%4 == 2) overall = -overall;

    coeff_class x = (mu_ab[0]*mu_ab[0]) / (mu_ab[1]*mu_ab[1]);

    coeff_class groupA = 0;
    groupA -= (mu_ab[1]-mu_ab[0])*(mu_ab[1]+mu_ab[0])*Hypergeometric3F2_Reg(
            {{-0.5, 0.5+a, -c}}, {{0.5, 1.5+a+b}}, 1);
    groupA += mu_ab[1]*mu_ab[1]*Hypergeometric3F2_Reg(
            {{-0.5, 0.5+a, -c}}, {{1.5, 1.5+a+b}}, 1);
    groupA -= mu_ab[0]*mu_ab[1]*Hypergeometric3F2_Reg(
            {{-0.5, 0.5+a, -c}}, {{1.5, 1.5+a+b}}, x);
    groupA *= 2 * std::sqrt(M_PI);
    groupA *= std::tgamma(static_cast<builtin_class>(1.0 + b));
    if (std::isnan(static_cast<builtin_class>(groupA))) {
        std::cerr << "Error: InteractionWindow_Equal(" << a << ", " << b << ", "
            << c << ", " << mu_ab << ") has a NaN in groupA." << std::endl;
    }
    
    coeff_class groupB = 0;
    groupB -= mu_ab[0]*mu_ab[0]*Hypergeometric3F2_Reg(
            {{a, 0.5+a, -b}}, {{2.0+a, 1.5+a+c}}, 1);
    // or {{a, 1.5, 2.0+a+b}}, {{2.0+a, 3.0+b+c}}, 1
    // or {{2.0, 1.0+c, 3.0+b+c}}, {{3.5+a+b+c, 3.0+c}}, 1
    groupB += a*(mu_ab[1]-mu_ab[0])*(mu_ab[1]+mu_ab[0])*Hypergeometric3F2_Reg(
            {{0.5+a, 1.0+a, -b}}, {{2.0+a, 1.5+a+c}}, 1);
    groupB += mu_ab[0]*mu_ab[0]*std::pow(static_cast<builtin_class>(x), a)
        * Hypergeometric3F2_Reg({{a, 0.5+a, -b}}, {{2.0+a, 1.5+a+c}}, x);
    // FIXME: the below gamma function diverges when a=0. Is this valid input?
    groupB *= std::tgamma(static_cast<builtin_class>(a));
    groupB *= std::tgamma(static_cast<builtin_class>(1.0 + c));
    if (std::isnan(static_cast<builtin_class>(groupB))) {
        std::cerr << "Error: InteractionWindow_Equal(" << a << ", " << b << ", "
            << c << ", " << mu_ab << ") has a NaN in groupB." << std::endl;
    }

    return overall * (groupA + groupB);
}

coeff_class InteractionWindow_Greater(const builtin_class a, const builtin_class b, 
        const builtin_class c, const std::array<builtin_class,2>& mu1_ab,
        const std::array<builtin_class,2>& mu2_ab) {
    if (a == 0) return InteractionWindow_Greater_Zero(b, c, mu1_ab, mu2_ab);
    coeff_class overall = std::pow(4, static_cast<builtin_class>(-a-1.0));
    overall *= std::sqrt(M_PI);
    overall *= std::tgamma(static_cast<builtin_class>(2*a));
    overall *= std::tgamma(static_cast<builtin_class>(1.0 + c));

    coeff_class hypergeos = 0;
    for (std::size_t i = 0; i < 2; ++i) {
        coeff_class mu1 = mu1_ab[i];
        for (std::size_t j = 0; j < 2; ++j) {
            coeff_class mu2 = mu2_ab[j];
            coeff_class x = (mu2*mu2)/(mu1*mu1);
            int sign = (i+j)%2 == 0 ? 1 : -1;
            hypergeos += sign * mu2*mu2 
                * std::pow(static_cast<builtin_class>(x), a)
                * Hypergeometric3F2_Reg({{0.5+a, -b, a}}, {{1.5+a+c, 2.0+a}}, x);
        }
    }

    if (std::isnan(static_cast<builtin_class>(hypergeos))) {
        std::cerr << "Error: InteractionWindow_Greater(" << a << ", " << b 
            << ", " << c << ", " << mu1_ab << ", " << mu2_ab 
            << ") returns a NaN." << std::endl;
    }
    return overall * hypergeos;
}

coeff_class InteractionWindow_Greater_Zero(builtin_class b, 
        const builtin_class c, const std::array<builtin_class,2>& mu1_ab,
        const std::array<builtin_class,2>& mu2_ab) {
    coeff_class overall = M_PI/8;
    overall *= std::tgamma(1.0 + c) / std::tgamma(1.5 + b + c);
    // *= (-1)^1+b
    if (static_cast<int>(std::round(2*b)) % 2 == 1) {
        return 0; // pure imaginary so it's going to have to cancel
    } else if (static_cast<int>(std::round(b)) % 2 == 0) {
        overall = -overall;
    }

    // *= gamma(-1-b)/cos(\pi b) diverges when 2b is an integer (i.e. always)
    // but probably cancels against a zero in the sum of the hypergeometrics. 
    // We can't properly represent this right now, so here's a hack: (FIXME)
    b += 0.05;
    overall *= std::tgamma(-1-b) / std::cos(M_PI*b);
    
    coeff_class hypergeos = 0;
    for (std::size_t i = 0; i < 2; ++i) {
        coeff_class mu1 = mu1_ab[i];
        for (std::size_t j = 0; j < 2; ++j) {
            coeff_class mu2 = mu2_ab[j];
            coeff_class x = (mu1*mu1)/(mu2*mu2);
            int sign = (i+j)%2 == 0 ? 1 : -1;
            hypergeos += sign * mu2*mu2 
                * std::pow(static_cast<builtin_class>(x), -b)
                * Hypergeometric3F2_Reg(
                        {{-1.0-b, -b, -0.5-b-c}}, {{0.5-b, 1.0-b}}, x);
        }
    }

    // FIXME?? There's also another term that looks purely imaginary so I think
    // it should have to cancel out, right??

    return overall * hypergeos;
}

// coeff_class NPlus2Window(const char n, const char r,
        // const std::array<builtin_class,2>& mu1_ab,
        // const std::array<builtin_class,2>& mu2_ab) {
    // builtin_class a = r/2.0;
// 
    // if (mu1_ab[0] < mu2_ab[0]) {
        // return NPlus2Window_Less(n, a, mu1_ab, mu2_ab);
    // } else if (mu1_ab[0] == mu2_ab[0]) {
        // return NPlus2Window_Equal(n, a, mu1_ab);
    // } else {
        // return NPlus2Window_Greater(n, a, mu1_ab, mu2_ab);
    // }
// }

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

            coeff_class term = sign * std::pow(mu1, 0.25) * std::pow(mu2, 1.25)
                * std::pow(x, 0.25*n);
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

const DMatrix& MuPart(const std::array<char,2>& nr, const std::size_t partitions) {
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
