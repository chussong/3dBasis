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
// this works for every matrix computation except the same-n interaction,
// because that one passes different arguments
DMatrix MuPart(const std::size_t partitions, const coeff_class partWidth, 
        const MATRIX_TYPE type) {
    if (type == MAT_INNER || type == MAT_MASS) {
        return DMatrix::Identity(partitions, partitions);
    } else if (type == MAT_KINETIC) {
        return MuPart_Kinetic(partitions, partWidth);
    } else {
        std::cerr << "Error: the requested MuPart type has not yet been "
            << "implemented." << std::endl;
        throw std::logic_error("MuPart: type not implemented");
    }
}

DMatrix MuPart_Kinetic(const std::size_t partitions, 
        const coeff_class partWidth) {
    DMatrix output = DMatrix::Zero(partitions, partitions);
    for (std::size_t k = 0; k < partitions; ++k) {
        output(k, k) = partWidth*partWidth/2;
        output(k, k) *= k*k + (k+1)*(k+1);
    }
    return output;
}

// interaction (same n) matrix computations -----------------------------------

namespace {
    std::unordered_map<std::array<char,3>, DMatrix, 
        boost::hash<std::array<char,3>> > intCache;
    std::unordered_map<std::size_t,DMatrix> zeroMatrix;
    std::unordered_map<std::array<builtin_class,6>,coeff_class,
        boost::hash<std::array<builtin_class,6>> > hgfrCache;
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
        const std::size_t partitions, const coeff_class partWidth) {
    if (r[0]%2 == 1) {
        if (zeroMatrix.count(partitions) == 0) {
            zeroMatrix.emplace(partitions, DMatrix::Zero(partitions, partitions));
        }
        return zeroMatrix[partitions];
    }
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
    
    coeff_class groupB = 0;
    groupB -= mu_ab[0]*mu_ab[0]*Hypergeometric3F2_Reg(
            {{a, 0.5+a, -b}}, {{2.0+a, 1.5+a+c}}, 1);
    groupB += a*(mu_ab[1]-mu_ab[0])*(mu_ab[1]+mu_ab[0])*Hypergeometric3F2_Reg(
            {{0.5+a, 1.0+a, -b}}, {{2.0+a, 1.5+a+c}}, 1);
    groupB += mu_ab[0]*mu_ab[0]*std::pow(static_cast<builtin_class>(x), a)
        * Hypergeometric3F2_Reg({{a, 0.5+a, -b}}, {{2.0+a, 1.5+a+c}}, x);
    groupB *= std::tgamma(static_cast<builtin_class>(a));
    groupB *= std::tgamma(static_cast<builtin_class>(1.0 + c));

    return overall * (groupA + groupB);
}

coeff_class InteractionWindow_Greater(const builtin_class a, const builtin_class b, 
        const builtin_class c, const std::array<builtin_class,2>& mu1_ab,
        const std::array<builtin_class,2>& mu2_ab) {
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

    return overall * hypergeos;
}

// this 3F2 is based on GSL's series approximation for the 2F1
coeff_class Hypergeometric3F2(const std::array<builtin_class,3>& a, 
        const std::array<builtin_class,2>& b, const builtin_class x) {
    coeff_class sum_pos = 1.0;
    coeff_class sum_neg = 0.0;
    coeff_class del_pos = 1.0;
    coeff_class del_neg = 0.0;
    coeff_class del = 1.0;
    coeff_class del_prev;
    coeff_class k = 0.0;
    int i = 0;

    coeff_class result;

    do {
        if(++i > 30000) {
            result = sum_pos - sum_neg;
            GSL_ERROR ("didn't converge", GSL_EMAXITER);
        }
        del_prev = del;
        del *= (a[0]+k)*(a[1]+k)*(a[2]+k) * x / ((b[0]+k)*(b[1]+k) * (k+1.0));

        if(del > 0.0) {
            del_pos  =  del;
            sum_pos +=  del;
        }
        else if(del == 0.0) {
            /* Exact termination (a, b, or c was a negative integer).
            */
            del_pos = 0.0;
            del_neg = 0.0;
            break;
        }
        else {
            del_neg  = -del;
            sum_neg -=  del;
        }

        /*
         * This stopping criteria is taken from the thesis
         * "Computation of Hypergeometic Functions" by J. Pearson, pg. 31
         * (http://people.maths.ox.ac.uk/porterm/research/pearson_final.pdf)
         * and fixes bug #45926
         */
        if (std::abs<builtin_class>(del_prev / (sum_pos - sum_neg)) < EPSILON &&
            std::abs<builtin_class>(del / (sum_pos - sum_neg)) < EPSILON)
            break;

        k += 1.0;
    } while(std::abs<builtin_class>((del_pos + del_neg)/(sum_pos-sum_neg)) > EPSILON);

    result = sum_pos - sum_neg;

    return result;
}

coeff_class Hypergeometric3F2_Reg(const std::array<builtin_class,3>& a, 
        const std::array<builtin_class,2>& b, const builtin_class x) {
    std::array<builtin_class,6> params = {{a[0], a[1], a[2], b[0], b[1], x}};
    if (hgfrCache.count(params) == 0) {
        coeff_class reg = std::tgamma(static_cast<builtin_class>(b[0]))
            * std::tgamma(static_cast<builtin_class>(b[1]));
        hgfrCache.emplace(params, Hypergeometric3F2(a, b, x) / reg);
    }

    return hgfrCache[params];
}

const DMatrix& MuPart(const char r, const std::size_t partitions,
        const coeff_class partWidth) {
    throw std::logic_error("MuPart for n+2: this function not yet implemented");
}
