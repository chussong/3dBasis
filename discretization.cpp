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
DMatrix MuPart(const Mono& A, const Mono& B, const std::size_t partitions, 
        const coeff_class partWidth, const MATRIX_TYPE type) {
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
        output(k, k) = partWidth*partWidth;
        output(k, k) *= (k*k + (k+1)*(k+1))/2;
    }
    return output;
}

// interaction (same n) matrix computations -----------------------------------

void InteractionCache::SetPartitions(const std::size_t partitions, 
        const coeff_class partWidth) {
    this->partitions = partitions;
    this->partWidth = partWidth;
}

bool InteractionCache::HasPartitions(const std::size_t partitions,
        const coeff_class partWidth) {
    return (partitions == this->partitions && partWidth == this->partWidth);
}

void InteractionCache::Emplace(const std::array<char,3>& key, 
        DMatrix value) {
    cache.emplace(key, std::move(value));
}

bool InteractionCache::Contains(const std::array<char,3>& key) {
    return (cache.count(key) != 0);
}

const DMatrix& InteractionCache::operator[](const std::array<char,3>& key)
    const {
    return cache.at(key);
}

void InteractionCache::Clear() {
    cache.clear();
}

namespace {
    InteractionCache intCache;
    std::unordered_map<std::size_t,DMatrix> zeroMatrix;
} // anonymous namespace

// this is the MuPart for the same-n interaction matrix
//
// Interestingly, there appears to be an analytic answer to this: (-1)^a \pi
// ( (2 b! pFq({0.5, 0.5+a, -c}, {1.5, 1.5+a+b}, 1)) +
// (c! Gamma[a] Gamma[1.5+a+b] pFq_reg({a, 0.5+a, -b}, {1+a, 1.5+a+c}, 1)) )
// / (4 Gamma[0.5-a] Gamma[1.5+a+b])
// NOTE: the a,b,c given above are the exponents of r^2, 1-r^2, 1-(ar)^2,
// so they're half of the entries of the r array
//
// the analytic answer has a couple of interesting properties: first, for 
// half-integer a, it's 0 due to the gamma function in the denominator, but most
// importantly both of those pFq functions are just finite polynomials due to
// the negative integers in the p lists
//
// this function returns the block of the interaction matrix corresponding to
// one pair of monomials. The sum of all of these blocks will be the block's
// actual value.
const DMatrix& MuPart(const std::array<char,3>& r, 
        const std::size_t partitions, const coeff_class partWidth) {
    if (r[0]%2 == 1) {
        if (zeroMatrix.count(partitions) == 0) {
            zeroMatrix.emplace(partitions, DMatrix::Zero(partitions, partitions));
        }
        return zeroMatrix[partitions];
    }
    if (!intCache.Contains(r)) {
        // 2d trapezoid rule to approximate the integral's values
        //
        // FIXME: trapezoid rule only gives ~2 digits of precision. Use
        // analytic answer.
        DMatrix cornerVals = DMatrix::Identity(partitions+1, partitions+1);
        cornerVals *= RIntegral(r[0], r[1], r[2], 1);
        for (std::size_t winA = 0; winA <= partitions; ++winA) {
            for (std::size_t winB = 0; winB <= partitions; ++winB) {
                coeff_class alpha = (winA*partWidth) / (winB*partWidth);
                cornerVals(winA, winB) = RIntegral(r[0], r[1], r[2], alpha);
                // Below is the contribution from (winA, winB) <-> (winB, winA);
                // this relation appears to be valid in general, irrespective of
                // the min(1/a, 1) in the upper bound of the integral. That is,
                // I(a, b, c, 1/alpha) = alpha^(2a+1) I(a, c, b, alpha)
                cornerVals(winB, winA) = 
                    std::pow<builtin_class>(alpha, r[0]+1)
                    * RIntegral(r[0], r[2], r[1], alpha);
            }
        }
        DMatrix block(partitions, partitions);
        for (std::size_t winA = 0; winA < partitions; ++winA) {
            for (std::size_t winB = 0; winB < partitions; ++winB) {
                block(winA, winB) = partWidth*partWidth
                    * cornerVals.block(winA, winB, 2, 2).mean();
            }
        }
        intCache.Emplace(r, std::move(block));
    }

    return intCache[r];
}

// FIXME: the 3F2s from here can be reused for 4 windows each
//
// FIXME: this is only the correct answer for mu1 <= mu2
coeff_class InteractionWindow(const std::array<char,3>& r,
        const std::array<coeff_class,2>& mu1_ab,
        const std::array<coeff_class,2>& mu2_ab) {
    coeff_class a = r[0]/2.0;
    coeff_class b = r[1]/2.0;
    coeff_class c = r[2]/2.0;

    coeff_class output = std::tgamma(static_cast<builtin_class>(1.0 + b));
    output /= std::tgamma(static_cast<builtin_class>(0.5 - a));
    output *= std::pow(M_PI, 1.5) / 4;
    if (r[0]%4 == 2) output = -output;

    coeff_class hypergeos = 0;
    for (coeff_class mu1 : mu1_ab) {
        for (std::size_t i = 0; i < 2; ++i) {
            coeff_class mu2 = mu2_ab[i];
            int sign = i%2 == 0 ? -1 : 1;
            hypergeos += sign * mu1*mu2 * Hypergeometric3F2_Reg(
                    {{0.5+a, -c, -0.5}}, {{1.5+a+b, 1.5}}, (mu1*mu1)/(mu2*mu2));
        }
    }

    return output * hypergeos;
}

// r integral if \alpha < 1; {a, b, c} are the respective exponents of 
// {r^2, 1 - r^2, 1 - \alpha^2 r^2}
//
// WARNING: IF ALPHA > 1, THIS IS WRONG; THIS MAY LEAD TO TROUBLE WITH THE UPPER
// BOUND OF THE INTEGRAL BEING MIN(1/ALPHA, 1)
//
// FIXME: IT SEEMS THERE IS ACTUALLY AN ANALYTIC ANSWER FOR THE INTEGRAL OF THIS
// BETWEEN x1 AND x2; IT'S A COUPLE OF pFq FUNCTIONS WHICH ARE A FINITE-DEGREE
// POLYNOMIAL. SHOULD USE THAT INSTEAD OF AVERAGING THIS
coeff_class RIntegral(const char twoA, const char twoB, const char twoC, 
        const coeff_class alpha) {
    if (twoA%2 == 1) return 0;

    coeff_class a = static_cast<coeff_class>(twoA)/2;
    coeff_class b = static_cast<coeff_class>(twoB)/2;
    coeff_class c = static_cast<coeff_class>(twoC)/2;

    coeff_class logOutput = std::lgamma(static_cast<builtin_class>(1.0 + b));
    logOutput -= std::lgamma(static_cast<builtin_class>(0.5 - a));
    logOutput -= std::lgamma(static_cast<builtin_class>(1.5 + a + b));
    coeff_class output = M_PI*std::exp<builtin_class>(logOutput).real()/2;
    if (twoA%4 == 2) output = -output;

    // gsl_sf_result hyper;
    // int errorCode = gsl_sf_hyperg_2F1_e(-b, -0.5-a-b-c, 0.5-a-b, alpha*alpha, &hyper);
    // if (errorCode != 0) {
        // std::cerr << "Error in RIntegral->gsl_sf_hyperg_2F1\n";
        // errorCode = gsl_sf_hyperg_2F1_renorm_e(-b, -0.5-a-b-c, 0.5-a-b, alpha*alpha, &hyper);
        // if (errorCode != 0) {
            // std::cerr << "Renormalized one returned error too: " 
                // << gsl_strerror(errorCode) << std::endl;
        // }
    // }
    // return output*(hyper.val);

    return output*Hypergeometric2F1(0.5 + a, -c, 1.5 + a + b, alpha*alpha);
}

// ------------------------- TEMPORARY, JUST FOR TESTING ----------------------
// this is hyperg_2F1_series from GSL, copied here because it has static linkage
coeff_class Hypergeometric2F1(const coeff_class a, const coeff_class b,
        const coeff_class c, const coeff_class x) {
    coeff_class sum_pos = 1.0;
    coeff_class sum_neg = 0.0;
    coeff_class del_pos = 1.0;
    coeff_class del_neg = 0.0;
    coeff_class del = 1.0;
    coeff_class del_prev;
    coeff_class k = 0.0;
    int i = 0;

    coeff_class result;

    if(std::abs<builtin_class>(c) < EPSILON) {
        GSL_ERROR ("c singularity", GSL_EDOM);
    }

    do {
        if(++i > 30000) {
            result = sum_pos - sum_neg;
            GSL_ERROR ("didn't converge", GSL_EMAXITER);
        }
        del_prev = del;
        del *= (a+k)*(b+k) * x / ((c+k) * (k+1.0));  /* Gauss series */

        if(del > 0.0) {
            del_pos  =  del;
            sum_pos +=  del;
        }
        else if(del == 0.0) {
            /* Exact termination (a or b was a negative integer).
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

// this 3F2 is based on GSL's series approximation for the 2F1
//
// FIXME?? in this problem, a[2] and b[1] are always -1/2 and 3/2, respectively,
// so we don't actually need to be passing them as arguments
coeff_class Hypergeometric3F2(const std::array<coeff_class,3>& a, 
        const std::array<coeff_class,2>& b, const coeff_class x) {
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

// FIXME?? b[1] is always 3/2 in this problem, so we can simplify this a bit
coeff_class Hypergeometric3F2_Reg(const std::array<coeff_class,3>& a, 
        const std::array<coeff_class,2>& b, const coeff_class x) {
    coeff_class reg = std::tgamma(static_cast<builtin_class>(b[0]))
        * std::tgamma(static_cast<builtin_class>(b[1]));
    return Hypergeometric3F2(a, b, x) / reg;
}
