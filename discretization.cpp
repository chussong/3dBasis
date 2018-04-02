#include "discretization.hpp"

// return exponent P of \mu^P to deal with when multiplying A and B; this is
// N-3 + total number of perp derivates + 1 because the integrals are originally
// over \mu^2 instead of \mu (and I think only one of them survives so not +2??)
//
// I believe this corresponds to the exponent of \mu found in the inner product
char MuExponent(const Mono& A, const Mono& B) {
	return A.NParticles() - 3 + A.TotalPt() + B.TotalPt() + 1;
}

DVector MuIntegral(const Mono& A, const Mono& B, const std::size_t partitions,
        const coeff_class partWidth, const MATRIX_TYPE calculationType) {
    if (calculationType == MAT_INNER || calculationType == MAT_MASS) {
        DVector output = MuIntegral_Body(MuExponent(A, B), partitions, partWidth); 
        for (Eigen::Index k = 0; k < output.size(); ++k) {
            output(k) *= MuNorm(A, k, partWidth)*MuNorm(B, k, partWidth)/M_PI;
        }
        return output;
    } else {
        return DVector::Zero(partitions);
    }
}

// this only returns a DVector instead of a DMatrix because the mass term has a
// \delta_{ij} so only the diagonal is nonzero for each pair of monomials
//
// the 4 is from converting two integrals from d\mu^2 to d\mu
DVector MuIntegral_Body(const char muExp, const std::size_t partitions,
                const coeff_class partitionWidth) {
    DVector output(partitions);
    coeff_class left = 0;
    coeff_class right = 0;
    for (std::size_t i = 0; i < partitions; ++i) {
        left = right;
        right += partitionWidth;
        output[i] = (std::pow<double>(right, muExp+1) 
                - std::pow<double>(left, muExp+1)) / (muExp + 1);
        // std::cout << "Integral(" << (int)muExp << ", " << i << ", " 
            // << partitionWidth << ") = " << output[i] << std::endl;
    }
    return output;
}

// normalization for mu integrals; this is the script N_{n, \lambda, k} in
// Zuhair's notes
coeff_class MuNorm(const Mono& A, const std::size_t k, 
        const coeff_class partWidth) {
    // auto exponent = A.NParticles() - 1 + 2*A.TotalPt();
    auto exponent = MuExponent(A, A) + 1;
    coeff_class left = k*partWidth;
    coeff_class right = (k+1)*partWidth;

    // coeff_class denom = std::pow<builtin_class>(right, exponent);
    // denom -= std::pow<builtin_class>(left, exponent);
    // denom /= M_PI*exponent;

    coeff_class square = M_PI*exponent/(std::pow<builtin_class>(right, exponent)
            - std::pow<builtin_class>(left, exponent));

    auto output = std::sqrt<builtin_class>(square);
    if (output.imag() != 0) {
        std::cerr << "Warning: MuNorm(" << A << ", " << k << ", " << partWidth
            << ") is complex!?" << std::endl;
    }
    // std::cout << "MuNorm(" << A << ", " << k << ", " << partWidth << ") = "
        // << output.real() << std::endl;
    return output.real()/2;
}

// expands a mass matrix with the mu factors elided into one with them partially
// integrated over, by taking each entry of the input matrix and turning it into
// a block. Each block is diagonal because this is a mass matrix; it's not quite
// so simple for the interaction term
DMatrix PartitionMu_Mass(const Basis<Mono>& minimalBasis, const DMatrix& mass,
		const std::size_t partitions, const coeff_class partWidth) {
    DMatrix output = DMatrix::Zero(mass.rows()*partitions, mass.cols()*partitions);
    for (Eigen::Index row = 0; row < mass.rows(); ++row) {
        for (Eigen::Index col = 0; col < mass.cols(); ++col) {
            DVector integrated = MuIntegral(minimalBasis[row], 
                            minimalBasis[col], partitions, partWidth,
                            MAT_MASS);
            for (std::size_t p = 0; p < partitions; ++p) {
                    output(row*partitions + p, col*partitions + p) 
                            = mass(row, col) * integrated(p);
            }
        }
    }
    return output;
}

// for each monomial in the basis, create a vector with 1/sqrt(IPintegral) over
// each mu window. Return these as a matrix; the full basis state will then be
// [\sum_b X_{ib} Y_{bp} M_b], where X is the polysOnMinBasis matrix, Y is the
// output of this function, and M is the minBasis
DMatrix DiscretizeMonos(const Basis<Mono>& minBasis, 
        const std::size_t partitions, const coeff_class partWidth) {
    DMatrix output(minBasis.size(), partitions);
    /**************************************************************************/
    /* DANGER! PRETTY SURE THIS SHOULD BE DELETED AND IS NOT USEFUL AT ALL!!! */
    /**************************************************************************/
    for (std::size_t i = 0; i < minBasis.size(); ++i) {
        DVector IP = MuIntegral(
                    minBasis[i], minBasis[i], partitions, partWidth, MAT_INNER);
        for (std::size_t p = 0; p < partitions; ++p) {
            auto root = std::sqrt<builtin_class>(IP(p));
            if (root.imag() != 0) {
                throw std::logic_error("DiscretizeMonos: negative IP given.");
            }
            output(i, p) = 1/root.real();
        }
    }
    return output;
}

// given monomials A and B, give the PxP tensor representing block A,B of the
// specified matrix
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

// return the mu integral, summed over all partitions, for each pair of
// monomials in the given basis
DMatrix MuTotal(const Basis<Mono>& minBasis, const std::size_t partitions,
        const coeff_class partWidth, const MATRIX_TYPE calculationType) {
    DMatrix output(minBasis.size(), minBasis.size());
    for (Eigen::Index row = 0; row < output.rows(); ++row) {
        DVector windows = MuIntegral(minBasis[row], minBasis[row], partitions,
                partWidth, calculationType);
        output(row, row) = windows.sum();
        for (Eigen::Index col = row+1; col < output.cols(); ++col) {
            DVector windows = MuIntegral(minBasis[row], minBasis[col], 
                    partitions, partWidth, calculationType);
            output(row, col) = windows.sum();
            output(col, row) = output(row, col);
        }
    }
    return output;
}

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
} // anonymous namespace

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
DMatrix InteractionMu(const std::array<char,3> r, 
        const std::size_t partitions, const coeff_class partWidth) {
    if (!intCache.Contains(r)) {
        // for trapezoid rule, need values at exterior points without "*2",
        // meaning, values at alpha = 0 and alpha = inf
        DMatrix block(partitions, partitions);
        for (std::size_t winA = 0; winA <= partitions; ++winA) {
            block(winA, winA) = partWidth*RIntegral(r[0], r[1], r[2], 1);
            for (std::size_t winB = winA+1; winB < partitions; ++winB) {
                coeff_class alpha = (winA*partWidth) / (winB*partWidth);

                block(winA, winB) = partWidth*RIntegral(r[0], r[1], r[2], alpha);
                // Below is the contribution from (winA, winB) <-> (winB, winA);
                // this relation appears to be valid in general, irrespective of
                // the min(1/a, 1) in the upper bound of the integral. That is,
                // I(a, b, c, 1/alpha) = alpha^(2a+1) I(a, c, b, alpha)
                block(winB, winA) = partWidth
                    * std::pow<builtin_class>(alpha, 2*r[0]+1)
                    * RIntegral(r[0], r[2], r[1], alpha);
            }
        }
        intCache.Emplace(r, block);
    }

    return intCache[r];
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
    coeff_class a = static_cast<coeff_class>(twoA)/2;
    coeff_class b = static_cast<coeff_class>(twoB)/2;
    coeff_class c = static_cast<coeff_class>(twoC)/2;
    coeff_class logOutput = (2*c)*std::log<builtin_class>(alpha).real();
    logOutput += std::lgamma(static_cast<builtin_class>(0.5 + a + c));
    logOutput += std::lgamma(static_cast<builtin_class>(1.0 + b));
    logOutput -= std::lgamma(static_cast<builtin_class>(1.5 + a + b + c));
    coeff_class output = std::exp<builtin_class>(logOutput).real()/2;
    if (twoC%2 == 1) output = -output;

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

    return output*Hypergeometric2F1(-c, -0.5-a-b-c, 0.5-a-c, 1/(alpha*alpha));
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
