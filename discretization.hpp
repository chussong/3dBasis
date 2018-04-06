#ifndef DISCRETIZATION_HPP
#define DISCRETIZATION_HPP

#include <array>
#include <vector>
#include <cmath>
#include <unordered_map>

#include <boost/functional/hash.hpp>

#include "constants.hpp"
#include "mono.hpp"
#include "basis.hpp"

DMatrix DiscretizePolys(const DMatrix& polysOnMinBasis, 
        const std::size_t partitions);

// direct matrices ------------------------------------------------------------

DMatrix MuPart(const std::size_t partitions, const coeff_class partWidth, 
        const MATRIX_TYPE type);
DMatrix MuPart_Kinetic(const std::size_t partitions, 
        const coeff_class partWidth);

// same-n interactions --------------------------------------------------------

const DMatrix& MuPart(const std::array<char,3>& r, 
        const std::size_t partitions, const coeff_class partWidth);
coeff_class InteractionWindow(const std::array<char,3>& r,
        const std::array<builtin_class,2>& mu1_ab,
        const std::array<builtin_class,2>& mu2_ab);
coeff_class InteractionWindow_Less(const builtin_class a, const builtin_class b,
        const builtin_class c, const std::array<builtin_class,2>& mu1_ab,
        const std::array<builtin_class,2>& mu2_ab);
coeff_class InteractionWindow_Equal(const builtin_class a, const builtin_class b,
        const builtin_class c, const std::array<builtin_class,2>& mu_ab);
coeff_class InteractionWindow_Greater(const builtin_class a, const builtin_class b,
        const builtin_class c, const std::array<builtin_class,2>& mu1_ab,
        const std::array<builtin_class,2>& mu2_ab);
coeff_class Hypergeometric3F2(const std::array<builtin_class,3>& a, 
        const std::array<builtin_class,2>& b, const builtin_class x);
coeff_class Hypergeometric3F2_Reg(const std::array<builtin_class,3>& a, 
        const std::array<builtin_class,2>& b, const builtin_class x);

// n+2 interactions -----------------------------------------------------------

const DMatrix& MuPart(const char r, const std::size_t partitions,
        const coeff_class partWidth);

// generic hypergeometric function --------------------------------------------

constexpr int ITERATION_LIMIT = 30000;
constexpr coeff_class PRECISION_LIMIT = 1e-10;

// this is an adaptation of the series representation of GSL's Hypergeometric2F1
//
// FIXME: this function should just call the one below it and de-regularize the
// result by multiplying by Prod[Gamma[b]]
template<int P, int Q>
coeff_class HypergeometricPFQ(const std::array<builtin_class,P>& a, 
        const std::array<builtin_class,Q>& b, const builtin_class x) {
    coeff_class sum_pos = 1.0;
    coeff_class sum_neg = 0.0;
    coeff_class del_pos = 1.0;
    coeff_class del_neg = 0.0;
    coeff_class del = 1.0;
    coeff_class del_prev;
    coeff_class k = 0.0;
    int i = 0;

    do {
        if(++i > ITERATION_LIMIT) {
            // return sum_pos - sum_neg;
            throw (std::runtime_error("HypergeometricPFQ did not converge."));
        }
        del_prev = del;
        for (coeff_class a_i : a) del *= (a_i + k);
        for (coeff_class b_i : b) del /= (b_i + k);
        del *= x / (k + 1.0);

        if(del > 0.0) {
            del_pos  =  del;
            sum_pos +=  del;
        }
        else if(del == 0.0) {
            /* Exact termination (some a[i] was a negative integer).
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
        if (std::abs<builtin_class>(del_prev / (sum_pos - sum_neg)) 
                    < PRECISION_LIMIT 
                && std::abs<builtin_class>(del / (sum_pos - sum_neg)) 
                    < PRECISION_LIMIT) {
            break;
        }

        k += 1.0;
    } while(std::abs<builtin_class>((del_pos + del_neg)/(sum_pos-sum_neg)) 
            > PRECISION_LIMIT);

    return sum_pos - sum_neg;
}

template<int P, int Q>
coeff_class HypergeometricPFQ_Reg(const std::array<builtin_class,P>& a, 
        const std::array<builtin_class,Q>& b, const builtin_class x) {
    coeff_class del = 1.0;
    coeff_class del_prev;
    builtin_class k = 0.0; // k is the index of the most recent COMPLETED term
    int i = 0;

    std::set<coeff_class> zeros;
    for (int i = 0; i < P; ++i) {
        if (a[i] < 0 && std::abs(std::round(a[i]) - a[i]) < EPSILON) {
            zeros.insert(std::round(-a[i]));
        }
    }
    std::set<coeff_class> divergences;
    for (int i = 0; i < Q; ++i) {
        if (b[i] < 0 && std::abs(std::round(b[i]) - b[i]) < EPSILON) {
            divergences.insert(std::round(-b[i]));
        }
    }

    if (divergences.size() > 0) {
        // if there is a zero before the final divergence, all terms will be 0
        if (zeros.size() > 0 && *zeros.begin() <= *divergences.rbegin()) {
            return 0;
        } else {
            // first nonzero term of regularized series
            k = *divergences.rbegin() + 1;
            del = std::pow(x, k) / std::tgamma(k+1);
            for (builtin_class a_i : a) {
                if (zeros.count(a_i) == 0) {
                    del *= std::tgamma(a_i + k) / std::tgamma(a_i);
                } else {
                    coeff_class prod = 1;
                    for (builtin_class n = 1; n <= std::round(a_i + k); ++n) {
                        prod *= a_i;
                    }
                }
            }
            for (builtin_class b_i : b) {
                if (divergences.count(b_i) == 0) del /= std::tgamma(b_i + k);
            }
        }
    } else {
        for (builtin_class b_i : b) del /= std::tgamma(b_i);
    }

    coeff_class sum_pos = 0.0;
    coeff_class sum_neg = 0.0;
    coeff_class del_pos = 0.0;
    coeff_class del_neg = 0.0;

    if (del >= 0.0) {
        sum_pos = del;
        del_pos = del;
    } else {
        sum_neg = del;
        del_neg = del;
    }

    do {
        if(++i > ITERATION_LIMIT) {
            // return sum_pos - sum_neg;
            throw (std::runtime_error("HypergeometricPFQ_Reg did not converge."));
        }
        del_prev = del;
        for (coeff_class a_i : a) del *= (a_i + k);
        for (coeff_class b_i : b) del /= (b_i + k);
        del *= x / (k + 1.0);

        if(del > 0.0) {
            del_pos  =  del;
            sum_pos +=  del;
        }
        else if(del == 0.0) {
            /* Exact termination (some a[i] was a negative integer).
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
        if (std::abs<builtin_class>(del_prev / (sum_pos - sum_neg)) 
                    < PRECISION_LIMIT 
                && std::abs<builtin_class>(del / (sum_pos - sum_neg)) 
                    < PRECISION_LIMIT) {
            break;
        }

        k += 1.0;
    } while(std::abs<builtin_class>((del_pos + del_neg)/(sum_pos-sum_neg)) 
            > PRECISION_LIMIT);

    return sum_pos - sum_neg;
}

#endif
