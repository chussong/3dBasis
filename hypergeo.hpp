#ifndef HYPERGEO_HPP
#define HYPERGEO_HPP

#include <cmath>
#include <array>
#include <set>
#include <unordered_map>
#include <iostream>

#include <boost/functional/hash.hpp>

#include "constants.hpp"
#include "io.hpp"

// generalized hypergeometric function via templates. Also includes memoized
// versions of 2F1, 3F2, and 4F3 in cpp file, but the templates are not memoized
// and work with the header by itself
//
// warning: my application uses small numbers (< 10 or so) for a and b, so I 
// just use std::tgamma; if you want this to be usable for larger arguments,
// change std::tgamma to std::lgamma and adjust arithmetic accordingly

// parameters -----------------------------------------------------------------

constexpr int ITERATION_LIMIT = 1e8;
constexpr coeff_class PRECISION_LIMIT = 1e-10;

// memoized functions for specific cases --------------------------------------

coeff_class Hypergeometric2F1(const builtin_class a, const builtin_class b,
                              const builtin_class c, const builtin_class x);
coeff_class Hypergeometric2F1_Reg(const builtin_class a, const builtin_class b,
                                  const builtin_class c, const builtin_class x);
coeff_class Hypergeometric3F2(const builtin_class a1, const builtin_class a2,
                              const builtin_class a3, const builtin_class b1,
                              const builtin_class b2, const builtin_class x);
coeff_class Hypergeometric3F2(const std::array<builtin_class,3>& a, 
                              const std::array<builtin_class,2>& b, 
                              const builtin_class x);
coeff_class Hypergeometric3F2(const std::array<builtin_class,6>& params);
coeff_class Hypergeometric3F2_Reg(const builtin_class a1, 
                                  const builtin_class a2,
                                  const builtin_class a3,
                                  const builtin_class b1,
                                  const builtin_class b2,
                                  const builtin_class x);
coeff_class Hypergeometric3F2_Reg(const std::array<builtin_class,3>& a, 
        const std::array<builtin_class,2>& b, const builtin_class x);
coeff_class Hypergeometric3F2_Reg(const std::array<builtin_class,6>& params);
coeff_class Hypergeometric4F3(const builtin_class a1, const builtin_class a2,
                              const builtin_class a3, const builtin_class a4,
                              const builtin_class b1, const builtin_class b2,
                              const builtin_class b3, const builtin_class x);
coeff_class Hypergeometric4F3(const std::array<builtin_class,4>& a, 
                              const std::array<builtin_class,3>& b, 
                              const builtin_class x);
coeff_class Hypergeometric4F3(const std::array<builtin_class,8>& params);
coeff_class Hypergeometric4F3_Reg(const builtin_class a1, const builtin_class a2,
                                  const builtin_class a3, const builtin_class a4,
                                  const builtin_class b1, const builtin_class b2,
                                  const builtin_class b3, const builtin_class x);
coeff_class Hypergeometric4F3_Reg(const std::array<builtin_class,4>& a, 
                                  const std::array<builtin_class,3>& b, 
                                  const builtin_class x);
coeff_class Hypergeometric4F3_Reg(const std::array<builtin_class,8>& params);

// general templates ----------------------------------------------------------

inline bool IsNegInt(const builtin_class x) {
    return std::round(x) < 0 && std::abs(std::round(x) - x) < EPSILON;
}

// this is an adaptation of the series representation of GSL's Hypergeometric2F1
//
// CAREFUL: this is the RENORMALIZED generalized hypergeometric function, i.e.
// it has been divided by Product(Gamma(b[i])); the ordinary PFQ is implemented 
// in terms of the renormalized one, although it would work very slightly better
// to have it as a separate, nearly identical function
template<std::size_t P, std::size_t Q>
coeff_class HypergeometricPFQ_Body(const std::array<builtin_class,P>& a, 
        const std::array<builtin_class,Q>& b, const builtin_class x) {

    coeff_class del = 1.0;
    coeff_class del_prev;
    builtin_class k = 0.0; // k is the index of the most recent COMPLETED term
    int i = 0;

    // if any a is a negative integer, all terms AFTER -a[i] will be 0
    std::set<coeff_class> zeros;
    for (std::size_t i = 0; i < P; ++i) {
        if (IsNegInt(a[i])) {
            zeros.insert(std::round(-a[i]));
        }
    }
    // if any b is a negative integer, all terms BEFORE -b[i] will be 0
    // note: these are not actually divergences because PFQ is renormalized
    std::set<coeff_class> divergences;
    for (std::size_t i = 0; i < Q; ++i) {
        if (IsNegInt(b[i])) {
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

template<std::size_t P, std::size_t Q>
coeff_class HypergeoUnitArgument_Reg(const std::array<builtin_class,P>& a, 
                                     const std::array<builtin_class,Q>& b) {
    // for lack of better ideas, just try the series and hope it converges
    return HypergeometricPFQ_Body<P,Q>(a, b, 1);
}

template<>
inline coeff_class HypergeoUnitArgument_Reg(const std::array<builtin_class,2>& a, 
        const std::array<builtin_class,1>& b) {
    // builtin_class num = std::tgamma(b[0] - a[0] - a[1]);
    // builtin_class den = std::tgamma(b[0] - a[0]) * std::tgamma(b[0] - a[1]);
    // return num / den;

    builtin_class num = std::lgamma(b[0] - a[0] - a[1]);
    builtin_class den = std::lgamma(b[0] - a[0]) + std::lgamma(b[0] - a[1]);
    return std::exp(num - den);
}

template<>
inline coeff_class HypergeoUnitArgument_Reg(
        const std::array<builtin_class,3>& aRef, 
        const std::array<builtin_class,2>& bRef) {

    // try the series, see if it converges; if not, try a transformation
    coeff_class output = std::nan("");
    try {
        output = HypergeometricPFQ_Body<3,2>(aRef, bRef, 1);
    }
    catch (const std::runtime_error& err) {
        // try Thomae's transformation
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 2; ++j) {
                if (IsNegInt(bRef[j] - aRef[i])) {
                    // std::cout << "THOMAE SUCCESSFUL" << std::endl;
                    auto a = aRef;
                    auto b = bRef;
                    if (i != 0) std::swap(a[i], a[0]);
                    if (j != 0) std::swap(b[j], b[0]);
                    std::array<builtin_class,3> a2;
                    a2[0] = b[0] - a[0];
                    a2[1] = b[0] - a[1];
                    a2[2] = b[0] + b[1] - a[0] - a[1] - a[2];
                    std::array<builtin_class,2> b2;
                    b2[0] = b[0] + b[1] - a[0] - a[2];
                    b2[0] = b[0] + b[1] - a[0] - a[1];
                    try {
                        output = std::tgamma(a2[2])
                                 * HypergeometricPFQ_Body<3,2>(a2, b2, 1);
                    }
                    catch (const std::runtime_error& err) {
                        output = std::nan("");
                    }
                }
            }
        }
    }

    if (std::isfinite(static_cast<builtin_class>(output))) {
        return output;
    } else {
        throw std::runtime_error("HypergeoUnitArgument_Reg<3,2> did not "
                                 "converge.");
    }
}

template<std::size_t P, std::size_t Q>
coeff_class HypergeoSpecialCase_Reg(const std::array<builtin_class,P>& a,
        const std::array<builtin_class,Q>& b, const builtin_class x) {
    for (std::size_t i = 0; i < P; ++i) {
        if (a[i] == 0) {
            // std::cout << P << 'F' << Q << '(' << a << "; " << b << "; " << x
                // << ") has a 0 in." << std::endl;
            coeff_class gammaProd = 1;
            for (auto b_i : b) {
                if (std::abs(b_i) < EPSILON || IsNegInt(b_i)) {
                    gammaProd = 0;
                    break;
                } else {
                    gammaProd *= std::tgamma(b_i);
                }
            }
            if (gammaProd == 0) {
                return 0;
            } else {
                return 1/gammaProd;
            }
        }

        // if any of the a are negative integers, don't bother with special case
        if (IsNegInt(a[i])) {
            return std::nan("");
        }

        // check the b to see if any match this a (and cancel them if so)
        for (std::size_t j = 0; j < Q; ++j) {
            if (a[i] == b[j]) {
                // this declaration is necessary so that it doesn't try to 
                // instantiate the template with P or Q < 0, even though this 
                // code path is inaccessible if P or Q is 0
                constexpr std::size_t lowerP = P > 0 ? P-1 : P;
                std::array<builtin_class,lowerP> a2;
                for (std::size_t k = 0; k < P; ++k) {
                    if (k < i) {
                        a2[k] = a[k];
                    } else if (k > i) {
                        a2[k-1] = a[k];
                    }
                }

                constexpr std::size_t lowerQ = Q > 0 ? Q-1 : Q;
                std::array<builtin_class,lowerQ> b2;
                for (std::size_t k = 0; k < Q; ++k) {
                    if (k < j) {
                        b2[k] = b[k];
                    } else if (k > j) {
                        b2[k-1] = b[k];
                    }
                }

                coeff_class value = HypergeoSpecialCase_Reg<lowerP, 
                                                            lowerQ>(a2, b2, x);
                if (std::isnan(static_cast<builtin_class>(value))) {
                    // the reduced one isn't a special case, so do series
                    value = HypergeometricPFQ_Body<lowerP, lowerQ>(a2, b2, x);
                }

                return value / std::tgamma(b[j]);
            }
        }
    }

    if (std::abs(x - 1) < EPSILON) return HypergeoUnitArgument_Reg<P,Q>(a, b);

    // returning this nan causes the body to run (i.e. this case isn't special)
    return std::nan("");
}

template<std::size_t P, std::size_t Q>
coeff_class HypergeometricPFQ_Reg(const std::array<builtin_class,P>& a, 
        const std::array<builtin_class,Q>& b, const builtin_class x) {
    coeff_class specialCase = HypergeoSpecialCase_Reg<P,Q>(a, b, x);
    if (!std::isnan(static_cast<builtin_class>(specialCase))) return specialCase;

    return HypergeometricPFQ_Body<P,Q>(a, b, x);
}

template<>
inline coeff_class HypergeometricPFQ_Reg<2,1>(const std::array<builtin_class,2>& a,
        const std::array<builtin_class,1>& b, const builtin_class x) {
    // previously I called the gsl implementation of this, but as of May 2018 it
    // actually has a bug in gsl_sf_hyperg_2F1_renorm (renorm version only) 
    // where it gives the wrong answer if c (our b[0]) is negative, so we'll use
    // the implementation in this file

    // builtin_class gslResult = gsl_sf_hyperg_2F1_renorm(a[0], a[1], b[0], x);
    // if (std::isfinite(gslResult)) {
        // return gslResult;
    // } else {
        coeff_class specialCase = HypergeoSpecialCase_Reg<2,1>(a, b, x);
        if (!std::isnan(static_cast<builtin_class>(specialCase))) return specialCase;

        return HypergeometricPFQ_Body<2,1>(a, b, x);
    // }
}

template<std::size_t P, std::size_t Q>
coeff_class HypergeometricPFQ(const std::array<builtin_class,P>& a, 
        const std::array<builtin_class,Q>& b, const builtin_class x) {
    coeff_class reg =  HypergeometricPFQ_Reg<P,Q>(a, b, x);
    for (auto b_i : b) reg *= std::tgamma(b_i);
    return reg;
}

#endif
