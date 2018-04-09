#ifndef HYPERGEO_HPP
#define HYPERGEO_HPP

#include <cmath>
#include <array>
#include <iostream>

#include "constants.hpp"

// generic hypergeometric function --------------------------------------------

// FIXME: sort a & b parameter lists before calling these functions
// (they don't assume the lists are sorted, but it will make it easier to hash
// the results)

// forward declaration will be necessary for a few things
template<std::size_t P, std::size_t Q>
coeff_class HypergeometricPFQ_Reg(const std::array<builtin_class,P>&, 
        const std::array<builtin_class,Q>&, const builtin_class);

constexpr int ITERATION_LIMIT = 30000;
constexpr coeff_class PRECISION_LIMIT = 1e-10;

inline bool IsNegInt(const builtin_class x) {
    return std::round(x) < 0 && std::abs(std::round(x) - x) < EPSILON;
}

// this is an adaptation of the series representation of GSL's Hypergeometric2F1
//
// FIXME: this function should just call the one below it and de-regularize the
// result by multiplying by Prod[Gamma[b]]
/*template<int P, int Q>
coeff_class HypergeometricPFQ(const std::array<builtin_class,P>& a, 
        const std::array<builtin_class,Q>& b, const builtin_class x) {
    coeff_class reg =  HypergeometricPFQ_Reg<P,Q>(a, b, x);
    for (auto b_i : b) reg *= std::tgamma(b_i);
    return reg;

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
            del_pos = 0.0;
            del_neg = 0.0;
            break;
        }
        else {
            del_neg  = -del;
            sum_neg -=  del;
        }

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
}*/

template<std::size_t P, std::size_t Q>
coeff_class HypergeoUnitArgument_Reg(const std::array<builtin_class,P>&, 
        const std::array<builtin_class,Q>&) {
    throw std::logic_error("HypergeoUnitArgument_Reg: unspecialized (P,Q).");
}

template<>
inline coeff_class HypergeoUnitArgument_Reg(const std::array<builtin_class,2>& a, 
        const std::array<builtin_class,1>& b) {
    builtin_class num = std::tgamma(b[0] - a[0] - a[1]);
    builtin_class den = std::tgamma(b[0] - a[0]) * std::lgamma(b[0] - a[1]);
    return std::exp(num / den);
}

template<>
inline coeff_class HypergeoUnitArgument_Reg(
        const std::array<builtin_class,3>& aRef, 
        const std::array<builtin_class,2>& bRef) {
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
                return std::tgamma(a2[2])*HypergeometricPFQ_Reg<3,2>(a2, b2, 1);
            }
        }
    }

    // try the series anyway, see if it converges (if not, it'll throw)
    std::cerr << "Warning: unsolvable HypergeoUnitArgument_Reg<3,2>(" << aRef 
        << ", " << bRef << ")" << std::endl;
    return std::nan("");
}

template<std::size_t P, std::size_t Q>
coeff_class HypergeoSpecialCase_Reg(const std::array<builtin_class,P>& a,
        const std::array<builtin_class,Q>& b, const builtin_class x) {
    for (std::size_t i = 0; i < P; ++i) {
        if (a[i] == 0) {
            coeff_class gammaProd = 1;
            for (auto b_i : b) {
                if (std::abs(b_i) < EPSILON) {
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

                return HypergeometricPFQ_Reg<lowerP, lowerQ>(a2, b2, x);
            }
        }
    }

    if (std::abs(x - 1) < EPSILON) return HypergeoUnitArgument_Reg<P,Q>(a, b);

    return std::nan("");
}

template<std::size_t P, std::size_t Q>
coeff_class HypergeometricPFQ_Reg(const std::array<builtin_class,P>& a, 
        const std::array<builtin_class,Q>& b, const builtin_class x) {

    coeff_class specialCase = HypergeoSpecialCase_Reg<P,Q>(a, b, x);
    if (!std::isnan(static_cast<builtin_class>(specialCase))) return specialCase;

    coeff_class del = 1.0;
    coeff_class del_prev;
    builtin_class k = 0.0; // k is the index of the most recent COMPLETED term
    int i = 0;

    std::set<coeff_class> zeros;
    for (std::size_t i = 0; i < P; ++i) {
        if (IsNegInt(a[i])) {
            zeros.insert(std::round(-a[i]));
        }
    }
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
coeff_class HypergeometricPFQ(const std::array<builtin_class,P>& a, 
        const std::array<builtin_class,Q>& b, const builtin_class x) {
    coeff_class reg =  HypergeometricPFQ_Reg<P,Q>(a, b, x);
    for (auto b_i : b) reg *= std::tgamma(b_i);
    return reg;
}

#endif
