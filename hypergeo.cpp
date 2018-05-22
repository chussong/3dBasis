#include "hypergeo.hpp"

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

coeff_class Hypergeometric2F1_Reg(const builtin_class a, const builtin_class b,
                                  const builtin_class c, const builtin_class x){
    static std::unordered_map< std::array<builtin_class,4>,coeff_class,
                          boost::hash<std::array<builtin_class,4>> > cache;

    const std::array<builtin_class,4> params = {{a, b, c, x}};
    if (cache.count(params) == 0) {
        cache.emplace(params, HypergeometricPFQ_Reg<2,1>({{a,b}}, {{c}}, x));
    }

    return cache[params];
}

coeff_class Hypergeometric3F2(const builtin_class a1, const builtin_class a2,
                              const builtin_class a3, const builtin_class b1,
                              const builtin_class b2, const builtin_class x) {
    return Hypergeometric3F2({{a1, a2, a3, b1, b2, x}});
}

coeff_class Hypergeometric3F2(const std::array<builtin_class,3>& a, 
                              const std::array<builtin_class,2>& b, 
                              const builtin_class x) {
    return Hypergeometric3F2({{a[0], a[1], a[2], b[0], b[1], x}});
}

coeff_class Hypergeometric3F2(const std::array<builtin_class,6>& params) {
    return Hypergeometric3F2_Reg(params) 
           * std::tgamma(params[3]) * std::tgamma(params[4]);
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

coeff_class Hypergeometric4F3(const builtin_class a1, const builtin_class a2,
                              const builtin_class a3, const builtin_class a4,
                              const builtin_class b1, const builtin_class b2,
                              const builtin_class b3, const builtin_class x) {
    return Hypergeometric4F3({{a1, a2, a3, a4, b1, b2, b3, x}});
}

coeff_class Hypergeometric4F3(const std::array<builtin_class,4>& a, 
                              const std::array<builtin_class,3>& b, 
                              const builtin_class x) {
    return Hypergeometric4F3({{a[0], a[1], a[2], a[3], b[0], b[1], b[2], x}});
}

coeff_class Hypergeometric4F3(const std::array<builtin_class,8>& params) {
    static std::unordered_map<std::array<builtin_class,8>,coeff_class,
                          boost::hash<std::array<builtin_class,8>> > cache;

    if (cache.count(params) == 0) {
        coeff_class value;
        try {
            value = HypergeometricPFQ<4,3>({{params[0],   params[1], 
                                             params[2],   params[3]}},
                                           {{params[4],   params[5],
                                             params[6]}}, params[7]);
            if (!std::isfinite(static_cast<builtin_class>(value))) {
                std::cerr << "Error: 4F3(" << params << ") = " << value << '\n';
            }
        }
        catch (const std::runtime_error& err) {
            std::cerr << "Error: 4F3(" << params << ") did not converge.\n";
            value = 0.0/0.0;
        }
        // std::cout << "Hypergeometric4F3(" << params << ") = " <<  value 
            // << '\n';
        cache.emplace(params, value);
    }

    return cache[params];
}

coeff_class Hypergeometric4F3_Reg(const builtin_class a1, const builtin_class a2,
                                  const builtin_class a3, const builtin_class a4,
                                  const builtin_class b1, const builtin_class b2,
                                  const builtin_class b3, const builtin_class x) {
    return Hypergeometric4F3_Reg({{a1, a2, a3, a4, b1, b2, b3, x}});
}

coeff_class Hypergeometric4F3_Reg(const std::array<builtin_class,4>& a, 
                                  const std::array<builtin_class,3>& b, 
                                  const builtin_class x) {
    return Hypergeometric4F3_Reg({{a[0], a[1], a[2], a[3], b[0], b[1], b[2], x}});
}

coeff_class Hypergeometric4F3_Reg(const std::array<builtin_class,8>& params) {
    static std::unordered_map<std::array<builtin_class,8>,coeff_class,
                          boost::hash<std::array<builtin_class,8>> > cache;

    if (cache.count(params) == 0) {
        coeff_class value;
        try {
            value = HypergeometricPFQ_Reg<4,3>({{params[0],   params[1], 
                                                 params[2],   params[3]}},
                                               {{params[4],   params[5],
                                                 params[6]}}, params[7]);
            if (!std::isfinite(static_cast<builtin_class>(value))) {
                std::cerr << "Error: 4F3_Reg(" << params << ") = " << value 
                    << '\n';
            }
        }
        catch (const std::runtime_error& err) {
            std::cerr << "Error: 4F3_Reg(" << params << ") did not converge.\n";
            value = 0.0/0.0;
        }
        // std::cout << "Hypergeometric4F3_Reg(" << params << ") = " <<  value 
            // << '\n';
        cache.emplace(params, value);
    }

    return cache[params];
}
