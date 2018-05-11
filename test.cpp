#include "test.hpp"

namespace Test {

bool RunAllTests(const Arguments& args) {
    Multinomial::Initialize(1, 6);
    Multinomial::Initialize(2, 6);
    Multinomial::Initialize(3, 6);
    Multinomial::Initialize(4, 6);
    Multinomial::Initialize(5, 6);
    Multinomial::Initialize(6, 6);

    OStream& console = *args.console;

    console << "----- PERFORMING ALL AVAILABLE UNIT TESTS -----" << endl;
    bool result = true;
    result &= MatrixInternal::PermuteXY(console);
    result &= MatrixInternal::InteractionTermsFromXY(console);
    result &= MatrixInternal::CombineInteractionFs(console);
    result &= MatrixInternal::UPlusIntegral(console);
    // result &= RIntegral(console);
    result &= Hypergeometric(console);

    int numP = 3;
    int degree = 7;
    // std::size_t partitions = 4;
    std::vector<Basis<Mono>> allEvenBases;
    std::vector<Basis<Mono>> allOddBases;
    for (int deg = numP; deg <= degree; ++deg) {
            splitBasis<Mono> degBasis(numP, deg, args);
            allEvenBases.push_back(degBasis.EvenBasis());
            allOddBases.push_back(degBasis.OddBasis());
    }
    std::vector<Poly> evenStates = ::Orthogonalize(allEvenBases, console, false);
    std::vector<Poly> oddStates = ::Orthogonalize(allOddBases, console, true);
    Basis<Mono> minBasis = ::MinimalBasis(evenStates);
    // result &= Test::InteractionMatrix(minBasis, args);

    result &= MuPart_NtoN(args);

    return result;
}

namespace MatrixInternal {

bool PermuteXY(OStream& console) {
    console << "----- MatrixInternal::PermuteXY -----" << endl;
    std::vector<std::string> xAndy{"10", "1000", "1010", "1100", "210000",
        "222210", "111111", "210012", "221001"};
    for (auto& xy : xAndy) {
        console << "TEST CASE: " << xy << endl;
        do {
            console << xy << endl;
        } while (::MatrixInternal::PermuteXY(xy));
    }

    console << "----- PASSED -----" << endl;
    return true;
}

bool InteractionTermsFromXY(OStream& console) {
    console << "----- MatrixInternal::InteractionTermsFromXY -----" << endl;
    std::vector<std::string> testCases {
        {2, 1, 0, 1, 0, 0},
        {2, 1, 0, 0, 1, 0},
        {2, 1, 0, 0, 0, 1},
        {0, 1, 2, 2, 0, 0},
        {0, 1, 2, 1, 1, 0},
        {0, 1, 2, 0, 0, 2}
    };
    for (const auto& xy : testCases) {
        console << "CASE: " << MVectorOut(xy) << endl;
        for (const auto& term : ::MatrixInternal::InteractionTermsFromXY(xy)) {
            console << term << endl;
        }
    }

    console << "----- PASSED -----" << endl;
    return true;
}

bool CombineInteractionFs(OStream& console) {
    console << "----- MatrixInternal::CombineInteractionFs -----" << endl;

    // 3 particle
    // std::vector<std::vector<char>> uPlusCases {
        // {5, 2}, {5, 4}, {4, 3}, {5, 2}, {5, 4}, {4, 3}, {2, 2}, {2, 4},
        // {1, 3}, {2, 2}, {2, 6}, {1, 5}, {0, 4}, {2, 4}, {1, 3}
    // };
    // std::vector<std::vector<char>> uMinusCases {
        // {3, 0}, {3, 0}, {3, 1}, {3, 0}, {3, 0}, {3, 1}, {8, 4}, {8, 4},
        // {8, 5}, {8, 4}, {8, 4}, {8, 5}, {8, 6}, {8, 4}, {8, 5}
    // };
    // std::vector<std::vector<char>> yTildeCases {
        // {1, 0}, {1, 0}, {0, 1}, {1, 0}, {1, 0}, {0, 1}, {2, 0}, {2, 0},
        // {1, 1}, {2, 0}, {2, 0}, {1, 1}, {0, 2}, {2, 0}, {1, 1}
    // };

    // 4 particle
    std::vector<std::vector<char>> uPlusCases {
        {5, 2, 2}, {5, 4, 0}, {4, 3, 2}, {5, 2, 2}, {5, 4, 0}, {4, 3, 2}, 
        {2, 2, 4}, {2, 4, 2}, {1, 3, 4}, {2, 2, 4}, {2, 6, 0}, {1, 5, 2}, 
        {0, 4, 4}, {2, 4, 2}, {1, 3, 4}
    };
    std::vector<std::vector<char>> uMinusCases {
        {3, 0, 0}, {3, 0, 0}, {3, 1, 1}, {3, 0, 2}, {3, 0, 2}, {3, 1, 1}, 
        {8, 4, 4}, {8, 4, 2}, {8, 5, 3}, {8, 4, 2}, {8, 4, 4}, {8, 5, 3}, 
        {8, 6, 0}, {8, 4, 1}, {8, 5, 2}
    };
    std::vector<std::vector<char>> yTildeCases {
        {1, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 0, 1}, {1, 0, 1}, {0, 1, 1}, 
        {2, 0, 2}, {2, 0, 0}, {1, 1, 2}, {2, 0, 1}, {2, 0, 1}, {1, 1, 0}, 
        {0, 2, 2}, {2, 0, 1}, {1, 1, 1}
    };

    std::vector<std::size_t> indices(uPlusCases.size());
    for (std::size_t i = 0; i < indices.size(); ++i) indices[i] = i;
    std::random_shuffle(indices.begin(), indices.end());

    for (std::size_t i = 0; i < indices.size()-1; ++i) {
        ::MatrixInternal::MatrixTerm_Intermediate f1, f2;
        f1.uPlus  = uPlusCases [i];
        f1.uMinus = uMinusCases[i];
        f1.yTilde = yTildeCases[i];
        f2.uPlus  = uPlusCases [i+1];
        f2.uMinus = uMinusCases[i+1];
        f2.yTilde = yTildeCases[i+1];

        console << "CASE:\n" << f1 << " *\n" << f2 << " =\n";
        auto results = ::MatrixInternal::CombineInteractionFs({f1}, {f2});
        for (const auto& res : results) {
            console << res << endl;
        }
    }

    console << "----- PASSED -----" << endl;
    return true;
}

bool UPlusIntegral(OStream& console) {
    console << "----- ::UPlusIntegral -----" << endl;
    bool passed = true;
    passed &= UPlusIntegral_Case(5, 2, 0.0634921, console);
    passed &= UPlusIntegral_Case(2, 4, 0.0833333, console);
    passed &= UPlusIntegral_Case(8, 8, 0.0015873, console);

    if (passed) {
        console << "----- PASSED -----" << endl;
        return true;
    } else {
        console << "----- FAILED -----" << endl;
        return false;
    }
}

bool UPlusIntegral_Case(const builtin_class a, const builtin_class b, 
        const builtin_class expected, OStream& console) {
    constexpr builtin_class tol = 1e-5;
    builtin_class answer = ::MatrixInternal::UPlusIntegral(a, b);
    console << "UPlusIntegral(" << a << ", " << b << ") == " << answer;
    if (std::abs(answer - expected) <= tol*answer) {
        console << " == " << expected << " (PASS)" << endl;
        return true;
    } else {
        console << " != " << expected << " (FAIL)" << endl;
        // std::cerr << " (beta: " << gsl_sf_beta(a/2 + 1, b/2 + 1) << ")" << std::endl;
        return false;
    }
}

} // namespace MatrixInternal

// bool RIntegral(OStream& console) {
    // console << "----- ::RIntegral -----" << std::endl;
    // for (coeff_class a = 1; a < 5; ++a) {
        // for (coeff_class b = 1; b < 5; ++b) {
            // for (coeff_class c = 1; c < 5; ++c) {
                // for (coeff_class alpha = 0.125; alpha <= 1.0; alpha += 0.125) {
                    // console << std::vector<coeff_class>{a, b, c, alpha}
                        // << " = " << ::RIntegral(a, b, c, alpha) << std::endl;
                // }
            // }
        // }
    // }
    // console << "----- PASSED -----" << std::endl;
    // return true;
// }

bool Hypergeometric(OStream& console) {
    console << "----- ::HypergeometricPFQ -----" << endl;
    bool passed = true;

    passed &= HypergeometricPFQ_Case<2,1>({{1,2}}, {{3}}, 0.4, 1.38532, console);
    passed &= HypergeometricPFQ_Reg_Case<2,1>({{1,2}}, {{3}}, 0.4, 0.69266, console);
    passed &= HypergeometricPFQ_Case<2,1>({{1,-2}}, {{3}}, 0.4, 0.76, console);
    passed &= HypergeometricPFQ_Reg_Case<2,1>({{1,-2}}, {{3}}, 0.4, 0.38, console);
    passed &= HypergeometricPFQ_Reg_Case<2,1>({{1,-2}}, {{-3}}, 0.4, 0.0, console);
    passed &= HypergeometricPFQ_Reg_Case<2,1>({{1,2}}, {{-3}}, 0.4, 65.8436, console);

    passed &= HypergeometricPFQ_Case<3,2>({{1,2,3}}, {{4,5}}, 0.6, 1.24177, console);
    passed &= HypergeometricPFQ_Reg_Case<3,2>({{1,2,3}}, {{4,5}}, 0.6, 0.00862338, console);
    passed &= HypergeometricPFQ_Case<3,2>({{1,-2,3}}, {{4,5}}, 0.6, 0.8344, console);
    passed &= HypergeometricPFQ_Reg_Case<3,2>({{1,-2,3}}, {{4,5}}, 0.6, 0.00579444, console);
    passed &= HypergeometricPFQ_Reg_Case<3,2>({{1,-2,3}}, {{-4,5}}, 0.6, 0.0, console);
    passed &= HypergeometricPFQ_Reg_Case<3,2>({{1,2,3}}, {{-4,5}}, 0.6, 58.4183, console);
    passed &= HypergeometricPFQ_Reg_Case<3,2>({{1,2,-3}}, {{-4,-5}}, 0.6, 0.0, console);
    passed &= HypergeometricPFQ_Reg_Case<3,2>({{1,2,3}}, {{-4,-5}}, 0.6, 4.61311e14, console);

    // argument x=1 requires special treatment that's not implemented yet
    // passed &= HypergeometricPFQ_Case<2,1>({{1,2}}, {{4}}, 1.0, 3.0, console);
    // passed &= HypergeometricPFQ_Reg_Case<2,1>({{1,2}}, {{4}}, 1.0, 0.5, console);
    // passed &= HypergeometricPFQ_Case<3,2>({{1,2,3}}, {{4,5}}, 1.0, 1.56475, console);
    // passed &= HypergeometricPFQ_Reg_Case<3,2>({{1,2,3}}, {{4,5}}, 1.0, 0.0108663, console);

    if (passed) {
        console << "----- PASSED -----" << endl;
    } else {
        console << "----- FAILED -----" << endl;
    }
    return passed;
}

bool InteractionMatrix(const Basis<Mono>& basis, const Arguments& args) {
    OStream& console = *args.console;
    console << "----- ::InteractionMatrix -----" << endl;
    console << ::InteractionMatrix(basis, args.partitions) << endl;
    console << "----- PASSED -----" << endl;
    return true;
}

bool MuPart_NtoN(const Arguments& args) {
    OStream& console = *args.console;
    console << "----- ::MuPart_NtoN -----" << endl;
    DMatrix muPart = ::MuPart_NtoN(3, {{0, 0}}, 5);
    DMatrix reference(5, 5);
    reference << 0.2385140, 0.1321650, 0.0947868, 0.0783601, 0.0683649,
                 0.1321650, 0.1717750, 0.1140410, 0.0866316, 0.0733793,
                 0.0947868, 0.1140410, 0.1470790, 0.1032710, 0.0808491,
                 0.0783601, 0.0866316, 0.1032710, 0.1322320, 0.0957497,
                 0.0683649, 0.0733793, 0.0808491, 0.0957497, 0.1218690;
    DMatrix quotient = muPart.array() / reference.array();
    if (quotient.isOnes(1e-5)) {
        console << "----- PASSED -----" << endl;
        return true;
    } else {
        console << muPart << "\n^muPart not equal to reference_\n" << reference
            << "\nratio:\n" << quotient << '\n';
        console << "----- FAILED -----" << endl;
        return false;
    }
}

} // namespace Test
