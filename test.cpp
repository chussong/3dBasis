#include "test.hpp"

namespace Test {

bool RunAllTests() {
    Multinomial::Initialize(1, 6);
    Multinomial::Initialize(2, 6);
    Multinomial::Initialize(3, 6);
    Multinomial::Initialize(4, 6);
    Multinomial::Initialize(5, 6);
    Multinomial::Initialize(6, 6);

    std::cout << "----- PERFORMING ALL AVAILABLE UNIT TESTS -----" << std::endl;
    bool result = true;
    result &= MatrixInternal::PermuteXY();
    result &= MatrixInternal::InteractionTermsFromXY();
    result &= MatrixInternal::CombineInteractionFs();
    result &= RIntegral();

    int numP = 3;
    int degree = 7;
    std::size_t partitions = 4;
    coeff_class partWidth = 0.5;
    std::vector<Basis<Mono>> allEvenBases;
    std::vector<Basis<Mono>> allOddBases;
    for (int deg = numP; deg <= degree; ++deg) {
            splitBasis<Mono> degBasis(numP, deg, 0);
            allEvenBases.push_back(degBasis.EvenBasis());
            allOddBases.push_back(degBasis.OddBasis());
    }
    std::vector<Poly> evenStates = ::Orthogonalize(allEvenBases, 
            partitions, partWidth, std::cout);
    std::vector<Poly> oddStates = ::Orthogonalize(allOddBases, 
            partitions, partWidth, std::cout);
    Basis<Mono> minBasis = ::MinimalBasis(evenStates);
    result &= Test::InteractionMatrix(minBasis, partitions, partWidth);

    return result;
}

namespace MatrixInternal {

bool PermuteXY() {
    std::cout << "----- MatrixInternal::PermuteXY -----" << std::endl;
    std::vector<std::string> xAndy{"10", "1000", "1010", "1100", "210000",
        "222210", "111111", "210012", "221001"};
    for (auto& xy : xAndy) {
        std::cout << "TEST CASE: " << xy << std::endl;
        do {
            std::cout << xy << std::endl;
        } while (::MatrixInternal::PermuteXY(xy));
    }

    std::cout << "----- PASSED -----" << std::endl;
    return true;
}

bool InteractionTermsFromXY() {
    std::cout << "----- MatrixInternal::InteractionTermsFromXY -----" << std::endl;
    std::vector<std::string> testCases {
        {2, 1, 0, 1, 0, 0},
        {2, 1, 0, 0, 1, 0},
        {2, 1, 0, 0, 0, 1},
        {0, 1, 2, 2, 0, 0},
        {0, 1, 2, 1, 1, 0},
        {0, 1, 2, 0, 0, 2}
    };
    for (const auto& xy : testCases) {
        std::cout << "CASE: " << MVectorOut(xy) << std::endl;
        for (const auto& term : ::MatrixInternal::InteractionTermsFromXY(xy)) {
            std::cout << term << std::endl;
        }
    }

    std::cout << "----- PASSED -----" << std::endl;
    return true;
}

bool CombineInteractionFs() {
    std::cout << "----- MatrixInternal::CombineInteractionFs -----" << std::endl;

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

        std::cout << "CASE:\n" << f1 << " *\n" << f2 << " =\n";
        auto results = ::MatrixInternal::CombineInteractionFs({f1}, {f2});
        for (const auto& res : results) {
            std::cout << res << std::endl;
        }
    }

    std::cout << "----- PASSED -----" << std::endl;
    return true;
}

} // namespace MatrixInternal

bool RIntegral() {
    std::cout << "----- ::RIntegral -----" << std::endl;
    for (coeff_class a = 1; a < 5; ++a) {
        for (coeff_class b = 1; b < 5; ++b) {
            for (coeff_class c = 1; c < 5; ++c) {
                for (coeff_class alpha = 0.125; alpha <= 1.0; alpha += 0.125) {
                    std::cout << std::vector<coeff_class>{a, b, c, alpha}
                        << " = " << ::RIntegral(a, b, c, alpha) << std::endl;
                }
            }
        }
    }
    std::cout << "----- PASSED -----" << std::endl;
    return true;
}

bool InteractionMatrix(const Basis<Mono>& basis, const std::size_t partitions,
        const coeff_class partWidth) {
    std::cout << "----- ::InteractionMatrix -----" << std::endl;
    std::cout << ::InteractionMatrix(basis, partitions, partWidth) << std::endl;
    std::cout << "----- PASSED -----" << std::endl;
    return true;
}

} // namespace Test
