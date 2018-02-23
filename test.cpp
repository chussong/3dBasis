#include "test.hpp"

namespace Test {

bool RunAllTests() {
    bool result = true;
    result &= MatrixInternal::PermuteXY();

    return result;
}

namespace MatrixInternal {

bool PermuteXY() {
    std::cout << "TESTING MatrixInternal::PermuteXY" << std::endl;
    std::vector<std::string> xAndy{"10", "1000", "1010", "1100", "210000",
        "222210", "111111", "210012", "221001"};
    for (auto& xy : xAndy) {
        std::cout << "TEST CASE: " << xy << std::endl;
        do {
            std::cout << xy << std::endl;
        } while (::MatrixInternal::PermuteXY(xy));
    }
    return true;
}

} // namespace MatrixInternal
} // namespace Test
