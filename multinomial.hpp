#ifndef MULTINOMIAL_HPP
#define MULTINOMIAL_HPP

#include <vector>
#include <set>
#include <iostream>
#include <memory> // unique_ptr
#include <unordered_map>

#include "constants.hpp" // for coeff_class
#include "io.hpp" // MVectorOut

namespace Multinomial {

class MultinomialTable;

// sorts mVectors so their n goes up but the Ms go down
struct MVectorPrecedence {
    constexpr bool operator()(const std::string& A, const std::string& B) const;
};
typedef std::set<std::string, MVectorPrecedence> MVectorContainer;

void Initialize(const char particleNumber, const char highestN);
void Clear();
std::unique_ptr<MultinomialTable>& GetTable(const std::size_t n, const char d);
// void FillTo(const char particleNumber, const char newHighestN);
MVectorContainer GetMVectors(const unsigned char particleNumber, const char n);
coeff_class Choose(const char particleNumber, const char n, 
		const std::vector<char>& m);
coeff_class Lookup(const char particleNumber, const std::string& nAndm);

class MultinomialTable {
    public:
        explicit MultinomialTable(const unsigned char particleNumber);

        coeff_class Choose(const char n, const std::vector<char>& m);
        coeff_class Lookup(std::string nAndm);
        void FillTo(const char newHighestN);
        MVectorContainer GetMVectors(const unsigned char n);

        char HighestN() const { return highestN; }

    private:
        const unsigned char particleNumber;
        MVectorContainer mVectors;
        std::unordered_map<std::string, coeff_class> table;
        char highestN;

        void ComputeMVectors(const char newHighestN);
        static bool AdvanceMVector(std::string& mVector);
};

} // namespace Multinomial

#endif
