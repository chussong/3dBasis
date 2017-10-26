#ifndef BINOMIAL_HPP
#define BINOMIAL_HPP

#include <vector>
#include <iostream>

#include "constants.hpp" // for coeff_class

coeff_class ExactBinomial(const char n, const char k);
void ExactBinomial_FillTo(const char newHighestN);

class BinomialTable {
	public:
		BinomialTable();

		coeff_class Choose(const char n, char k);
		void FillTo(const char newHighestN);

	private:
		std::vector<coeff_class> table;
		char highestN;

		coeff_class Lookup(const char n, const char k) const;
};

#endif
