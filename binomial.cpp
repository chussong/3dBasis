#include "binomial.hpp"

namespace {
	BinomialTable binomialTable;
} // anonymous namespace

coeff_class ExactBinomial(const char n, const char k) {
	return binomialTable.Choose(n, k);
}

void ExactBinomial_FillTo(const char newHighestN) {
	binomialTable.FillTo(newHighestN);
}

// populates the first two rows automatically to make FillTo smoother
BinomialTable::BinomialTable(): table({6, 10}), highestN(5) {
}

coeff_class BinomialTable::Choose(const char n, char k) {
	if (k > n) return 0;
	if (k == 0 || k == n) return 1;
	if (k == 1 || k == n-1) return n;

	if (n > highestN) {
		std::cerr << "Warning: requested binomial(" << n << ", " << k << "),"
			<< " but this table was only filled up to " << highestN << "."
			<< std::endl;
		FillTo(n);
	}

	// we don't want to store identical results, so stick to left side of table
	k = std::min(k, static_cast<char>(n-k));

	return Lookup(n, k);
}

coeff_class BinomialTable::Lookup(const char n, const char k) const {
	// we're not storing the two outer layers of the triangle, so our indexing
	// starts from (4 choose 2) and we adjust n and k accordingly
	char q = n/2;
	char r = n%2;
	return table[q*q + r*q - 3*q - r + k];
}

// Recursively fill the table up to n = newHighestN by adding the previous rows
//
// NOTE: THIS FUNCTION DE FACTO REQUIRES THAT HIGHESTN BE INITIALIZED TO 4 AND
// THAT THE FIRST ROW BE POPULATED AUTOMATICALLY
void BinomialTable::FillTo(const char newHighestN) {
	if (newHighestN < 4) {
		std::cerr << "Error: asked to fill a BinomialTable to " 
			<< std::to_string(newHighestN) << ", which is less than the "
			<< "minimum, 4." << std::endl;
		return;
	}
	if (newHighestN < highestN) {
		std::cerr << "Warning: askes to fill a BinomialTable to " 
			<< std::to_string(newHighestN) << ", but it's already filled to " 
			<< std::to_string(highestN) << ". Will do nothing." << std::endl;
		return;
	}

	auto i = table.size(); // index of first entry to update
	auto j = i - highestN/2 + 1; // index of triangle entry above-right of i
	char q = (newHighestN+1)/2;
	char r = (newHighestN+1)%2;
	table.resize(q*q + r*q - 3*q - r);

	for (auto n = highestN+1; n <= newHighestN; ++n) {
		// the term on the left edge is special
		table[i] = (n-1) + table[j];
		++i; // advance i without advancing j to make j be the above-left term
		
		// do the terms in the middle (or all of them in odd-numbered rows)
		for (auto k = 3; k < (n-1)/2; ++k) {
			table[i] = table[j] + table[j+1];
			++i;
			++j;
		}

		// in even-numbered rows, the term on the right edge is special again
		if (n%2 == 0) {
			table[i] = 2*table[j];
			++i;
		}
		++j; // advance j on its own to return it to above-right for next row
	}
}
