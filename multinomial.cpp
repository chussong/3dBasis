#include "multinomial.hpp"

namespace Multinomial {

namespace {
	std::unique_ptr<MultinomialTable> multinomialTable;
	constexpr char hexMap[16] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', 
		'9', 'A', 'B', 'C', 'D', 'E', 'F'};
} // anonymous namespace

constexpr bool MVectorPrecedence::operator()(const std::string& A, 
		const std::string& B) const {
	if (A.size() != B.size()) {
		throw (std::runtime_error("MVectorPrecedence size mismatch: " + 
					std::to_string(A.size()) + " vs " +
					std::to_string(B.size())));
	}
	if (A.size() < 2 || B.size() < 2) {
		throw (std::runtime_error("MVectorPrecedence argument(s) too short."));
	}
	if (A[0] != B[0]) return A[0] < B[0];
	for (auto i = 1u; i < A.size(); ++i) {
		if (A[i] != B[i]) return A[i] > B[i];
	}
	return false;
}

std::string MVectorOut(std::string mVector) {
	for (auto i = 0u; i < mVector.size(); ++i) {
		if (mVector[i] > 15) throw(std::out_of_range("MVectorOut range error"));
		mVector[i] = hexMap[static_cast<int>(mVector[i])];
	}
	return mVector;
}

void Initialize(const char particleNumber, const char highestN) {
	multinomialTable = std::make_unique<MultinomialTable>(particleNumber);
	if (particleNumber >= 2) FillTo(highestN);
}

void Clear() {
	multinomialTable = nullptr;
}

// vector of all mVectors whose total "n" is exactly the supplied n
//
// if this turns out to be slow, we can avoid the copy by passing iterators to
// a slightly reorganized container for the mVectors
MVectorContainer GetMVectors(const char n) {
	return multinomialTable->GetMVectors(n);
}

coeff_class Choose(const char n, const std::vector<char>& m) {
	return multinomialTable->Choose(n, m);
}

coeff_class Lookup(const std::string& nAndm) {
	return multinomialTable->Lookup(nAndm);
}

void FillTo(const char newHighestN) {
	if (multinomialTable == nullptr) {
		std::cerr << "Error: asked to fill a nullptr MultinomialTable to n="
			<< std::to_string(newHighestN) << "." << std::endl;
		return;
	}
	multinomialTable->FillTo(newHighestN);
}

MultinomialTable::MultinomialTable(const char particleNumber): 
	particleNumber(particleNumber) {
	if (particleNumber < 2) {
		std::cerr << "Warning: MultinomialTable has been constructed with "
			<< "particleNumber " << std::to_string(particleNumber) 
			<< " < 2; your program is probably about to crash." << std::endl;
	}
}

// we fill this by constructing Pascal's simplex, in which each entry is the sum
// of all entries on the previous level which can reach it
//
// we will not store any entries with the terms out of order, which is a 
// generalization of the policy we used for binomials. We will also not store
// any whose first term is the degree, since these are all 1 
void MultinomialTable::FillTo(const char newHighestN) {
	ComputeMVectors(newHighestN);
	//std::cout << "mVectors constructed, now size " << mVectors.size() << std::endl;

	// each key is the sum of the values of all lower keys which can produce it
	for (std::string key : mVectors) {
		//std::cout << "Evaluating this key based on lower values: " 
			//<< MVectorOut(key) << std::endl;
		highestN = std::max(highestN, key[0]);
		coeff_class value = 0;
		--key[0];
		--key[1];
		value += Lookup(key);
		int i = 2;
		do {
			++key[i-1];
			--key[i];
			value += Lookup(key);
			++i;
		} while (i < particleNumber+1 && key[i] != 0);
		//--key[1];
		++key[0];
		++key[i-1];
		if (value == 0) value = 1;
		//std::cout << "Emplacing this key: " << MVectorOut(key) << 
			//" with value: " << value << std::endl;
		table.emplace(std::move(key), value);
	}
	highestN = newHighestN;
}

// vector of all mVectors whose total "n" is exactly the supplied n
//
// if this turns out to be slow, we can avoid the copy by passing iterators to
// a slightly reorganized container for the mVectors
MVectorContainer MultinomialTable::GetMVectors(const char n) {
	//std::cout << "Returning some of " << mVectors.size() << " mVectors." << std::endl;
	std::string lower(particleNumber+1, n);
	std::string upper(particleNumber+1, 0);
	upper[0] = n;
	//std::cout << "Returning mVectors between " << MVectorOut(lower) << " and " 
		//<< MVectorOut(upper) << std::endl;
	return MVectorContainer(mVectors.lower_bound(lower), 
			mVectors.upper_bound(upper));
}

// computes all mVectors with n up to newHighestN
void MultinomialTable::ComputeMVectors(const char newHighestN) {
	for (char n = newHighestN; n >= 0; --n) {
		std::string mVector(particleNumber+1, 0);
		mVector[0] = n;
		mVector[1] = n;
		do { 
			mVectors.emplace_hint(mVectors.begin(), mVector);
		} while(AdvanceMVector(mVector));
	}
}

// Advances mVector to the next configuration at the same n then returns true. 
// If the given mVector was the last one at this n, returns false.
bool MultinomialTable::AdvanceMVector(std::string& mVector) {
	// we're leaving entry 0 untouched because it contains n
	//std::cout << "Advancing this mVector: " << MVectorOut(mVector) << std::endl;
	// start at the end, go backward until you find an entry that can go down
	for(unsigned int i = mVector.size()-2; i > 0; --i){
		// go forward from i until you find an entry that can go up
		unsigned int j = i + 1;
		do{
			if(mVector[i] >= mVector[j] + 2){
				--mVector[i];
				++mVector[j];
				// create the lowest-sorting vector with the given mVector[i]
				for (unsigned int k = j; k < mVector.size(); ++k) {
					for (unsigned int l = mVector.size()-1; 
							mVector[k-1] > mVector[k] && l > k; --l) {
						if (mVector[l] == 0) continue;
						++mVector[k];
						--mVector[l];
					}
					if (mVector[k] == 0) break;
				}
				return true;
			} else if(mVector[i] == mVector[j] + 1) {
				++j;
			} else {
				break;
			}
		} while (j < mVector.size());
	}
	return false;
}

coeff_class MultinomialTable::Choose(const char n, const std::vector<char>& m)
	const {
	return Lookup(n + std::string(m.begin(), m.end()));
}

// the first entry of nAndm is n, followed by the m vector
coeff_class MultinomialTable::Lookup(std::string nAndm) const {
	if (nAndm.size() != particleNumber + 1) {
		std::cerr << "Error: asked to choose a multinomial with an m vector "
			<< "whose size (" << nAndm.size()-1 << ") is different from the "
			<< "number of particles (" << std::to_string(particleNumber) 
			<< ")." << std::endl;
		return 0;
	}

	// we're using the permutation symmetry to only store sorted coefficients
	std::sort(nAndm.begin()+1, nAndm.end(), std::greater<char>());
	if(nAndm[0] <  nAndm[1]) return 0;
	//if (nAndm[0] == nAndm[1]) return 1;
	//if (nAndm[0] == nAndm[1] + 1) return nAndm[0];

	if (nAndm[0] > highestN) {
		std::cerr << "Error: requested multinomial(" 
			<< std::to_string(nAndm[0]) << ", m), but this table was only "
			<< "filled up to " << std::to_string(highestN) << "." << std::endl;
		return -1;
		//FillTo(n);
	}

	return table.at(nAndm);
}

} // namespace Multinomial
