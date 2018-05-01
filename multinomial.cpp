#include "multinomial.hpp"

namespace Multinomial {

namespace {
    std::vector<std::unique_ptr<MultinomialTable>> multinomialTable;
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

void Initialize(const char particleNumber, const char highestN) {
    if (multinomialTable.size() <= static_cast<std::size_t>(particleNumber)) {
        multinomialTable.resize(particleNumber+1);
    }
    if (multinomialTable[particleNumber] == nullptr) {
        multinomialTable[particleNumber] = 
                            std::make_unique<MultinomialTable>(particleNumber);
    }
    /*if (particleNumber >= 2)*/ multinomialTable[particleNumber]->FillTo(highestN);
}

void Clear() {
    multinomialTable.clear();
}

std::unique_ptr<MultinomialTable>& GetTable(const std::size_t n, const char d) {
    if (multinomialTable.size() < n+1 
        || multinomialTable[n] == nullptr
        || multinomialTable[n]->HighestN() < d) {
        // throw std::logic_error("GetMVectors called on uninitialized pNumber.");
        Initialize(n, d);
    }
    return multinomialTable[n];
}

// vector of all mVectors whose total "n" is exactly the supplied n
//
// if this turns out to be slow, we can avoid the copy by passing iterators to
// a slightly reorganized container for the mVectors
MVectorContainer GetMVectors(const unsigned char particleNumber, const char n) {
    return GetTable(particleNumber, n)->GetMVectors(n);
}

// binomial coefficient (n, m)
coeff_class Choose(const char n, const char m) {
    // FIXME? if this is slow, special case for binomial without vector overhead
    return GetTable(2, n)->Lookup(std::string({{n, static_cast<char>(n-m), m}}));
}

// multinomial coefficient (n, \vec m)
coeff_class Choose(const char particleNumber, const char n, 
                   const std::vector<char>& m) {
    return GetTable(particleNumber, n)->Choose(n, m);
}

coeff_class Lookup(const char particleNumber, const std::string& nAndm) {
    return GetTable(particleNumber, nAndm[0])->Lookup(nAndm);
}

/*void FillTo(const char particleNumber, const char newHighestN) {
    if (GetTable(particleNumber, newHighestN) == nullptr) {
        std::cerr << "Error: asked to fill a nullptr MultinomialTable to n="
            << std::to_string(newHighestN) << "." << std::endl;
        return;
    }
    multinomialTable[particleNumber]->FillTo(newHighestN);
}*/

MultinomialTable::MultinomialTable(const unsigned char particleNumber): 
	particleNumber(particleNumber), highestN(-1) {
    // if (particleNumber < 0) {
        // std::cerr << "Warning: MultinomialTable has been constructed with "
            // << "particleNumber " << std::to_string(particleNumber) 
            // << " < 0; your program is probably about to crash." << std::endl;
    // }
}

// we fill this by constructing Pascal's simplex, in which each entry is the sum
// of all entries on the previous level which can reach it
//
// we will not store any entries with the terms out of order, which is a 
// generalization of the policy we used for binomials. We will also not store
// any whose first term is the degree, since these are all 1 
void MultinomialTable::FillTo(const char newHighestN) {
    if (highestN >= newHighestN) return;

    // first, generate the mVectors we need to compute the coefficients for
    ComputeMVectors(newHighestN);
    //std::cout << "mVectors constructed, now size " << mVectors.size() << std::endl;

    table.emplace(std::string(particleNumber+1, 0), 1);

    // the case with a single particle has to be handled differently because it
    // wrecks the main loop
    if (particleNumber == 1) {
        for (const std::string& key : mVectors) table.emplace(key, 1);
        highestN = newHighestN;
        return;
    }

    // each key is the sum of the values of all lower keys which can produce it
    for (std::string key : mVectors) {
        if (key[0] == 0) {
            //std::cout << "Emplacing " << MVectorOut(key) << " with value 1."
                    //<< std::endl;
            //table.emplace(std::move(key), 1);
            continue;
        }
        // std::cout << "Evaluating this key based on lower values: " 
                // << MVectorOut(key) << std::endl;
        highestN = std::max(highestN, key[0]);
        coeff_class value = 0;
        --key[0];
        --key[1];
        value += Lookup(key);
        int i = 2;
        while (i < particleNumber+1 && key[i] != 0) {
            ++key[i-1];
            --key[i];
            try {
                value += Lookup(key);
            }
            catch (std::out_of_range) {
                std::cerr << "Error: attempted to look up " << MVectorOut(key)
                        << " but it did not exist." << std::endl;
                throw;
            }
            ++i;
        }
        //--key[1];
        ++key[0];
        ++key[i-1];
        if (value == 0) value = 1;
        // std::cout << "Emplacing this key: " << MVectorOut(key) << 
            // " with value: " << value << std::endl;
        table.emplace(std::move(key), value);
    }
    highestN = newHighestN;
}

// vector of all mVectors whose total "n" is exactly the supplied n
//
// if this turns out to be slow, we can avoid the copy by passing iterators to
// a slightly reorganized container for the mVectors
MVectorContainer MultinomialTable::GetMVectors(const unsigned char n) {
    // std::cout << "Returning mVectors corresponding to n = " << (int)n
        // << " from among the " << mVectors.size() << " known ones." << std::endl;
    if (n == 0) {
        // std::cout << "Someone is requesting the 0 mVector in "
            // << int(particleNumber) << "particles." << std::endl;
        return {std::string(particleNumber+1, 0)};
    }
    std::string lower(particleNumber+1, n);
    std::string upper(particleNumber+1, 0);
    upper[0] = n;
    //std::cout << "Returning mVectors between " << MVectorOut(lower) << " and " 
            //<< MVectorOut(upper) << std::endl;
    // the below is obviously only for safety, remove it if needed
    if (MVectorContainer(mVectors.lower_bound(lower), 
                         mVectors.upper_bound(upper)).size() == 0) {
        std::cerr << "WARNING: request for mVectors at (" << (int)particleNumber
            << ", " << (int)n << ") returned an empty set drawn from the total "
            << "of " << mVectors.size() << " mVectors with n up to "
            << (int)highestN << "." << std::endl;
    }
    return MVectorContainer(mVectors.lower_bound(lower), 
                            mVectors.upper_bound(upper));
}

// computes all mVectors with n up to newHighestN
void MultinomialTable::ComputeMVectors(const char newHighestN) {
    // std::cout << "Computing mVectors for n = " << (int)newHighestN << std::endl;
    mVectors.clear();
    for (char n = newHighestN; n >= 0; --n) {
        std::string mVector(particleNumber+1, 0);
        mVector[0] = n;
        mVector[1] = n;
        do { 
            // std::cout << "Emplacing an mVector: " << MVectorOut(mVector) 
                    // << std::endl;
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
    // bool debug = false;
    // if (MVectorOut(mVector) == "A622") {
        // std::cout << "Incrementing A622." << std::endl;
        // debug = true;
    // }
    for(unsigned int i = mVector.size()-2; i > 0; --i){
        if (mVector[i] < 2) continue;
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
                        while (mVector[l] != 0 && mVector[k-1] > mVector[k]) {
                            ++mVector[k];
                            --mVector[l];
                        }
                    }
                    if (mVector[k] == 0) break;
                }
                // if (debug) {
                    // std::cout << "Changed it to " << MVectorOut(mVector) 
                        // << std::endl;
                // }
                return true;
            } else if(mVector[i] == mVector[j] + 1) {
                // mV[j] can't be increased but maybe a future one can
                ++j;
            } else {
                // mV[i] == mV[i+1] so we already know it can't be decreased
                break;
            }
        } while (j < mVector.size());
    }
    return false;
}

coeff_class MultinomialTable::Choose(const char n, const std::vector<char>& m) {
    return Lookup(n + std::string(m.begin(), m.end()));
}

// the first entry of nAndm is n, followed by the m vector
coeff_class MultinomialTable::Lookup(std::string nAndm) {
    // std::cout << "Lookup of mVector " << MVectorOut(nAndm) << " with length "
        // << nAndm.size() << std::endl;
    if (nAndm.size() != particleNumber + 1) {
        std::cerr << "Error: asked to choose a multinomial with an m vector "
            << "whose size (" << nAndm.size()-1 << ") is different from the "
            << "number of particles (" << std::to_string(particleNumber) 
            << ")." << std::endl;
        return 0;
    }

    // we're using the permutation symmetry to only store sorted coefficients
    std::sort(nAndm.begin()+1, nAndm.end(), std::greater<char>());
    if (nAndm[0] <  nAndm[1]) return 0;
    // the below line should not be necessary, since the zero vectors should be
    // emplaced in the table. I have no idea what the problem is
    // if (nAndm[0] == nAndm[1]) return 1;
    //if (nAndm[0] == nAndm[1] + 1) return nAndm[0];

    if (nAndm[0] > highestN) {
        // std::cerr << "Error: requested multinomial(" 
            // << std::to_string(nAndm[0]) << ", m), but this table was only "
            // << "filled up to " << std::to_string(highestN) << "." << std::endl;
        // return -1;
        FillTo(nAndm[0]);
    }

    // std::cout << "Looking up this mVector: " << MVectorOut(nAndm) << std::endl;
    return table.at(nAndm);
}

} // namespace Multinomial
