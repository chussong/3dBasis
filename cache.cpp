#include "cache.hpp"

// generates a cache of all math function outputs used in Zuhair's inner product
// (most of which are gamma functions)
//
// This could be smarter: we could just have a degree instead of maxPm and
// maxPt, which would allow us to e.g. not generate the (8, 6) coefficients for
// a degree 8 computation (since those are actually degree 14)
GammaCache::GammaCache(const int particleNumber, const int maxPm, 
		const int maxPt): particleNumber(particleNumber), highestKnownPm(0),
		highestKnownPt(0) {
	commonPrefactor = M_PI*M_PI/std::pow(4*M_PI, particleNumber);
	Update(maxPm, maxPt);
}

void GammaCache::Update(const int maxPm, const int maxPt) {
	if(maxPm <= highestKnownPm && maxPt <= highestKnownPt) return;
	highestKnownPm = maxPm;
	highestKnownPt = maxPt;

	UpdatePrefactors(maxPm, maxPt);
	UpdateMiddle(maxPm, maxPt);
	UpdateInner(maxPm, maxPt);
}

void GammaCache::UpdatePrefactors(const int maxPm, const int maxPt) {
	prefactorCache.resize((maxPm - particleNumber + 1)*(maxPt/2 + 1));
	//size_t DEBUG_NUMFILLED = 0;
	for(int totalPm = particleNumber; totalPm <= maxPm; ++totalPm){
		for(int totalPt = 0; totalPt <= maxPt; totalPt += 2){
			Prefactor(totalPm, totalPt) = GeneratePrefactor(totalPm, totalPt);
			//std::cout << "Prefactor(" << totalPm << ", " << totalPt << ") = "
				//<< Prefactor(totalPm, totalPt) << std::endl;
			//++DEBUG_NUMFILLED;
		}
	}
	//std::cout << "Prefactor cache updated. Filled " << DEBUG_NUMFILLED << "/"
		//<< prefactorCache.size() << " values." << std::endl;
}

coeff_class GammaCache::GeneratePrefactor(const int totalPm, const int totalPt)
																		const {
	// particleNumber is a coeff_class here, so none of these divisions round
	coeff_class current = -std::lgamma(totalPm + totalPt + particleNumber/2);
	current -= std::lgamma((totalPt + particleNumber - 1)/2);
	current -= std::lgamma(totalPm + (totalPt + particleNumber - 1)/2);
	current = std::exp(current);
	current *= std::pow(2, -2*totalPm - totalPt - particleNumber + 4);
	return commonPrefactor*current;
}

coeff_class GammaCache::Prefactor(const int totalPm, const int totalPt) const {
	return prefactorCache[(totalPm - particleNumber) 
		+ (highestKnownPm - particleNumber)*(totalPt/2)];
}

coeff_class& GammaCache::Prefactor(const int totalPm, const int totalPt) {
	return prefactorCache[(totalPm - particleNumber) 
		+ (highestKnownPm - particleNumber)*(totalPt/2)];
}

void GammaCache::UpdateMiddle(const int maxPm, const int maxPt) {
	middlePmCoeff = ((maxPt/2 + 1)*(maxPt/2 + 2))/2;

	middleGammaCache.resize((maxPm-particleNumber+1)*middlePmCoeff);
	//size_t DEBUG_NUMFILLED = 0;
	for(int totalPm = particleNumber; totalPm <= maxPm; ++totalPm){
		for(int totalPt = 0; totalPt <= maxPt; totalPt += 2){
			for(int totalK = 0; totalK <= totalPt/2; ++totalK){
				Middle(totalPm, totalPt, totalK) = 
									GenerateMiddle(totalPm, totalPt, totalK);
				//std::cout << "Middle(" << totalPm << ", " << totalPt << ", "
					//<< totalK << ") = " << Middle(totalPm, totalPt, totalK) 
					//<< std::endl;
				//++DEBUG_NUMFILLED;
			}
		}
	}
	//std::cout << "Middle gamma updated. Filled " << DEBUG_NUMFILLED << "/"
		//<< middleGammaCache.size() << " values." << std::endl;
}

coeff_class GammaCache::GenerateMiddle(const int totalPm, const int totalPt,
										const int totalK) const {
	coeff_class ret = std::lgamma(totalPm + (totalPt + particleNumber - 1)/2 
									+ totalK);
	ret += std::lgamma((totalPt + 1.0)/2 - totalK);
	ret = std::exp(ret) / (1 << 2*totalK);
	//ret /= 1 << 2*totalK;
	//if((totalPt/2-totalK) % 2 == 1) ret = -ret;
	return (totalPt/2 - totalK) % 2 == 1 ? -ret : ret;
}

coeff_class GammaCache::Middle(const int totalPm, const int totalPt,
		const int totalK) const {
	return middleGammaCache[(totalPm-particleNumber)*middlePmCoeff
		+ ((totalPt/2)*(1 + totalPt/2))/2
		+ totalK];
}

coeff_class& GammaCache::Middle(const int totalPm, const int totalPt,
		const int totalK) {
	//std::cout << "Accessing Middle[" << (totalPm-particleNumber)*middlePmCoeff
		//+ ((totalPt/2)*(1 + totalPt/2))/2 + totalK << "] from (" << totalPm
		//<< ", " << totalPt << ", " << totalK << ")" << std::endl;
	return middleGammaCache[(totalPm-particleNumber)*middlePmCoeff
		+ ((totalPt/2)*(1 + totalPt/2))/2
		+ totalK];
}

// pm goes up to maxPm, pt goes up to maxPt, k goes up to pt/2
void GammaCache::UpdateInner(const int maxPm, const int maxPt) {
	innerPmCoeff = 1 + maxPt + ( (maxPt-1)/2 + ((maxPt-1)/2)*((maxPt-1)/2) )/2
					+ ( maxPt/2 + (maxPt/2)*(maxPt/2) )/2;

	innerLogCache.resize((maxPm+1)*innerPmCoeff);
	//size_t DEBUG_NUMFILLED = 0;
	for(int pm = 0; pm <= maxPm; ++pm){
		for(int pt = 0; pt <= maxPt; ++pt){
			for(int k = 0; k <= pt/2; ++k){
				Inner(pm, pt, k) = GenerateInner(pm, pt, k);
				//++DEBUG_NUMFILLED;
			}
		}
	}
	//std::cout << "Inner gamma updated. Filled " << DEBUG_NUMFILLED << "/"
		//<< innerLogCache.size() << " values." << std::endl;
}

coeff_class GammaCache::GenerateInner(const int pm, const int pt, const int k) 
																		const {
	coeff_class ret = std::lgamma(2*pm + pt + 1);
	ret += std::lgamma(pt + 1);
	ret -= std::lgamma(k + 1);
	ret -= std::lgamma(pm + k + 1);
	ret -= std::lgamma(pt - 2*k + 1);
	return ret;
}

coeff_class GammaCache::Inner(const int pm, const int pt, const int k) const {
	return innerLogCache[pm*innerPmCoeff + pt - pt/2 
		+ ( (pt-1)/2 + ((pt-1)/2)*((pt-1)/2) )/2 + (pt/2 + (pt/2)*(pt/2))/2
		+ k];
}

coeff_class& GammaCache::Inner(const int pm, const int pt, const int k) {
	//std::cout << "Accessing Inner[" << pm*innerPmCoeff + pt - pt/2 
		//+ ( (pt-1)/2 + ((pt-1)/2)*((pt-1)/2) )/2 + (pt/2 + (pt/2)*(pt/2))/2
		//+ k << "] from (" << pm << ", " << pt << ", " << k << ")" << std::endl;
	return innerLogCache[pm*innerPmCoeff + pt - pt/2 
		+ ( (pt-1)/2 + ((pt-1)/2)*((pt-1)/2) )/2 + (pt/2 + (pt/2)*(pt/2))/2
		+ k];
}

// creates an empty KVectorBundle to represent an illegal configuration
KVectorBundle::KVectorBundle() {
}

// cfgVector is a vector of (Pt[i] + Pt[j])/2 for each contracted pair (i,j),
// except that its last entry is the totalK to be used
KVectorBundle::KVectorBundle(const std::vector<char>& cfgVector) {
	//std::cout << "cfgVector: " << cfgVector << std::endl;
	//std::cout << "Abbreviated cfgVector: "
		//<< std::vector<char>(cfgVector.begin(), cfgVector.end()-1) << std::endl;
	kVectors = Generate(cfgVector.back(), 
			std::vector<char>(cfgVector.begin(), cfgVector.end()-1), 0);
	//std::cout << kVectors.size() << " entries in kVectors:" << std::endl;
	//for(auto it = kVectors.begin(); it != kVectors.end(); it += cfgVector.size()-1){
		//std::cout << std::vector<char>(it, it+cfgVector.size()-1) << std::endl;
	//}
}

// cfgVector contains the (Pt[i] + Pt[j])/2 for each contraction, representing
// the maximum possible K that can be put there. The former final entry 
// containing totalK has been removed and passed separately.
std::vector<char> KVectorBundle::Generate(const char totalK, 
		const std::vector<char>& cfgVector, const size_t start) {
	std::vector<char> ret;
	if(start >= cfgVector.size()) return ret;
	const char maxK = cfgVector[start];
	// if we're filling the last spot, return empty if we can't fit totalK in it
	if(start == cfgVector.size()-1 && totalK >  maxK) return ret;
	if(start == cfgVector.size()-1 && totalK <= maxK) return {totalK};

	for(char thisK = 0; thisK <= std::min(maxK, totalK); ++thisK){
		std::vector<char> remainingCfgs = Generate(totalK - thisK, cfgVector, 
				start+1);
		size_t remainingSize = cfgVector.size() - start - 1;
		if(remainingCfgs.size() % remainingSize != 0){
			std::cerr << "Error: KVectorBundle::Generate has the wrong "
				<< "remainingSize (" << remainingSize << "), which conflicts "
				<< "with the size of next level's vector (" 
				<< remainingCfgs.size() <<")." << std::endl;
			return ret;
		}
		// if this is slow, we can resize ret beforehand and then fill instead
		// of pushing back and inserting
		for(auto it = remainingCfgs.begin(); it != remainingCfgs.end(); 
													it += remainingSize){
			ret.push_back(thisK);
			ret.insert(ret.end(), it, it + remainingSize);
		}
	}
	return ret;
}

char KVectorBundle::operator[] (const size_t index) const {
	return kVectors[index];
}

// totalPt is the total (Pt[i] + Pt[j])/2 of the configuration; the cfgVector
// is the arrangement of (Pt[i] + Pt[j])/2, except for the last entry which is
// totalK.
KVectorCache::KVectorCache(const int particleNumber, const int maxPt) {
	for(char totalPt = 0; totalPt <= maxPt; ++totalPt){
		for(int totalK = 0; totalK <= totalPt; ++totalK){
			std::vector<char> cfgVector(particleNumber + 1, 0);
			cfgVector[0] = totalPt;
			cfgVector[particleNumber] = totalK;
			do{
				bundles.emplace_front(cfgVector);
				do{
					bundleMap.emplace(cfgVector, bundles.front());
				}while(std::prev_permutation(cfgVector.begin(), cfgVector.end()-1));
			}while(NextSortedCfg(cfgVector));
		}
	}
}

// the last entry contains the totalK and should not be changed
bool KVectorCache::NextSortedCfg(std::vector<char>& cfgVector){
	for(size_t pos = cfgVector.size()-2; pos > 0; --pos){
		if(cfgVector[pos-1] >= cfgVector[pos] + 2){
			--cfgVector[pos-1];
			++cfgVector[pos];
			return true;
		}
	}
	for(size_t pos = 1u; pos < cfgVector.size()-1; ++pos){
		cfgVector[0] += cfgVector[pos];
		cfgVector[pos] = 0;
	}
	return false;
}

// ptVector is a vector containing the (Pt[i] + Pt[j])/2 for each contraction
const KVectorBundle& KVectorCache::FromPt(std::vector<char> ptVector,
		const char totalK) const {
	int maxK = 0;
	for(char contraction : ptVector) maxK += contraction;
	if(maxK < totalK) return nullBundle;

	ptVector.push_back(totalK);
	if(bundleMap.find(ptVector) == bundleMap.end()){
		std::cerr << "Error: attempted to look up ptVector = " << ptVector 
			<< ", but found nothing." << std::endl;
		throw std::runtime_error("KVectorCache::FromPt");
	}
	return bundleMap.find(ptVector)->second;
}

// compare two vectors of ints, returning true if A is "less than" B in a sense	
// equivalent to if A and B were written as digits of a number in base infinity,
// where A.front() is the most significant and A.back() is the least significant
bool KVectorCache::CfgLessThan::operator() (const std::vector<char>& A, 
		const std::vector<char>& B) const {
	if(A.size() != B.size()) throw std::runtime_error("CfgLessThan size error");
	for(size_t i = 0; i < A.size(); ++i){
		if(A[i] == B[i]) continue;
		return A[i] < B[i];
	}
	return false; // they're equal, so A is not less than B
}
