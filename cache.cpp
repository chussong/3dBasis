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
