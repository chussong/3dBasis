#ifndef CACHE_HPP
#define CACHE_HPP

#include <cmath>
//#include <iostream> // debug only

#include "constants.hpp"

class GammaCache {
	public:
		GammaCache(const int particleNumber, const int maxPm, const int maxPt);
		void Update(const int maxPm, const int maxPt);

		coeff_class Prefactor(const int totalPm, const int totalPt) const;
		coeff_class& Prefactor(const int totalPm, const int totalPt);
		coeff_class Middle(const int totalPm, const int totalPt,
				const int totalK) const;
		coeff_class& Middle(const int totalPm, const int totalPt,
				const int totalK);
		coeff_class Inner(const int pm, const int pt, const int k) const;
		coeff_class& Inner(const int pm, const int pt, const int k);

	private:
		const coeff_class particleNumber;
		int highestKnownPm;
		int highestKnownPt;
		coeff_class commonPrefactor;
		std::vector<coeff_class> prefactorCache;
		std::vector<coeff_class> middleGammaCache;
		std::vector<coeff_class> innerLogCache;
		int middlePmCoeff;
		int innerPmCoeff;

		void UpdatePrefactors(const int maxPm, const int maxPt);
		coeff_class GeneratePrefactor(const int maxPm, const int maxPt) const;

		void UpdateMiddle(const int maxPm, const int maxPt);
		coeff_class GenerateMiddle(const int maxPm, const int maxPt,
									const int totalK) const;

		void UpdateInner(const int maxPm, const int maxPt);
		coeff_class GenerateInner(const int pm, const int pt, const int k)const;
};

#endif
