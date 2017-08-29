#ifndef CACHE_HPP
#define CACHE_HPP

#include <cmath>
#include <iostream> // debug only
#include <forward_list>
#include <algorithm> // std::prev_permutation

#include "constants.hpp"
#include "io.hpp"

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

struct KVectorBundle {
	//explicit KVectorBundle(const std::vector<int>& cfgVector);
	KVectorBundle();
	KVectorBundle(const std::vector<char>& cfgVector);

	char operator() (const size_t kVec, const size_t particle) const;
	char operator[] (const size_t position) const;

	std::vector<char>::const_iterator		begin() const noexcept
			{ return kVectors.begin(); }
	std::vector<char>::iterator				begin() noexcept
			{ return kVectors.begin(); }
	std::vector<char>::const_iterator		end()	const noexcept
			{ return kVectors.end(); }
	std::vector<char>::iterator				end()	noexcept
			{ return kVectors.end(); }
	size_t size() const noexcept { return kVectors.size(); }

	private:
		std::vector<char> kVectors;
		static std::vector<char> Generate(const char totalK, 
				const std::vector<char>& cfgVector, const size_t start);
};

// maybe this should be implemented via a std::map from ptVectors to bundles?
class KVectorCache {
	public:
		KVectorCache(const int particleNumber, const int maxPt);
		const KVectorBundle& FromPt(std::vector<char> ptVector,
				const char totalK) const;
	private:
		struct CfgLessThan {
			bool operator()(const std::vector<char>& A, 
					const std::vector<char>& B) const;
		};
		//static bool CfgLessThan(const std::vector<int>& A, 
				//const std::vector<int>& B);

		//int totalPt;
		std::forward_list<KVectorBundle> bundles;
		std::map<std::vector<char>, const KVectorBundle&, CfgLessThan> bundleMap;
		KVectorBundle nullBundle;

		static bool NextSortedCfg(std::vector<char>& cfgVector);

		// we can get away with storing only one bundle per permutation class by
		// creating a map; when you generate the bundles, take the seed ptVector
		// and add all its permutations to the map, pointing back to the seed.
		// Then use the seed in the second map to look up the bundle.
		//
		// wtf just add a key for each permutation referring to the same bundle
};

#endif
