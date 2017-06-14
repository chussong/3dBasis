#ifndef BASIS_HPP
#define BASIS_HPP

#include <exception>
#include <algorithm>
#include <vector>
#include <array>
#include <memory>
#include <iostream>
#include <gmpxx.h>

class mono {
	mpq_class coeff;
	std::vector<int> pm;
	std::vector<int> pp;

	std::vector<int> IdentifyNodes() const;
	static std::vector<std::vector<int>> CfgsFromNodes(
			std::vector<int> nodes,
			std::vector<std::vector<int>> nodeEnergies);

	public:
		mono() = delete;
		mono(const std::vector<int> pm, const std::vector<int> pp,
				const mpq_class& coeff = 1):
			pm(pm), pp(pp), coeff(coeff) {}

		mpq_class& Coeff()		{ return coeff; }
		mpq_class Coeff() const	{ return coeff; }

		unsigned int NParticles() const { return pm.size(); }

		int Pm(const int i) const { return pm[i]; }
		int Pp(const int i) const { return pp[i]; }

		mono& operator*=(const mpq_class& x) { coeff *= x; return *this; }
		mono& operator/=(const mpq_class& x) { coeff /= x; return *this; }
		bool operator==(const mono& other) const;
		bool operator!=(const mono& other) const { return !(*this == other); }
		friend std::ostream& operator<<(std::ostream& os, const mono& out);

		mono DerivPm(const unsigned int particle);
		mono DerivPp(const unsigned int particle);
		std::vector<mono> DerivPm();
		std::vector<mono> DerivPp();

		static std::vector<std::shared_ptr<mono>> FinishPerpAtDegree(
				const std::vector<std::shared_ptr<mono>>& combined, const int degree);
};

class basis {
	int degree;
	int numP;
	std::vector<std::shared_ptr<mono>> basisMonos;

	static std::vector<std::shared_ptr<mono>> AddPlusUpToDegree(
			const std::vector<std::vector<int>>& minus, const int degree);
	static std::vector<std::shared_ptr<mono>> AddPlusUpToDegree(
			const std::vector<std::vector<int>>& minus, const int degree,
			const int M);

	public:
		basis() = delete;
		basis(const int numP, const int degree);
		basis(const int numP, const int degree, const int M);
		std::shared_ptr<mono> GetBasisMono(std::vector<int> pm,
				std::vector<int> pp);
		std::shared_ptr<mono> GetBasisMono(std::shared_ptr<mono> wildMono);
};

class poly {
	std::vector<mono> terms;

	public:
		std::vector<mono> Terms() { return terms; }

		poly& operator+=(const mono& x);
		poly& operator-=(const mono& x);
		poly& operator+=(const poly& x);
		poly& operator-=(const poly& x);

		poly DerivPm(const int);
		poly DerivPp(const int);
		poly DerivPm();
		poly DerivPp();

		poly L1() const;
		poly L2() const;
		poly L3() const;
};

std::ostream& operator<<(std::ostream& os, const std::vector<int>& out);
std::vector<std::vector<int>> Permute(const std::vector<std::vector<int>> ordered);
std::vector<std::vector<int>> GetStatesByDegree(const int numP, const int deg,
		const bool strict, const int min);
std::vector<std::vector<int>> GetStatesUpToDegree(const int numP, const int deg);
std::vector<std::vector<int>> GetStatesUpToDegree(const int numP, const int deg,
		const int M);
std::vector<std::vector<int>> GetStatesAtDegree(const int numP, const int deg);
std::vector<std::vector<int>> GetStatesAtDegree(const int numP, const int deg,
		const int M);

int SameQ(const poly& a, const poly& b); // if lin. dep., return a/b, else 0


#endif
