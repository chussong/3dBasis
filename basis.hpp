#ifndef BASIS_HPP
#define BASIS_HPP

#include <vector>
#include "mono.hpp"
#include "poly.hpp"

// a set of all ordered monomials of the given degree and number of particles
// * Can be accessed with [] like an array and iterated through with begin()
// and end().
// * Has a stream operator, which outputs all the component monos on new lines.
// * Most important function is ExpressPoly, which takes a polynomial and
// expresses it as a vector in this basis's space.
class basis {
	std::vector<mono> basisMonos;

	static std::vector<std::vector<int>> CfgsFromNodePartition(
			const std::vector<int> nodes,
			std::vector<std::vector<int>> nodeEnergies);
	static bool PpDominant(const mono& m);

	public:
		explicit basis(const std::vector<mono> basisMonos): basisMonos(basisMonos) {}
		basis(const int numP, const int degree, const int options);
		unsigned int FindInBasis(const std::vector<int>& pm,
				const std::vector<int>& pt, const std::vector<int>& pp);
		unsigned int FindInBasis(const mono& wildMono);

		//std::array<basis,2> ParitySplit() const;
		void DeleteOdd();
		void DeleteEven();
		void DeleteSymm();
		void DeleteAsymm();
		//void SortBasis();

		const	mono& operator[](size_t i)	const	{return basisMonos[i];     }
		friend std::ostream& operator<<(std::ostream& os, const basis& out);

		std::size_t size()					const	{return basisMonos.size(); }
		std::vector<mono>::const_iterator	begin() const noexcept
				{ return basisMonos.begin(); }
		std::vector<mono>::iterator			begin() noexcept
				{ return basisMonos.begin(); }
		std::vector<mono>::const_iterator	end()	const noexcept
				{ return basisMonos.end(); }
		std::vector<mono>::iterator			end()	noexcept
				{ return basisMonos.end(); }

		Triplet ExpressMono(const mono& toExpress, const int column,
				const int rowOffset) const;
		std::list<Triplet> ExpressPoly(const poly& toExpress, const int column, 
				const int rowOffset) const;

		static std::vector<std::vector<int>> CfgsFromNodes(const int remainingEnergy,
				const std::vector<int>& nodes, const bool exact);
		static std::vector<std::vector<particle>> CombinedCfgs(
				const std::vector<particle>& baseCfg,
				const std::vector<std::vector<int>>& newCfgs,
				const int componentToChange);
};

// Contains two bases and intelligently decides which one to use for various
// requests. Because this is in two pieces, it can not be iterated through.
class splitBasis {
	basis evenBasis;
	basis oddBasis;

	public:
		splitBasis(const int numP, const int degree, const int options);

		std::pair<unsigned int, basis*> FindInBasis(const std::vector<int>& pm,
				const std::vector<int>& pt, const std::vector<int>& pp);
		std::pair<unsigned int, basis*> FindInBasis(const mono& wildMono);

		basis& OddBasis() { return oddBasis; }
		const basis& OddBasis() const { return oddBasis; }
		basis& EvenBasis() { return evenBasis; }
		const basis& EvenBasis() const { return evenBasis; }

		splitBasis BecomeAsymmetric();
		void DeleteSymm();
		void DeleteAsymm();
		//splitPolyBasis AdditiveSymmetrization();

		std::list<Triplet> ExpressPoly(const poly& toExpress, const int column,
				const int row) const;

		static bool IsOdd (const mono& toTest);
		static bool IsEven(const mono& toTest);
		static bool IsSymm(const mono& toTest);
		static bool IsAsymm(const mono& toTest);
};

// Basis made of 2D+1 sub-bases stratified by M
class mBasis {
	std::vector<std::vector<poly>> multiplets;

	public:
		mBasis(const int numP, const int degree, const int options);

		basis BasisAtM(const int numP, const int degree, const int M,
				const int options);

		Matrix L3Matrix(const basis& startingMBasis, const basis& targetMBasis);
		std::vector<poly> CompleteMultiplet(const poly& topState);
};

// state generation -----------------------------------------------------------

std::vector<std::vector<int>> Permute(const std::vector<std::vector<int>> ordered);
std::vector<std::vector<int>> GetStatesByDegree(const int numP, const int deg,
		const bool exact, const int min);
std::vector<std::vector<int>> GetStatesUpToDegree(const int numP, const int deg,
		const int M = 0);
std::vector<std::vector<int>> GetStatesAtDegree(const int numP, const int deg,
		const int M = 0);
bool EoMAllowed(const std::vector<particle>& cfg);

#endif
