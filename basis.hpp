#ifndef BASIS_HPP
#define BASIS_HPP

#include <vector>
#include "constants.hpp"
#include "construction.hpp"
#include "mono.hpp"
#include "poly.hpp"

// a basis is anything which claims that it can express monomials and
// polynomials as vectors on itself. It can contain monos or polys as needed.
template<class T>
class Basis {
	std::vector<T> basisVectors;

	/*static std::vector<std::vector<int>> CfgsFromNodePartition(
			const std::vector<int> nodes,
			std::vector<std::vector<int>> nodeEnergies);*/
	static bool PpDominant(const mono& m);

	public:
		explicit Basis(const std::vector<T> basisVectors): basisVectors(basisVectors) {}
		Basis(const int numP, const int degree, const int options);

		unsigned int FindInBasis(const std::vector<int>& pm,
				const std::vector<int>& pt, const std::vector<int>& pp) const;
		unsigned int FindInBasis(const T& wildVector) const;

		void DeleteOdd();
		void DeleteEven();
		void DeleteSymm();
		void DeleteAsymm();

		//std::array<basis,2> ParitySplit() const;
		//virtual void DeleteOdd();
		//virtual void DeleteEven();
		//virtual void DeleteSymm();
		//virtual void DeleteAsymm();
		//void SortBasis();

		const T& operator[](size_t i) const;
		friend std::ostream& operator<<(std::ostream& os, const Basis<T>& out);

		std::size_t size() const;
		typename std::vector<T>::const_iterator		begin() const noexcept
				{ return basisVectors.begin(); }
		typename std::vector<T>::iterator			begin() noexcept
				{ return basisVectors.begin(); }
		typename std::vector<T>::const_iterator		end()	const noexcept
				{ return basisVectors.end(); }
		typename std::vector<T>::iterator			end()	noexcept
				{ return basisVectors.end(); }

		Triplet ExpressMono(const mono& toExpress, const int column,
				const int rowOffset) const;
		std::list<Triplet> ExpressPoly(const poly& toExpress, 
				const int column, const int rowOffset) const;

		/*std::vector<std::vector<int>> CfgsFromNodes(const int remainingEnergy,
				const std::vector<int>& nodes, const bool exact);
		std::vector<std::vector<particle>> CombinedCfgs(
				const std::vector<particle>& baseCfg,
				const std::vector<std::vector<int>>& newCfgs,
				const int componentToChange);*/
};

// a set of all ordered monomials of the given degree and number of particles
// * Can be accessed with [] like an array and iterated through with begin()
// and end().
// * Has a stream operator, which outputs all the component monos on new lines.
// * Most important function is ExpressPoly, which takes a polynomial and
// expresses it as a vector in this basis's space.
/*template<>
class monoBasis : public basis {
	std::vector<mono> basisMonos;

	public:
		explicit monoBasis(const std::vector<mono> basisMonos): basisMonos(basisMonos) {}
		monoBasis(const int numP, const int degree, const int options);

		unsigned int FindInBasis(const std::vector<int>& pm, 
				const std::vector<int>& pt, const std::vector<int>& pp) const;
		unsigned int FindInBasis(const mono& wildMono) const;
		void DeleteOdd();
		void DeleteEven();
		void DeleteSymm();
		void DeleteAsymm();

		const	mono& operator[](size_t i)	const	{return basisMonos[i];     }
		friend std::ostream& operator<<(std::ostream& os, const monoBasis& out);

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
		std::list<Triplet> ExpressPoly(const poly& toExpress, 
				const int column, const int rowOffset) const;
};

class polyBasis : public basis {
	std::vector<poly> basisPolys;

	public:
		explicit polyBasis(const std::vector<poly> basisPolys): basisPolys(basisPolys) {}
		polyBasis(const int numP, const int degree, const int options);

		unsigned int FindInBasis(const std::vector<int>& pm, 
				const std::vector<int>& pt, const std::vector<int>& pp) const;
		unsigned int FindInBasis(const mono& wildMono) const;
		unsigned int FindInBasis(const poly& wildPoly) const;
		void DeleteOdd();
		void DeleteEven();
		void DeleteSymm();
		void DeleteAsymm();

		const	poly& operator[](size_t i)	const	{return basisPolys[i];     }
		friend std::ostream& operator<<(std::ostream& os, const polyBasis& out);

		std::size_t size()					const	{return basisPolys.size(); }
		std::vector<poly>::const_iterator	begin() const noexcept
				{ return basisPolys.begin(); }
		std::vector<poly>::iterator			begin() noexcept
				{ return basisPolys.begin(); }
		std::vector<poly>::const_iterator	end()	const noexcept
				{ return basisPolys.end(); }
		std::vector<poly>::iterator			end()	noexcept
				{ return basisPolys.end(); }


		Triplet ExpressMono(const mono& toExpress, const int column,
				const int rowOffset) const;
		std::list<Triplet> ExpressPoly(const poly& toExpress, 
				const int column, const int rowOffset) const;
};*/

// Contains two bases and intelligently decides which one to use for various
// requests. Because this is in two pieces, it can not be iterated through.
template<class T>
class splitBasis {
	Basis<T> evenBasis;
	Basis<T> oddBasis;

	public:
		splitBasis(const int numP, const int degree, const int options);

		std::pair<unsigned int, Basis<T>*> FindInBasis(const std::vector<int>& pm,
				const std::vector<int>& pt, const std::vector<int>& pp);
		std::pair<unsigned int, Basis<T>*> FindInBasis(const mono& wildMono);

		Basis<T>& OddBasis() { return oddBasis; }
		const Basis<T>& OddBasis() const { return oddBasis; }
		Basis<T>& EvenBasis() { return evenBasis; }
		const Basis<T>& EvenBasis() const { return evenBasis; }

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
		static std::pair<poly, poly> OddEvenSplit  (const poly& toSplit);
		static std::pair<poly, poly> SymmAsymmSplit(const poly& toSplit);
};

// Basis made of 2D+1 sub-bases stratified by M
class mBasis {
	std::vector<Basis<mono>> mLevels;
	std::vector<std::vector<poly>> topStates; // L eigenstates with M = L = i

	public:
		mBasis(const int numP, const int degree, const int options);

		Basis<mono> BasisAtM(const int numP, const int degree, const int M,
				const int options);

		Matrix L3Matrix(const Basis<mono>& startingMBasis, const Basis<mono>& targetMBasis);
		std::vector<poly> CompleteMultiplet(const poly& topState);

		const Basis<mono>* Level(const int M) const;

		unsigned int Degree() const { return mLevels.size()-1; }
		unsigned int LevelSize(const int L) const;
		std::vector<poly>& TopStates(const int L) { return topStates[L]; }
		const std::vector<poly>& TopStates(const int L) const { return topStates[L]; }
};

// state generation -----------------------------------------------------------

/*std::vector<std::vector<int>> Permute(const std::vector<std::vector<int>> ordered);
std::vector<std::vector<int>> GetStatesByDegree(const int numP, const int deg,
		const bool exact, const int min);
std::vector<std::vector<int>> GetStatesUpToDegree(const int numP, const int deg,
		const int M = 0);
std::vector<std::vector<int>> GetStatesAtDegree(const int numP, const int deg,
		const int M = 0);*/
bool EoMAllowed(const std::vector<particle>& cfg);

template<class T>
inline Basis<T>::Basis(const int, const int, const int) {
	std::cerr << "Error: ordered to construct a basis by degree for an "
		<< "underlying type which has not been specialized. Please construct "
		<< "with a different type or write a specialization for this type."
		<< std::endl;
}

template<>
Basis<mono>::Basis(const int numP, const int degree, const int options) {
	// 1: generate all possibilities for P_-
	// 2: identify nodes in P_-, generate possible distributions of P_\perp to
	// 		the nodes and then P_\perp within each node, adding each at top level
	// 3: identify nodes in (P_- and P_\perp), repeating step 2 for P_+
	// 4: take list from step 3 and create a mono from each entry, then store
	// 		the list of these as basisMonos
	// NOTE: this does not make any use of the symmetries between the components
	// 		in constructing the basis; for instance, the states where every
	// 		pm = 0 are followed by a copy of the earlier parts of the basis
	const bool useEoM = (options & OPT_EQNMOTION) != 0;
	std::vector<std::vector<int>> minus = GetStatesUpToDegree(numP, degree);
	std::vector<std::vector<particle>> particleCfgs;
	for(auto& minusCfg : minus){
		std::vector<particle> newCfg(minusCfg.size());
		for(auto i = 0u; i < newCfg.size(); ++i) newCfg[i].pm = minusCfg[i];
		particleCfgs.push_back(newCfg);
	}

	std::vector<int> nodes;
	std::vector<std::vector<particle>> newCfgs;
	for(auto& configuration : particleCfgs){
		nodes = IdentifyNodes(configuration);
		int remainingEnergy = degree;
		for(auto& part : configuration) remainingEnergy -= part.pm;
		std::vector<std::vector<int>> perp(CfgsFromNodes(remainingEnergy, nodes,
															false));
		for(auto& newCfg : CombinedCfgs(configuration, perp, 2)){
			newCfgs.push_back(newCfg);
		}
	}

	particleCfgs.clear();
	for(auto& cfg : newCfgs){
		nodes = IdentifyNodes(cfg);
		int remainingEnergy = degree;
		for(auto& part : cfg) remainingEnergy -= part.pm + part.pt;
		std::vector<std::vector<int>> plus(CfgsFromNodes(remainingEnergy, nodes,
															true));
		for(auto& newCfg : CombinedCfgs(cfg, plus, 3)){
			if(!useEoM || EoMAllowed(newCfg)) particleCfgs.push_back(newCfg);
		}
	}

	for(auto& cfg : particleCfgs){
		basisVectors.emplace_back(cfg, useEoM);
		//std::cout << cfg << std::endl;
	}
}

template<class T>
inline unsigned int Basis<T>::FindInBasis(const std::vector<int>&,
		const std::vector<int>&, const std::vector<int>&) const{
	std::cout << "Unspecialized Basis::FindInBasis(vectors) should not be called. "
		<< "How was this basis object created in the first place?" << std::endl;
	return -1u;
}

template<class T>
inline Triplet Basis<T>::ExpressMono(const mono&, const int, const int) const{
	std::cout << "Unspecialized Basis::ExpressMono should never be called. "
		<< "How was a basis object created in the first place?" << std::endl;
	return Triplet(-1, -1, coeff_class(0));
}

template<class T>
inline std::list<Triplet> Basis<T>::ExpressPoly(const poly&, const int, const int) const{
	std::cout << "basis::ExpressMono should never be called. "
		<< "How was a basis object created in the first place?" << std::endl;
	return {{Triplet(-1, -1, coeff_class(0))}};
}

/*template<>
inline unsigned int Basis<mono>::FindInBasis(const mono& wildMono) const{
	for(auto i = 0u; i < basisVectors.size(); ++i){
		if(basisVectors[i] == wildMono) return i;
	}
	std::cout << "Warning! Failed to find the following mono in our basis: "
		<< wildMono << std::endl;
	return -1u;
}*/

template<class T>
inline unsigned int Basis<T>::FindInBasis(const T& wildVector) const{
	for(auto i = 0u; i < basisVectors.size(); ++i){
		if(basisVectors[i] == wildVector) return i;
	}
	std::cout << "Warning! Failed to find the following vector in our basis: "
		<< wildVector << std::endl;
	return -1u;
}

template<>
inline unsigned int Basis<mono>::FindInBasis(const std::vector<int>& pm, 
		const std::vector<int>& pt, const std::vector<int>& pp) const{
	return FindInBasis(mono(pm, pt, pp));
}

/*std::array<basis,2> Basis<mono>::ParitySplit() const{
	std::array<basis,2> ret = {{*this, *this}};
	for(int i = 0; i < 2; ++i){
		ret[i].basisVectors.erase(std::remove_if(ret[i].begin(), ret[i].end(), 
				[i](mono& m){ return m.TotalPt()%2 == i; } ));
	}
	return ret;
}*/

template<>
inline void Basis<mono>::DeleteOdd(){
	basisVectors.erase(std::remove_if(basisVectors.begin(), basisVectors.end(), 
				splitBasis<mono>::IsOdd), basisVectors.end());
}

template<>
inline void Basis<mono>::DeleteEven(){
	basisVectors.erase(std::remove_if(basisVectors.begin(), basisVectors.end(),
				splitBasis<mono>::IsEven), basisVectors.end());
}

// this deletes the +/- symmetric elements AND deletes the P+ dominant ones,
// defined to be those with P+ > P- or max(P+) > max(P-)
template<>
inline void Basis<mono>::DeleteSymm(){
	basisVectors.erase(std::remove_if(basisVectors.begin(), basisVectors.end(),
				splitBasis<mono>::IsSymm), basisVectors.end());
	basisVectors.erase(std::remove_if(basisVectors.begin(), basisVectors.end(),
				mono::MirrorIsBetter), basisVectors.end());
	std::cout << "Basis after deleting symmetric and unfavorable elements: \n";
	for(auto& m : basisVectors) std::cout << m << std::endl;
	std::cout << "--------------------" << std::endl;
}

template<>
inline void Basis<mono>::DeleteAsymm(){
	basisVectors.erase(std::remove_if(basisVectors.begin(), basisVectors.end(),
				splitBasis<mono>::IsAsymm), basisVectors.end());
}

/*void Basis<mono>::SortBasis(){
	std::sort(basisVectors.begin(), basisVectors.end(), mono::MonoPrecedence);
	std::cout << "Sorted basis to this:" << std::endl;
	for(auto i = 0u; i < basisVectors.size()/2; ++i){
		std::cout << basisVectors[i] << std::endl;
	}
	std::cout << "----- half way pt -----" << std::endl;
	for(auto i = basisVectors.size()/2; i < basisVectors.size(); ++i){
		std::cout << basisVectors[i] << std::endl;
	}
	std::cout << "-------------------" << std::endl;
}*/

template<class T>
inline std::ostream& operator<<(std::ostream& os, const Basis<T>& out){
	os << "{ ";
	for(auto& v : out.basisVectors) os << v.HumanReadable() << ", ";
	return os << "\b\b }";
}

/*template<>
inline std::ostream& operator<<(std::ostream& os, const Basis<mono>& out){
	os << "{ ";
	for(auto& m : out.basisVectors) os << m.HumanReadable() << ", ";
	return os << "\b\b }";
}*/

template<>
inline Triplet Basis<mono>::ExpressMono(const mono& toExpress, const int column,
		const int rowOffset) const{
	for(auto i = 0u; i < basisVectors.size(); ++i){
		if(toExpress == basisVectors[i])
			return Triplet(rowOffset+i, column, toExpress.Coeff());
	}
	std::cerr << "Error: tried to express the monomial " << toExpress
	<< " on the given basis but was not able to identify it." << std::endl;
	return Triplet(-1, -1, toExpress.Coeff());
}

template<>
inline std::list<Triplet> Basis<mono>::ExpressPoly(const poly& toExpress, 
		const int column, const int rowOffset) const{
	std::list<Triplet> ret;
	if(toExpress.size() == 0) return ret;
	unsigned int hits = 0u;
	for(auto i = 0u; i < basisVectors.size(); ++i){
		for(auto& term : toExpress){
			if(term == basisVectors[i]){
				ret.emplace_front(rowOffset+i, column, term.Coeff());
				++hits;
				break;
			}
		}
		if(hits == toExpress.size()) break;
	}
	if(hits < toExpress.size()){
		std::cerr << "Error: tried to express the polynomial " << toExpress
		<< " on the given basis but was only able to identify " << hits
		<< " of " << toExpress.size() << " terms." << std::endl;
	}
	//std::cout << "Expressed " << toExpress << " as the following triplets:" << std::endl;
	//for(auto& trip : ret) std::cout << trip << std::endl;
	return ret;
}

/*inline unsigned int Basis<poly>::FindInBasis(const poly& wildPoly) const{
	for(auto i = 0u; i < basisVectors.size(); ++i){
		if(basisVectors[i] == wildPoly) return i;
	}
	std::cout << "Warning! Failed to find the following poly in our basis: "
		<< wildPoly << std::endl;
	return -1u;
}

template<>
inline unsigned int Basis<poly>::FindInBasis(const mono& wildMono) const{
	return FindInBasis(poly(wildMono));
}

template<>
inline unsigned int Basis<poly>::FindInBasis(const std::vector<int>& pm, 
		const std::vector<int>& pt, const std::vector<int>& pp) const{
	return FindInBasis(mono(pm, pt, pp));
}*/

// This does not attempt to find a linear combination of known polynomials
// which would reproduce toExpress. Obviously it would be more correct if it
// did attempt to do so.
template<>
inline Triplet Basis<poly>::ExpressMono(const mono& toExpress, const int column,
		const int rowOffset) const{
	poly polyForm(toExpress);
	for(auto i = 0u; i < basisVectors.size(); ++i){
		if(polyForm == basisVectors[i])
			return Triplet(rowOffset+i, column, toExpress.Coeff());
	}
	std::cerr << "Error: tried to express the monomial " << toExpress
	<< " on the given basis but was not able to identify it." << std::endl;
	return Triplet(-1, -1, toExpress.Coeff());
}

// Again, no linear combinations are checked -- we just look to see if we have
// the poly as written. This is obviously not ideal. It also means we can only
// get exactly 1*the polynomial, which actually probably breaks stuff.
template<>
inline std::list<Triplet> Basis<poly>::ExpressPoly(const poly& toExpress, 
		const int column, const int rowOffset) const{
	for(auto i = 0u; i < basisVectors.size(); ++i){
		if(toExpress == basisVectors[i])
			return {{Triplet(rowOffset+i, column, 1)}};
	}
	std::cerr << "Error: tried to express the monomial " << toExpress
	<< " on the given basis but was not able to identify it." << std::endl;
	return {{Triplet(-1, -1, 0)}};
}

template<class T>
inline bool splitBasis<T>::IsOdd(const mono& toTest){
	return toTest.TotalPt()%2 == 1;
}

template<class T>
inline bool splitBasis<T>::IsEven(const mono& toTest){
	return !IsOdd(toTest);
}

template<class T>
inline bool splitBasis<T>::IsSymm(const mono& toTest){
	mono clone(toTest);
	clone.MirrorPM();
	/*std::cout << toTest;
	std::cout << (clone == toTest ? " == " : " != ");
	std::cout << clone << std::endl;*/
	return clone == toTest;
}

template<class T>
inline bool splitBasis<T>::IsAsymm(const mono& toTest){
	return !IsSymm(toTest);
}

// this could definitely be done more intelligently if speed were important
template<class T>
inline splitBasis<T>::splitBasis(const int numP, const int degree, const int options): 
	evenBasis(numP, degree, options), oddBasis(numP, degree, options){
	evenBasis.DeleteOdd();
	oddBasis.DeleteEven();
}

template<class T>
inline splitBasis<T> splitBasis<T>::BecomeAsymmetric(){
	splitBasis symm(*this);
	DeleteSymm();
	symm.DeleteAsymm();
	return symm;
}

template<class T>
inline void splitBasis<T>::DeleteSymm(){
	oddBasis.DeleteSymm();
	evenBasis.DeleteSymm();
}

template<class T>
inline void splitBasis<T>::DeleteAsymm(){
	oddBasis.DeleteAsymm();
	evenBasis.DeleteAsymm();
}

// create a copy of this basis which is symmetrized by adding mirror operators
// to all member monomials
//splitPolyBasis AdditiveSymmetrization(){
//}

template<class T>
inline std::list<Triplet> splitBasis<T>::ExpressPoly(const poly& toExpress, 
		const int column, const int rowOffset) const{
	std::list<Triplet> ret;
	for(auto& term : toExpress){
		if(IsOdd(term)){
			ret.push_front(oddBasis.ExpressMono(term, column, rowOffset));
		} else {
			ret.push_front(evenBasis.ExpressMono(term, column, rowOffset));
		}
	}
	return ret;
}

#endif
