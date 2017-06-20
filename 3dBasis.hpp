#ifndef BASIS_HPP
#define BASIS_HPP

#include <exception>
#include <algorithm>
#include <vector>
#include <array>
#include <memory>
#include <iostream>
#include <gmpxx.h>

struct particle{
	int pm;	// P_-
	int pt; // P_\perp
	int pp; // P_+
};

// a mono(mial) with coefficient. It should be impossible for an instance of
// this class to be out of order, so hopefully that's true!
class mono {
	mpq_class coeff;
	std::vector<particle> particles;

	std::vector<int> IdentifyNodes() const;
	template<typename T> std::vector<int> IdentifyNodes(T (*value)(particle)) const;
	std::vector<int> IdentifyPmNodes() const;
	std::vector<int> IdentifyPtNodes() const;
	std::vector<int> IdentifyPpNodes() const;
	/*static std::vector<std::vector<int>> CfgsFromNodes(
			std::vector<int> nodes,
			std::vector<std::vector<int>> nodeEnergies);*/

	/*static std::vector<mono> BruteGenerateMinus(const int numP, const int degree);
	static void BruteAddPerp(std::vector<mono>& haveMinus, const int numP,
			const int degree);
	static void BruteAddPlus(std::vector<mono>& havePerp, const int numP,
			const int degree);*/

	public:
		mono(): coeff(1) {}
		mono(const std::vector<int>& pm, const std::vector<int>& pt,
				const std::vector<int>& pp,	const mpq_class& coeff = 1);
		mono(const std::vector<particle>& particles, const mpq_class& coeff = 1):
			coeff(coeff), particles(particles) {}

		mpq_class& Coeff()		{ return coeff; }
		const mpq_class& Coeff() const	{ return coeff; }

		unsigned int NParticles() const { return particles.size(); }

		int Pm(const int i) const 	{ return particles[i].pm; }
		int& Pm(const int i) 		{ return particles[i].pm; }
		int Pt(const int i) const 	{ return particles[i].pt; }
		int& Pt(const int i) 		{ return particles[i].pt; }
		int Pp(const int i) const 	{ return particles[i].pp; }
		int& Pp(const int i) 	 	{ return particles[i].pp; }

		mono& operator*=(const mpq_class& x)	 { coeff *= x; return *this; }
		template<typename T>
			mono& operator*=(const T& x)		 { coeff *= x; return *this; }
		mono& operator/=(const mpq_class& x)	 { coeff /= x; return *this; }
		template<typename T>
			mono& operator/=(const T& x)		 { coeff /= x; return *this; }

		bool operator==(const mono& other) const;
		bool operator!=(const mono& other) const { return !(*this == other); }
		friend std::ostream& operator<<(std::ostream& os, const mono& out);

		template<typename T>
			friend mono operator*(mono x, const T&    y) { return x *= y; }
		template<typename T>
			friend mono operator/(mono x, const T&    y) { return x /= y; }
		template<typename T>
			friend mono operator*(const T& x,    mono y) { return y *= x; }
		template<typename T>
			friend mono operator/(const T& x,    mono y) { return y /= x; }
		mono operator-() const;

		bool IsOrdered() const;
		void Order();
		mono OrderCopy() const;

		mono DerivPm(const unsigned int targetParticle) const;
		mono DerivPt(const unsigned int targetParticle) const;
		mono DerivPp(const unsigned int targetParticle) const;
		std::vector<mono> DerivPm() const;
		std::vector<mono> DerivPt() const;
		std::vector<mono> DerivPp() const;

		mono MultPm(const unsigned int targetParticle) const;
		mono MultPt(const unsigned int targetParticle) const;
		mono MultPp(const unsigned int targetParticle) const;

		std::array<mono, 4> K1(const unsigned int targetParticle) const;
		std::array<mono, 5> K2(const unsigned int targetParticle) const;
		std::array<mono, 4> K3(const unsigned int targetParticle) const;

		std::vector<mono> K1() const;
		std::vector<mono> K2() const;
		std::vector<mono> K3() const;

		/*static std::vector<std::shared_ptr<mono>> FinishPerpAtDegree(
				const std::vector<std::shared_ptr<mono>>& combined, const int degree);*/

		//static std::vector<mono> BruteGenerateMonos(const int numP, const int degree);
};

// all monos in poly should be guaranteed to be ordered correctly!
class poly {
	std::vector<mono> terms;

	public:
		explicit poly() {}
		explicit poly(const mono starter) { terms.push_back(starter); }
		explicit poly(std::vector<mono> terms);

		poly& operator+=(const mono& x);
		poly& operator-=(const mono& x);
		poly& operator+=(const poly& x);
		poly& operator-=(const poly& x);

		friend poly operator+(poly x, const poly& y);
		friend poly operator-(poly x, const poly& y);
		friend poly operator+(poly x, const mono& y);
		friend poly operator-(poly x, const mono& y);
		friend poly operator+(const mono& x, poly y);
		friend poly operator-(const mono& x, poly y);

		const	mono& operator[](size_t i) const	{ return terms[i]; }
				mono& operator[](size_t i)			{ return terms[i]; }
		std::vector<mono>::const_iterator	begin() const noexcept
				{ return terms.begin(); }
		std::vector<mono>::iterator			begin() noexcept
				{ return terms.begin(); }
		std::vector<mono>::const_iterator	end()	const noexcept
				{ return terms.end(); }
		std::vector<mono>::iterator			end()	noexcept
				{ return terms.end(); }

		poly DerivPm(const int);
		poly DerivPp(const int);
		poly DerivPm();
		poly DerivPp();

		poly K1() const;
		poly K2() const;
		poly K3() const;

		poly L1() const;
		poly L2() const;
		poly L3() const;
};

// a set of all ordered monomials of the given degree and number of particles
class basis {
	int degree;
	int numP;
	std::vector<std::unique_ptr<mono>> basisMonos;

	/*static std::vector<std::shared_ptr<mono>> AddPlusUpToDegree(
			const std::vector<std::vector<int>>& minus, const int degree);*/
	/*static std::vector<std::shared_ptr<mono>> AddPlusUpToDegree(
			const std::vector<std::vector<int>>& minus, const int degree,
			const int M);*/

	static std::vector<std::vector<particle>> CombinedCfgs(
			const std::vector<particle>& baseCfg,
			const std::vector<std::vector<int>>& newCfgs,
			const int componentToChange);
	static std::vector<std::vector<int>> CfgsFromNodes(const int remainingEnergy,
			const std::vector<int>& nodes, const bool exact);
	static std::vector<std::vector<int>> CfgsFromNodePartition(
			const std::vector<int> nodes,
			std::vector<std::vector<int>> nodeEnergies);
	public:
		basis() = delete;
		basis(const int numP, const int degree);
		//basis(const int numP, const int degree, const int M);
		mono* GetBasisMono(const std::vector<int>& pm, const std::vector<int>&,
				const std::vector<int>& pp);
		mono* GetBasisMono(const mono& wildMono);

		std::vector<std::unique_ptr<mono>>& BasisMonos() { return basisMonos; }
		const std::vector<std::unique_ptr<mono>>& BasisMonos() const { return basisMonos; }

		std::vector<mpq_class> ExpressPoly(const poly& toExpress) const;
};

std::ostream& operator<<(std::ostream& os, const std::vector<int>& out);
std::ostream& operator<<(std::ostream& os, const std::vector<mpq_class>& out);
std::ostream& operator<<(std::ostream& os, const std::vector<particle>& out);
std::vector<std::vector<int>> Permute(const std::vector<std::vector<int>> ordered);
std::vector<std::vector<int>> GetStatesByDegree(const int numP, const int deg,
		const bool exact, const int min);
std::vector<std::vector<int>> GetStatesUpToDegree(const int numP, const int deg);
/*std::vector<std::vector<int>> GetStatesUpToDegree(const int numP, const int deg,
		const int M);*/
std::vector<std::vector<int>> GetStatesAtDegree(const int numP, const int deg);
/*std::vector<std::vector<int>> GetStatesAtDegree(const int numP, const int deg,
		const int M);*/

inline std::ostream& operator<<(std::ostream &o, const mpq_class &expr){
	return o << expr.get_str();
}

// T can be any type or class with an == operator; value indexes T by uint
template<typename Accessor>
//inline std::vector<int> IdentifyNodes(T (*value)(unsigned int), const size_t size){
inline std::vector<int> IdentifyNodes(Accessor A, const size_t size){
	std::vector<int> nodes;
	int newNode;
	newNode = 1;
	for(unsigned int i = 1; i < size; ++i){
		if(A(i) != A(i-1)){
			nodes.push_back(newNode);
			newNode = 1;
		} else {
			++newNode;
		}
	}
	nodes.push_back(newNode);
	return nodes;
}

// calls the generic IdentifyNodes using the class's particles and (*value)
template<typename T>
inline std::vector<int> mono::IdentifyNodes(T (*value)(particle)) const{
	return IdentifyNodes([this, value](unsigned int i){return value(this->particles[i]);},
			NParticles());
}

// should work for any containerish thing T with operator[]
template<typename T>
inline std::vector<int> IdentifyNodes(const T& container){
	return IdentifyNodes([container](unsigned int i){return container[i];},
			container.size());
}

// default behavior for a vector of particles including mono::IdentifyNodes()
template<>
inline std::vector<int> IdentifyNodes(const std::vector<particle>& particles){
	return IdentifyNodes([particles](unsigned int i){return 
			std::array<int, 3>({{particles[i].pm, particles[i].pt, particles[i].pp}}); },
			particles.size());
}

/*template<typename T, typename U>
inline std::vector<int> IdentifyNodes(T& container, U (*value)(
*/
#endif
