#ifndef MONO_HPP
#define MONO_HPP

#include <iostream>
#include <array>
#include <vector>
#include <ostream>
#include <cmath>		// lgamma
#include <algorithm>	// next_permutation

#include "constants.hpp"
#include "construction.hpp"

// a mono(mial) with coefficient. It should be impossible for an instance of
// this class to be out of order, so hopefully that's true! This class in its
// current form ONLY represents Dirichlet monomials using the equations of 
// motion.
// Some notes:
// * You CAN NOT add or subtract monos from each other. Instead, add or subtract
// them from a poly.
// * You CAN multiply or divide a mono by anything which can multiply or divide
// coeff_class. If you do, coeff is changed while the momenta are unaffected.
// * To multiply or divide by momenta, use MultPm(i), etc. These return a new
// mono, however; momenta of an existing mono can only be changed by directly
// accessing them through Pm(i), etc.
// * Two monos are equal (i.e. a == b is true) if their MOMENTA are the same.
// The coefficients do not matter in comparisons, only the momenta.
// * std::cout << someMono; will print the mono in this format:
// coeff * {p1_m, p2_m, ... }{p1_t, p2_t, ...}. You may prefer
// someMono.HumanReadable(), which more resembles how you'd write it on a board.
class Mono {
	coeff_class coeff;
	std::vector<particle> particles;

	std::vector<int> IdentifyNodes() const;
	template<typename T> std::vector<int> IdentifyNodes(T (*value)(particle)) const;
	std::vector<int> IdentifyPmNodes() const;
	std::vector<int> IdentifyPtNodes() const;

	public:
		Mono(): coeff(1) {}
		Mono(const std::vector<int>& pm, const std::vector<int>& pt, 
				const coeff_class& coeff = 1);
		Mono(const std::vector<particle>& particles, 
				const coeff_class& coeff = 1);

		coeff_class& Coeff()		{ return coeff; }
		const coeff_class& Coeff() const	{ return coeff; }

		unsigned int NParticles() const { return particles.size(); }

		const char& Pm(const int i) const;
		//char& Pm(const int i);
		void ChangePm(const int i, const char newValue);
		const char& Pt(const int i) const;
		//char& Pt(const int i);
		void ChangePt(const int i, const char newValue);
		int TotalPm() const;
		int TotalPt() const;
		int MaxPm() const;
		int MaxPt() const;
		int Degree() const { return TotalPm() + TotalPt(); }
		std::vector<size_t> CountIdentical() const;
		std::vector<size_t> PermutationVector() const;

		Mono& operator*=(const coeff_class& x)	 { coeff *= x; return *this; }
		template<typename T>
			Mono& operator*=(const T& x)		 { coeff *= x; return *this; }
		Mono& operator/=(const coeff_class& x)	 { coeff /= x; return *this; }
		template<typename T>
			Mono& operator/=(const T& x)		 { coeff /= x; return *this; }

		bool operator==(const Mono& other) const;
		bool operator!=(const Mono& other) const { return !(*this == other); }
		friend std::ostream& operator<<(std::ostream& os, const Mono& out);
        friend std::string MathematicaOutput(const Mono& out);
		std::string HumanReadable() const;

		template<typename T>
			friend Mono operator*(Mono x, const T&    y) { return x *= y; }
		template<typename T>
			friend Mono operator/(Mono x, const T&    y) { return x /= y; }
		template<typename T>
			friend Mono operator*(const T& x,    Mono y) { return y *= x; }
		Mono operator-() const;

		static bool ParticlePrecedence(const particle& a, const particle& b);
		bool IsOrdered() const;
		void Order();
		Mono OrderCopy() const;

		bool IsNull() const;
		void NullIfIllegal();

		Mono DerivPm(const unsigned int targetParticle) const;
		Mono DerivPt(const unsigned int targetParticle) const;
		Mono DerivPp(const unsigned int targetParticle) const;
		std::vector<Mono> DerivPm() const;
		std::vector<Mono> DerivPt() const;
		std::vector<Mono> DerivPp() const;

		Mono MultPm(const unsigned int targetParticle) const;
		Mono MultPt(const unsigned int targetParticle) const;
		Mono MultPp(const unsigned int targetParticle) const;

		// static coeff_class InnerProduct(const Mono& A, const Mono& B,
										// const GammaCache& cache,
										// const KVectorCache& kCache);
// 
		// static coeff_class IPZuhair(const Mono& A, const Mono& B,
									// const GammaCache& cache,
									// const KVectorCache& kCache);
		// static std::vector<std::vector<char>> VectorsAtK(const char totalK, 
				// const std::vector<size_t>& perm, const Mono& A, const Mono& B,
				// const size_t start);
};

// calls the generic IdentifyNodes using the class's particles and (*value)
template<typename T>
inline std::vector<int> Mono::IdentifyNodes(T (*value)(particle)) const{
	return IdentifyNodes([this, value](unsigned int i){return value(this->particles[i]);},
			NParticles());
}

// specialization of vector ostream template
template<>
inline std::ostream& operator<<(std::ostream& os, const std::vector<Mono>& out){
	os << "{ ";
	int written = 0;
	for(auto& element : out){
		if(element.Coeff() == 0) continue;
		if(element.Coeff() > 0) os << " ";
		os << element.HumanReadable() << ", ";
		++written;
	}
	if(written > 0) os << "\b\b ";
	return os << "\b }";
}

// partial specialization of array ostream template
template<int N>
inline std::ostream& operator<<(std::ostream& os, const std::array<Mono,N>& out){
	os << "{ ";
	for(auto& element : out){
		if(element.Coeff() > 0) os << " ";
		if(element.Coeff() != 0) os << element.HumanReadable() << ", ";
	}
	return os << "\b }";
}

#endif
