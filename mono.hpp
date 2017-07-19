#ifndef MONO_HPP
#define MONO_HPP

#include <iostream>
#include <array>
#include <vector>
#include <ostream>

#include "constants.hpp"
#include "construction.hpp"

// a mono(mial) with coefficient. It should be impossible for an instance of
// this class to be out of order, so hopefully that's true!
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
// coeff * {p1_m, p2_m, ... }{p1_t, p2_t, ...}{p1_p, p2_p, ...}. You may prefer
// someMono.HumanReadable(), which more resembles how you'd write it on a board.
class mono {
	coeff_class coeff;
	std::vector<particle> particles;
	bool usingEoM;

	std::vector<int> IdentifyNodes() const;
	template<typename T> std::vector<int> IdentifyNodes(T (*value)(particle)) const;
	std::vector<int> IdentifyPmNodes() const;
	std::vector<int> IdentifyPtNodes() const;
	std::vector<int> IdentifyPpNodes() const;

	public:
		mono(const bool usingEoM = false): coeff(1), usingEoM(usingEoM) {}
		mono(const std::vector<int>& pm, const std::vector<int>& pt,
				const std::vector<int>& pp,	const bool usingEoM = false,
				const coeff_class& coeff = 1);
		mono(const std::vector<particle>& particles, const bool usingEoM = false,
				const coeff_class& coeff = 1);

		coeff_class& Coeff()		{ return coeff; }
		const coeff_class& Coeff() const	{ return coeff; }

		unsigned int NParticles() const { return particles.size(); }
		bool UsingEoM() const { return usingEoM; }

		int Pm(const int i) const 	{ return particles[i].pm; }
		int& Pm(const int i) 		{ return particles[i].pm; }
		int Pt(const int i) const 	{ return particles[i].pt; }
		int& Pt(const int i) 		{ return particles[i].pt; }
		int Pp(const int i) const 	{ return particles[i].pp; }
		int& Pp(const int i) 	 	{ return particles[i].pp; }
		int TotalPm() const;
		int TotalPt() const;
		int TotalPp() const;
		int MaxPm() const;
		int MaxPt() const;
		int MaxPp() const;
		int Degree() const { return TotalPm() + TotalPt() + TotalPp(); }

		mono& operator*=(const coeff_class& x)	 { coeff *= x; return *this; }
		template<typename T>
			mono& operator*=(const T& x)		 { coeff *= x; return *this; }
		mono& operator/=(const coeff_class& x)	 { coeff /= x; return *this; }
		template<typename T>
			mono& operator/=(const T& x)		 { coeff /= x; return *this; }

		bool operator==(const mono& other) const;
		bool operator!=(const mono& other) const { return !(*this == other); }
		friend std::ostream& operator<<(std::ostream& os, const mono& out);
		std::string HumanReadable() const;

		template<typename T>
			friend mono operator*(mono x, const T&    y) { return x *= y; }
		template<typename T>
			friend mono operator/(mono x, const T&    y) { return x /= y; }
		template<typename T>
			friend mono operator*(const T& x,    mono y) { return y *= x; }
		mono operator-() const;

		static bool ParticlePrecedence(const particle& a, const particle& b);
		bool IsOrdered() const;
		void Order();
		mono OrderCopy() const;

		mono MirrorPM() const;
		static bool MirrorIsBetter(const mono& m);

		void NullIfIllegal();
		mono DerivPm(const unsigned int targetParticle) const;
		mono DerivPt(const unsigned int targetParticle) const;
		mono DerivPp(const unsigned int targetParticle) const;
		std::vector<mono> DerivPm() const;
		std::vector<mono> DerivPt() const;
		std::vector<mono> DerivPp() const;

		mono MultPm(const unsigned int targetParticle) const;
		mono MultPt(const unsigned int targetParticle) const;
		mono MultPp(const unsigned int targetParticle) const;

		// K on specific particle WITHOUT EoM
		std::array<mono, 4> K1(const unsigned int targetParticle,
				const coeff_class delta) const;
		std::array<mono, 5> K2(const unsigned int targetParticle,
				const coeff_class delta) const;
		std::array<mono, 4> K3(const unsigned int targetParticle,
				const coeff_class delta) const;

		// K on specific particle WITH EoM
		std::array<mono, 1> K1_EoM(const unsigned int targetParticle,
				const coeff_class delta) const;
		std::array<mono, 3> K2_EoM(const unsigned int targetParticle,
				const coeff_class delta) const;
		std::array<mono, 4> K3_EoM(const unsigned int targetParticle,
				const coeff_class delta) const;

		// K on entire mono
		std::vector<mono> K1(const coeff_class delta) const;
		std::vector<mono> K2(const coeff_class delta) const;
		std::vector<mono> K3(const coeff_class delta) const;

		// L on specific particle
		std::array<mono,2> L1(const unsigned int targetParticle) const;
		std::array<mono,1> L2(const unsigned int targetParticle) const;
		std::array<mono,1> L3(const unsigned int targetParticle) const;

		// L on entire mono
		std::vector<mono> L1() const;
		std::vector<mono> L2() const;
		std::vector<mono> L3() const;


};

// calls the generic IdentifyNodes using the class's particles and (*value)
template<typename T>
inline std::vector<int> mono::IdentifyNodes(T (*value)(particle)) const{
	return IdentifyNodes([this, value](unsigned int i){return value(this->particles[i]);},
			NParticles());
}

// specialization of vector ostream template
template<>
inline std::ostream& operator<<(std::ostream& os, const std::vector<mono>& out){
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
inline std::ostream& operator<<(std::ostream& os, const std::array<mono,N>& out){
	os << "{ ";
	for(auto& element : out){
		if(element.Coeff() > 0) os << " ";
		if(element.Coeff() != 0) os << element.HumanReadable() << ", ";
	}
	return os << "\b }";
}

#endif
