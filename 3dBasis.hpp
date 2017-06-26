#ifndef BASIS_HPP
#define BASIS_HPP

#include <iostream>
#include <exception>
#include <vector>
#include <array>
#include <list>
#include <utility>		// std::pair
#include <algorithm>
#include "Eigen/Core"
#include "Eigen/Sparse"
#include "Eigen/SPQRSupport"

/******************************************************************************/
/***** Define the type used to store the coefficients of the monomials    *****/
/***** and thus also the entries of the matrix whose kernel we want.      *****/
/***** WARNING: SparseQR (at least) seems to always fail with integers!   *****/
/***** I'm not sure why this is; my code never divides coeff_class.       *****/
/******************************************************************************/

//typedef mpq_class coeff_class;				// arbitrary precision rational
//typedef mpfr::mpreal coeff_class;				// arbitrary precision float
typedef double coeff_class;						// ordinary double-width float
//typedef long coeff_class;						// ordinary double-width integer

typedef Eigen::SparseMatrix<coeff_class> Matrix;
typedef Eigen::Matrix<coeff_class, Eigen::Dynamic, Eigen::Dynamic> DMatrix;
typedef Eigen::SparseVector<coeff_class> Vector;
typedef Eigen::Matrix<coeff_class, Eigen::Dynamic, 1> DVector;
typedef Eigen::Triplet<coeff_class> Triplet;

/******************************************************************************/
/***** Define which solver to use to find the kernel of the matrix.       *****/
/***** SparseQR and SPQR are both direct (i.e. exact, non-iterative)      *****/
/***** solvers. SPQR requires an external dependency, suitesparse;        *****/
/***** SparseQR is built into Eigen, but appears to be a lot slower.      *****/
/******************************************************************************/
//typedef Eigen::SparseQR<Matrix, Eigen::COLAMDOrdering<int>> QRSolver;
typedef Eigen::SPQR<Matrix> QRSolver; // direct interface to SPQR also possible


struct particle{
	int pm;	// P_-
	int pt; // P_\perp
	int pp; // P_+
};

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
// coeff * {p1_m, p2_m, ... }{p1_t, p2_t, ...}{p1_p, p2_p, ...}
class mono {
	coeff_class coeff;
	std::vector<particle> particles;

	std::vector<int> IdentifyNodes() const;
	template<typename T> std::vector<int> IdentifyNodes(T (*value)(particle)) const;
	std::vector<int> IdentifyPmNodes() const;
	std::vector<int> IdentifyPtNodes() const;
	std::vector<int> IdentifyPpNodes() const;

	public:
		mono(): coeff(1) {}
		mono(const std::vector<int>& pm, const std::vector<int>& pt,
				const std::vector<int>& pp,	const coeff_class& coeff = 1);
		mono(const std::vector<particle>& particles, const coeff_class& coeff = 1):
			coeff(coeff), particles(particles) {}

		coeff_class& Coeff()		{ return coeff; }
		const coeff_class& Coeff() const	{ return coeff; }

		unsigned int NParticles() const { return particles.size(); }

		int Pm(const int i) const 	{ return particles[i].pm; }
		int& Pm(const int i) 		{ return particles[i].pm; }
		int Pt(const int i) const 	{ return particles[i].pt; }
		int& Pt(const int i) 		{ return particles[i].pt; }
		int Pp(const int i) const 	{ return particles[i].pp; }
		int& Pp(const int i) 	 	{ return particles[i].pp; }
		int TotalPm() const;
		int TotalPt() const;
		int TotalPp() const;
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
};

// all monos in a poly(nomial) should be guaranteed to be ordered correctly!
// Notes on arithmetic:
// * You can add or subtract monos or polys from a poly. If you do, the poly
// will check if it already has each mono; if it does, it just increases the
// coefficient, while if it doesn't it appends that mono to terms.
// * Individual component monos can be accessed with [] just like an ordinary
// array. This can also be iterated through with begin() and end().
// * The output stream operator std::cout << somePoly prints the poly as a
// sum of its constituent monos.
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

		const	mono& operator[](size_t i)	const	{ return terms[i]; }
				mono& operator[](size_t i)			{ return terms[i]; }
		size_t	size()						const	{return terms.size(); }
		std::vector<mono>::const_iterator	begin() const noexcept
				{ return terms.begin(); }
		std::vector<mono>::iterator			begin() noexcept
				{ return terms.begin(); }
		std::vector<mono>::const_iterator	end()	const noexcept
				{ return terms.end(); }
		std::vector<mono>::iterator			end()	noexcept
				{ return terms.end(); }

		friend std::ostream& operator<<(std::ostream& os, const poly& out);

		poly DerivPm(const int);
		poly DerivPp(const int);
		poly DerivPm();
		poly DerivPp();

		poly K1() const;
		poly K2() const;
		poly K3() const;

		static poly K1(const mono& inputMono);
		static poly K2(const mono& inputMono);
		static poly K3(const mono& inputMono);
};

// a set of all ordered monomials of the given degree and number of particles
// * Can be accessed with [] like an array and iterated through with begin()
// and end().
// * Has a stream operator, which outputs all the component monos on new lines.
// * Most important function is ExpressPoly, which takes a polynomial and
// expresses it as a vector in this basis's space.
class basis {
	std::vector<mono> basisMonos;

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
		unsigned int FindInBasis(const std::vector<int>& pm,
				const std::vector<int>& pt, const std::vector<int>& pp);
		unsigned int FindInBasis(const mono& wildMono);

		//std::array<basis,2> ParitySplit() const;
		void DeleteOdd();
		void DeleteEven();

		const	mono& operator[](size_t i)	const	{return basisMonos[i];     }
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
};

// Contains two bases and intelligently decides which one to use for various
// requests. Because this is in two pieces, it can not be iterated through.
class splitBasis {
	basis evenBasis;
	basis oddBasis;

	public:
		splitBasis() = delete;
		splitBasis(const int numP, const int degree);

		std::pair<unsigned int, basis*> FindInBasis(const std::vector<int>& pm,
				const std::vector<int>& pt, const std::vector<int>& pp);
		std::pair<unsigned int, basis*> FindInBasis(const mono& wildMono);

		basis& OddBasis() { return oddBasis; }
		basis& EvenBasis() { return evenBasis; }

		std::list<Triplet> ExpressPoly(const poly& toExpress, const int column,
				const int row) const;

		static bool IsOdd(const mono& toTest);
		static bool IsEven(const mono& toTest);
};

std::ostream& operator<<(std::ostream& os, const Triplet& out);
std::vector<std::vector<int>> Permute(const std::vector<std::vector<int>> ordered);
std::vector<std::vector<int>> GetStatesByDegree(const int numP, const int deg,
		const bool exact, const int min);
std::vector<std::vector<int>> GetStatesUpToDegree(const int numP, const int deg);
std::vector<std::vector<int>> GetStatesAtDegree(const int numP, const int deg);


// functions interfacing with Eigen ------------------------------------------

Matrix KMatrix(const basis& startingBasis, const basis& targetBasis);
std::list<Triplet> ConvertToRows(const std::vector<poly>& polyForms, 
		const basis& targetBasis, const Eigen::Index rowOffset);
std::vector<poly> Kernel(const Matrix& KActions, const basis& startBasis);
poly VectorToPoly(const Vector& kernelVector, const basis& startBasis);
poly VectorToPoly(const DVector& kernelVector, const basis& startBasis);
poly ColumnToPoly(const Matrix& kernelMatrix, const Eigen::Index col, 
		const basis& startBasis);
poly ColumnToPoly(const DMatrix& kernelMatrix, const Eigen::Index col, 
		const basis& startBasis);


// templates -----------------------------------------------------------------

// T can be any type or class with an == operator; value indexes T by uint
template<typename Accessor>
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

// default behavior for a vector of particles; used by mono::IdentifyNodes()
template<>
inline std::vector<int> IdentifyNodes(const std::vector<particle>& particles){
	return IdentifyNodes([particles](unsigned int i){return 
			std::array<int, 3>({{particles[i].pm, particles[i].pt, particles[i].pp}}); },
			particles.size());
}

// stream output operator template for vectors
template<typename T>
inline std::ostream& operator<<(std::ostream& os, const std::vector<T>& out){
	os << "{";
	for(auto& element : out){
		if(element >= 0) os << " ";
		os << element << ",";
	}
	os << "\b }";
	return os;
}

// specialization of above template which "transposes" particle vectors
template<>
inline std::ostream& operator<<(std::ostream& os, const std::vector<particle>& out){
	os << "{";
	for(auto& p : out){
		if(p.pm >= 0) os << " ";
		os << p.pm << ",";
	}
	os << "\b }{";
	for(auto& p : out){
		if(p.pt >= 0) os << " ";
		os << p.pt << ",";
	}
	os << "\b }{";
	for(auto& p : out){
		if(p.pp >= 0) os << " ";
		os << p.pp << ",";
	}
	os << "\b }";
	return os;
}

#endif
