#ifndef BASIS_HPP
#define BASIS_HPP

#include <exception>
#include <algorithm>
#include <vector>
#include <array>
#include <list>
#include <utility>		// std::pair
#include <iostream>
#include <gmpxx.h>
#include "mpreal.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <unsupported/Eigen/MPRealSupport>

// Eigen extension for mpq_class plus a stream function
/*
inline std::ostream& operator<<(std::ostream &o, const coeff_class &expr){
	return o << expr.get_str();
}

namespace Eigen {
	template<> struct NumTraits<mpq_class> : GenericNumTraits<mpq_class>{
		typedef mpq_class Real;
		typedef mpq_class NonInteger;
		typedef mpq_class Literal;
		typedef mpq_class Nested;

		enum {
			IsComplex = 0,
			IsInteger = 0,
			ReadCost = 6,
			AddCost = 150,
			MulCost = 100,
			IsSigned = 1,
			RequireInitialization = 1
		};

		static inline Real epsilon() { return 0; }
		static inline Real dummy_precision() { return 0; }
		static inline Real digits10() { return 0; }
		static inline Real Pi() { return mpq_class(355,113); }
	};

	 namespace internal {

	template<> inline mpq_class random<mpq_class>()
	{
	  gmp_randclass rr(gmp_randinit_default);
	  return mpq_class(rr.get_z_bits(125),rr.get_z_bits(125));
	}

	template<> inline mpq_class random<mpq_class>(const mpq_class& a, const mpq_class& b)
	{
	  return a + (b-a) * random<mpq_class>();
	}

	inline bool isMuchSmallerThan(const mpq_class& a, const mpq_class& b, const mpq_class& eps)
	{
	  return ::abs(a) <= ::abs(b) * eps;
	}

	inline bool isApprox(const mpq_class& a, const mpq_class& b, const mpq_class& eps)
	{
		return ::abs(a-b)<eps;
	}

	inline bool isApproxOrLessThan(const mpq_class& a, const mpq_class& b, const mpq_class& eps)
	{
	  return a <= b;
	}

	template<> inline long double cast<mpq_class,long double>(const mpq_class& x)
	{ return x.get_d(); }

	template<> inline double cast<mpq_class,double>(const mpq_class& x)
	{ return x.get_d(); }

	template<> inline long cast<mpq_class,long>(const mpq_class& x)
	{ return x.get_d(); }

	template<> inline int cast<mpq_class,int>(const mpq_class& x)
	{ return x.get_d(); }

	// G+Smo
	template<> inline size_t cast<mpq_class,size_t>(const mpq_class& x)
	{ return x.get_d(); }

	template<> inline unsigned cast<mpq_class,unsigned>(const mpq_class& x)
	{ return x.get_d(); }

	  // Specialize GEBP kernel and traits
	  template<>
	  class gebp_traits<mpq_class, mpq_class, false, false>
	  {
	  public:
		typedef mpq_class ResScalar;
		enum {
		  Vectorizable = false,
		  LhsPacketSize = 1,
		  RhsPacketSize = 1,
		  ResPacketSize = 1,
		  NumberOfRegisters = 1,
		  nr = 1,
		  mr = 1,
		  LhsProgress = 1,
		  RhsProgress = 1
		};
		typedef ResScalar LhsPacket;
		typedef ResScalar RhsPacket;
		typedef ResScalar ResPacket;

	  };

	  template<typename Index, typename DataMapper, bool ConjugateLhs, bool ConjugateRhs>
	  struct gebp_kernel<mpq_class,mpq_class,Index,DataMapper,1,1,ConjugateLhs,ConjugateRhs>
	  {
		typedef mpq_class num_t;

		EIGEN_DONT_INLINE
		void operator()(const DataMapper& res, const num_t* blockA, const num_t* blockB, 
						Index rows, Index depth, Index cols, const num_t& alpha,
						Index strideA=-1, Index strideB=-1, Index offsetA=0, Index offsetB=0)
		{
		  if(rows==0 || cols==0 || depth==0)
			return;

		  num_t  acc1(0), tmp(0);        

		  if(strideA==-1) strideA = depth;
		  if(strideB==-1) strideB = depth;

		  for(Index i=0; i<rows; ++i)
		  {
			for(Index j=0; j<cols; ++j)
			{
			  const num_t *A = blockA + i*strideA + offsetA;
			  const num_t *B = blockB + j*strideB + offsetB;
			  acc1 = 0;
			  for(Index k=0; k<depth; k++)
			  {
				mpq_mul(tmp.__get_mp() , A[k].__get_mp(), B[0].__get_mp());
				mpq_add(acc1.__get_mp(), acc1.__get_mp(), tmp.__get_mp() );
			  }
			  
			  mpq_mul(acc1.__get_mp()    , acc1.__get_mp() , alpha.__get_mp());
			  mpq_add(res(i,j).__get_mp(), res(i,j).__get_mp(), acc1.__get_mp() );
			}
		  }
		}
	  };

	} // end namespace internal
}*/

/***** define the type used to store the coefficients of the monomials    *****/
/***** and thus also the entries of the matrix whose kernel we want.      *****/
//typedef mpq_class coeff_class;				// arbitrary precision rational
//typedef mpfr::mpreal coeff_class;				// arbitrary precision float
typedef double coeff_class;					// ordinary double-width float
//typedef int coeff_class;						// ordinary single-width integer

typedef Eigen::SparseMatrix<coeff_class> Matrix;
typedef Eigen::Triplet<coeff_class> Triplet;
typedef Eigen::SparseQR<Matrix, Eigen::COLAMDOrdering<int>> QRSolver;
// SPQR is a threaded alternative to SparseQR that requires external linking

struct particle{
	int pm;	// P_-
	int pt; // P_\perp
	int pp; // P_+
};

// a mono(mial) with coefficient. It should be impossible for an instance of
// this class to be out of order, so hopefully that's true!
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

		std::vector<mono>& BasisMonos() { return basisMonos; }
		const std::vector<mono>& BasisMonos() const { return basisMonos; }

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

		std::list<Triplet> ExpressPoly(const poly& toExpress, const int column, 
				const int row) const;
};

class splitBasis {
	basis oddBasis;
	basis evenBasis;

	public:
		splitBasis() = delete;
		splitBasis(const int numP, const int degree);

		std::pair<unsigned int, basis*> FindInBasis(const std::vector<int>& pm,
				const std::vector<int>& pt, const std::vector<int>& pp);
		std::pair<unsigned int, basis*> FindInBasis(const mono& wildMono);

		basis& Odd() { return oddBasis; }
		basis& Even() { return evenBasis; }

		std::list<Triplet> ExpressPoly(const poly& toExpress, const int column,
				const int row) const;

		static bool IsOdd(const mono& toTest);
};

/*std::ostream& operator<<(std::ostream& os, const std::vector<int>& out);
std::ostream& operator<<(std::ostream& os, const std::vector<coeff_class>& out);
std::ostream& operator<<(std::ostream& os, const std::vector<particle>& out);*/
std::ostream& operator<<(std::ostream& os, const Triplet& out);
//std::ostream& operator<<(std::ostream& os, const Matrix& out);
std::vector<std::vector<int>> Permute(const std::vector<std::vector<int>> ordered);
std::vector<std::vector<int>> GetStatesByDegree(const int numP, const int deg,
		const bool exact, const int min);
std::vector<std::vector<int>> GetStatesUpToDegree(const int numP, const int deg);
std::vector<std::vector<int>> GetStatesAtDegree(const int numP, const int deg);

// functions interfacing with Eigen
Matrix KMatrix(const basis& startingBasis, const basis& targetBasis);
std::list<Triplet> ConvertToRows(const std::vector<poly>& polyForms, 
		const basis& targetBasis, const Eigen::Index rowOffset);
std::vector<poly> Kernel(const Matrix& KActions, const basis& startBasis);
poly ColumnToPoly(const Matrix& Q, const Eigen::Index col, const basis& startBasis);

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
