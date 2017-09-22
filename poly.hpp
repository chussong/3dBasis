#ifndef POLY_HPP
#define POLY_HPP

#include <vector>
#include <string>
#include "mono.hpp"
#include "io.hpp"
#include "cache.hpp"

// all monos in a poly(nomial) should be guaranteed to be ordered correctly!
// Notes on arithmetic:
// * You can add or subtract monos or polys from a poly. If you do, the poly
// will check if it already has each mono; if it does, it just increases the
// coefficient, while if it doesn't it appends that mono to terms.
// * You can multiply and divide by basic arithmetic types but not mono or poly.
// * Individual component monos can be accessed with [] just like an ordinary
// array. This can also be iterated through with begin() and end().
// * The output stream operator std::cout << somePoly prints the poly as a
// sum of its constituent monos.
class poly {
	std::vector<mono> terms;

	public:
		explicit poly() {}
		explicit poly(const mono starter) { terms.push_back(starter); }
		explicit poly(const std::vector<mono>& terms);

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

		bool operator==(const poly& other) const;
		bool operator!=(const poly& other) const { return !(*this == other); }

		template<typename T> poly& operator*=(const T& y);
		template<typename T> poly& operator/=(const T& y);
		template<typename T>
			friend poly operator*(poly x, const T&    y) { return x*=y; }
		template<typename T>
			friend poly operator/(poly x, const T&    y) { return x/=y; }
		template<typename T>
			friend poly operator*(const T& x,    poly y) { return y*x; }
		poly operator-() const;

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
		std::string HumanReadable() const;

		poly DerivPm(const int);
		poly DerivPp(const int);
		poly DerivPm();
		poly DerivPp();

		void MirrorPM() { for(auto& t : terms) t.MirrorPM(); }

		poly K1(const coeff_class delta) const;
		poly K2(const coeff_class delta) const;
		poly K3(const coeff_class delta) const;

		static poly K1(const mono& inputMono, const coeff_class delta);
		static poly K2(const mono& inputMono, const coeff_class delta);
		static poly K3(const mono& inputMono, const coeff_class delta);

		static poly DeleteNonDirichlet(const poly& inputPoly);
		static poly DeleteNonDirichlet(const std::vector<mono>& inputMonos);

		static coeff_class InnerProduct(const poly& A, const poly& B,
						const GammaCache& cache, const KVectorCache& kCache);
};

template<typename T>
poly& poly::operator*=(const T& y){
	for(mono& m : terms){
		m *= y;
	}
	return *this;
}

template<typename T>
poly& poly::operator/=(const T& y){
	for(mono& m : terms){
		m /= y;
	}
	return *this;
}

/*template<typename T>
poly operator*(poly x, const T& y){
	poly ret(x);
	for(mono& m : ret){
		m *= y;
	}
	return ret;
}*/

/*template<typename T>
poly operator/(poly x, const T& y){
	poly ret(x);
	for(mono& m : ret){
		m /= y;
	}
	return ret;
}*/

// specialization of vector output template for vectors of polynomials
template<>
inline std::ostream& operator<<(std::ostream& os, const std::vector<poly>& out){
	if(out.size() == 0) return os << "{ }";
	os << "{ ";
	for(auto& element : out){
		os << element.HumanReadable() << " | ";
	}
	os << "\b\b \b}";
	return os;
}

#endif
