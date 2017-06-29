#ifndef POLY_HPP
#define POLY_HPP

#include <vector>
#include <string>
#include "mono.hpp"

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

		bool operator==(const poly& other) const;
		bool operator!=(const poly& other) const { return !(*this == other); }

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
		std::string HumanReadable();

		poly DerivPm(const int);
		poly DerivPp(const int);
		poly DerivPm();
		poly DerivPp();

		void MirrorPM() { for(auto& t : terms) t.MirrorPM(); }

		poly K1() const;
		poly K2() const;
		poly K3() const;

		static poly K1(const mono& inputMono);
		static poly K2(const mono& inputMono);
		static poly K3(const mono& inputMono);
};

#endif
