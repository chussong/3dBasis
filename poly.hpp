#ifndef POLY_HPP
#define POLY_HPP

#include <vector>
#include <string>
#include "mono.hpp"
#include "io.hpp"

// all monos in a poly(nomial) should be guaranteed to be ordered correctly!
// Notes on arithmetic:
// * You can add or subtract monos or Polys from a Poly. If you do, the Poly
// will check if it already has each mono; if it does, it just increases the
// coefficient, while if it doesn't it appends that mono to terms.
// * You can multiply and divide by basic arithmetic types but not mono or Poly.
// * Individual component monos can be accessed with [] just like an ordinary
// array. This can also be iterated through with begin() and end().
// * The output stream operator std::cout << somePoly prints the Poly as a
// sum of its constituent monos.
class Poly {
    std::vector<Mono> terms;

    public:
        explicit Poly() {}
        explicit Poly(const Mono starter) { terms.push_back(starter); }
        explicit Poly(const std::vector<Mono>& terms);

        Poly& operator+=(const Mono& x);
        Poly& operator-=(const Mono& x);
        Poly& operator+=(const Poly& x);
        Poly& operator-=(const Poly& x);

        friend Poly operator+(Poly x, const Poly& y);
        friend Poly operator-(Poly x, const Poly& y);
        friend Poly operator+(Poly x, const Mono& y);
        friend Poly operator-(Poly x, const Mono& y);
        friend Poly operator+(const Mono& x, Poly y);
        friend Poly operator-(const Mono& x, Poly y);

        bool operator==(const Poly& other) const;
        bool operator!=(const Poly& other) const { return !(*this == other); }

        template<typename T> Poly& operator*=(const T& y);
        template<typename T> Poly& operator/=(const T& y);
        template<typename T>
                friend Poly operator*(Poly x, const T&    y) { return x*=y; }
        template<typename T>
                friend Poly operator/(Poly x, const T&    y) { return x/=y; }
        template<typename T>
                friend Poly operator*(const T& x,    Poly y) { return y*x; }
        Poly operator-() const;

        const	Mono& operator[](size_t i)      const   { return terms[i]; }
        Mono& operator[](size_t i)                      { return terms[i]; }
        size_t	size()                          const   {return terms.size(); }
        std::vector<Mono>::const_iterator               begin() const noexcept
                        { return terms.begin(); }
        std::vector<Mono>::iterator                     begin() noexcept
                        { return terms.begin(); }
        std::vector<Mono>::const_iterator               end()   const noexcept
                        { return terms.end(); }
        std::vector<Mono>::iterator                     end()   noexcept
                        { return terms.end(); }

        friend std::ostream& operator<<(std::ostream& os, const Poly& out);
        std::string HumanReadable() const;

        Poly DerivPm(const int);
        Poly DerivPm();

};

template<typename T>
Poly& Poly::operator*=(const T& y){
    for(Mono& m : terms){
        m *= y;
    }
    return *this;
}

template<typename T>
Poly& Poly::operator/=(const T& y){
    for(Mono& m : terms){
        m /= y;
    }
    return *this;
}

/*template<typename T>
Poly operator*(Poly x, const T& y){
    Poly ret(x);
    for(Mono& m : ret){
        m *= y;
    }
    return ret;
}*/

/*template<typename T>
Poly operator/(Poly x, const T& y){
    Poly ret(x);
    for(Mono& m : ret){
        m /= y;
    }
    return ret;
}*/

// specialization of vector output template for vectors of polynomials
template<>
inline std::ostream& operator<<(std::ostream& os, const std::vector<Poly>& out){
    if(out.size() == 0) return os << "{ }";
    os << "{ ";
    for(auto& element : out){
        os << element.HumanReadable() << " | ";
    }
    os << "\b\b \b}";
    return os;
}

#endif
