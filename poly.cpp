#include "poly.hpp"

// the reason this uses += rather than a copy is so that coefficients can be
// added for monomials which appear more than once
Poly::Poly(const std::vector<Mono>& terms){
    if(terms.size() == 0){
        std::cerr << "Warning! A Poly has been constructed from an empty Mono."
            << std::endl;
    }
    for(auto& newTerm : terms){
        *this += newTerm;
    }
}

Poly& Poly::operator+=(const Mono& x){
    if(std::abs<builtin_class>(x.Coeff()) < EPSILON) return *this;

    for(auto it = terms.begin(); it != terms.end(); ++it) {
        if(*it == x){
            it->Coeff() += x.Coeff();
            if(std::abs<builtin_class>(it->Coeff()) < EPSILON) terms.erase(it);
            return *this;
        }
    }
    terms.push_back(x);
    return *this;
}

Poly& Poly::operator-=(const Mono& x){
    return *this += -x;
}

Poly& Poly::operator+=(const Poly& x){
    for(auto& mn : x.terms){
        *this += mn;
    }
    return *this;
}

Poly& Poly::operator-=(const Poly& x){
    for(auto& mn : x.terms){
        *this -= mn;
    }
    return *this;
}

Poly operator+(Poly x, const Poly& y){
    return x += y;
}

Poly operator-(Poly x, const Poly& y){
    return x -= y;
}

Poly operator+(Poly x, const Mono& y){
    return x += y;
}

Poly operator-(Poly x, const Mono& y){
    return x -= y;
}

Poly operator+(const Mono& x, Poly y){
    return y += x;
}

Poly operator-(const Mono& x, Poly y){
    return y -= x;
}

Poly Poly::operator-() const{
    Poly ret(*this);
    for(auto& m : ret) m = -m;
    return ret;
}

bool Poly::operator==(const Poly& other) const{
    bool found;
    for(auto& term1 : terms){
        found = false;
        for(auto& term2 : other.terms){
            if(term1 == term2){
                if(std::abs<builtin_class>(term1.Coeff() - term2.Coeff()) > EPSILON) {
                    return false;
                }
                found = true;
                break;
            }
        }
        if(!found) return false;
    }
    return true;
}

/*OStream& operator<<(OStream& os, const Poly& out){
    for(auto& component : out) os << component << " + ";
    return out.size() > 0 ? os << "\b\b \b\b" : os;
}*/

std::ostream& operator<<(std::ostream& os, const Poly& out){
    return out.size() > 0 ? os << out.HumanReadable() : os;
}

std::string Poly::HumanReadable() const {
    if(terms.size() == 0) return "";
    std::ostringstream os;
    if (terms[0].Coeff() < 0) os << "- ";
    os << terms[0].HumanReadable();
    for (std::size_t i = 1; i < terms.size(); ++i) {
        if (terms[i].Coeff() < 0) {
            os << " - ";
        } else {
            os << " + ";
        }
        os << terms[i].HumanReadable();
    }
    return os.str();
}
