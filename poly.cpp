#include "poly.hpp"

// the reason this uses += rather than a copy is so that coefficients can be
// added for monomials which appear more than once
poly::poly(const std::vector<mono>& terms){
	if(terms.size() == 0){
		std::cout << "Warning! A poly has been constructed from an empty mono."
			<< std::endl;
	}
	for(auto& newTerm : terms){
		*this += newTerm;
	}
}

poly& poly::operator+=(const mono& x){
	if(x.Coeff() == 0) return *this;

	for(auto& tm : terms){
		if(tm == x){
			tm.Coeff() += x.Coeff();
			return *this;
		}
	}
	terms.push_back(x);
	return *this;
}

poly& poly::operator-=(const mono& x){
	return *this += -x;
}

poly& poly::operator+=(const poly& x){
	for(auto& mn : x.terms){
		*this += mn;
	}
	return *this;
}

poly& poly::operator-=(const poly& x){
	for(auto& mn : x.terms){
		*this -= mn;
	}
	return *this;
}

poly operator+(poly x, const poly& y){
	return x += y;
}

poly operator-(poly x, const poly& y){
	return x -= y;
}

poly operator+(poly x, const mono& y){
	return x += y;
}

poly operator-(poly x, const mono& y){
	return x -= y;
}

poly operator+(const mono& x, poly y){
	return y += x;
}

poly operator-(const mono& x, poly y){
	return y -= x;
}

bool poly::operator==(const poly& other) const{
	bool found;
	for(auto& term1 : terms){
		found = false;
		for(auto& term2 : other.terms){
			if(term1 == term2){
				if(term1.Coeff() != term2.Coeff()) return false;
				found = true;
				break;
			}
		}
		if(!found) return false;
	}
	return true;
}

std::ostream& operator<<(std::ostream& os, const poly& out){
	for(auto& component : out) os << component << " + ";
	return out.size() > 0 ? os << "\b\b \b\b" : os;
}

/*std::string poly::HumanReadable(){
	if(terms.size() == 0) return "";
	std::string ret = "";
	if(terms[0].Coeff() < 0) ret.append("  ");
	for(auto& term : terms){
		if(term.Coeff() < 0) ret.append("\b\b- ");
		if(std::abs(term.Coeff()) != 1){
			ret.append(std::to_string(std::abs(term.Coeff())));
		}
		ret.append(term.HumanReadable() + " + ");
	}
	ret.erase(ret.end()-3, ret.end());
	return ret;
}*/

std::string poly::HumanReadable() const{
	if(terms.size() == 0) return "";
	std::ostringstream os;
	if(terms[0].Coeff() < 0) os << "  ";
	for(auto& term : terms){
		if(term.Coeff() < 0) os << "\b\b- ";
		if(std::abs(term.Coeff()) != 1) os << std::abs(term.Coeff());
		os << term.HumanReadable() << " + ";
	}
	os << "\b\b\b"; // this will leave a '+' around if <2 more chars are written
	return os.str();
}

poly poly::K1(const coeff_class delta) const{
	poly ret;
	for(auto& term : terms){
		for(auto& newTerm : term.K1(delta)) ret += newTerm;
	}
	return ret;
}

poly poly::K2(const coeff_class delta) const{
	poly ret;
	for(auto& term : terms){
		for(auto& newTerm : term.K2(delta)) ret += newTerm;
	}
	return ret;
}

poly poly::K3(const coeff_class delta) const{
	poly ret;
	for(auto& term : terms){
		for(auto& newTerm : term.K3(delta)) ret += newTerm;
	}
	return ret;
}

