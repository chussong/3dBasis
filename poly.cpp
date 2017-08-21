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
	if(std::abs(x.Coeff()) < EPSILON) return *this;

	for(auto it = terms.begin(); it != terms.end(); ++it){
		if(*it == x){
			it->Coeff() += x.Coeff();
			if(std::abs(it->Coeff()) < EPSILON) terms.erase(it);
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

poly poly::operator-() const{
	poly ret(*this);
	for(auto& m : ret) m = -m;
	return ret;
}

bool poly::operator==(const poly& other) const{
	bool found;
	for(auto& term1 : terms){
		found = false;
		for(auto& term2 : other.terms){
			if(term1 == term2){
				if(std::abs(term1.Coeff() - term2.Coeff()) > EPSILON) return false;
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
		//if(std::abs(term.Coeff()) != 1) os << std::abs(term.Coeff());
		os << term.HumanReadable() << " + ";
	}
	os << "\b\b \b\b"; // this will leave a '+' around if <2 more chars are written
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

poly poly::DeleteNonDirichlet(const poly& inputPoly){
	poly ret;
	for(auto& term : inputPoly){
		if(term.IsDirichlet()) ret += term;
	}
	return ret;
}

poly poly::DeleteNonDirichlet(const std::vector<mono>& inputMonos){
	poly ret;
	for(auto& term : inputMonos){
		if(term.IsDirichlet()) ret += term;
	}
	return ret;
}

coeff_class poly::InnerProduct(const poly& A, const poly& B){
	coeff_class ret = 0;
	for(auto& monoA : A){
		for(auto& monoB : B){
			ret += mono::InnerProduct(monoA, monoB);
		}
	}
	return ret;
}
