#include "mono.hpp"

mono::mono(const std::vector<particle>& particles, const bool usingEoM, 
		const coeff_class& coeff): coeff(coeff), particles(particles), 
		usingEoM(usingEoM){
	if(usingEoM) Order();
}

mono::mono(const std::vector<int>& pm, const std::vector<int>& pt,
		const std::vector<int>& pp,	const bool usingEoM,
		const coeff_class& coeff): coeff(coeff), usingEoM(usingEoM){
	if(pm.size() != pt.size() || pt.size() != pp.size()){
		std::cerr << "Error: attempted to construct a monomial out of particle "
			<< "data with different sizes: {" << pm.size() << "," << pt.size()
			<< "," << pp.size() << "}. It will be blank instead." << std::endl;
		return;
	}
	particles.resize(pm.size());
	for(auto i = 0u; i < pm.size(); ++i){
		particles[i].pm = pm[i];
		particles[i].pt = pt[i];
		particles[i].pp = pp[i];
	}
	if(usingEoM) Order();
}

// note: like all operations on completed monos, this assumes both are ordered
bool mono::operator==(const mono& other) const{
	if(particles.size() != other.particles.size()){
		std::cerr << "Error: asked to compare two monomials with different "
			<< "numbers of particles." << std::endl;
		return false;
	}
	for(unsigned int i = 0; i < particles.size(); ++i){
		if(particles[i].pm != other.particles[i].pm
				|| particles[i].pt != other.particles[i].pt
				|| particles[i].pp != other.particles[i].pp) return false;
	}
	return true;
}

std::ostream& operator<<(std::ostream& os, const mono& out){
	if(out.coeff == 1) return os << out.particles;
	return os << out.coeff << " * " << out.particles;
}

/*std::string mono::HumanReadable(){
	std::string ret = "";
	for(auto& p : particles){
		if(p.pm != 0){
			ret.append("M");
			if(p.pm != 1) ret.append("^" + std::to_string(p.pm));
		}
		if(p.pt != 0){
			ret.append("T");
			if(p.pt != 1) ret.append("^" + std::to_string(p.pt));
		}
		if(p.pp != 0){
			ret.append("P");
			if(p.pp != 1) ret.append("^" + std::to_string(p.pp));
		}
		ret.append("Φ");
	}
	return ret;
}*/

std::string mono::HumanReadable() const{
	std::ostringstream os;
	for(auto& p : particles){
		if(p.pm != 0){
			os << "M";
			if(p.pm != 1) os << "^" << p.pm;
		}
		if(p.pt != 0){
			os << "T";
			if(p.pt != 1) os << "^" << p.pt;
		}
		if(p.pp != 0){
			os << "P";
			if(p.pp != 1) os << "^" << p.pp;
		}
		os << "Φ";
	}
	return os.str();
}

mono mono::operator-() const{
	mono ret(*this);
	ret.coeff *= -1;
	return ret;
}

int mono::TotalPm() const{
	int total = 0;
	for(auto& p : particles) total += p.pm;
	return total;
}

int mono::TotalPt() const{
	int total = 0;
	for(auto& p : particles) total += p.pt;
	return total;
}
int mono::TotalPp() const{
	int total = 0;
	for(auto& p : particles) total += p.pp;
	return total;
}

// if the mono is ordered, particles[0] is guaranteed to have the highest Pm
int mono::MaxPm() const{
	if(particles.size() < 1) return -1;
	return particles[0].pm;
}

// for this and MaxPp(), we have to actually check
int mono::MaxPt() const{
	int max = -1;
	for(auto& p : particles) if(p.pt > max) max = p.pt;
	return max;
}

// technically we could interrupt this when p.pp < max
int mono::MaxPp() const{
	int max = -1;
	for(auto& p : particles) if(p.pp > max) max = p.pp;
	return max;
}

// The following two functions might actually work just as well as old Order():
bool mono::ParticlePrecedence(const particle& a, const particle& b){
	if(a.pm != b.pm) return a.pm > b.pm;
	if(a.pt != b.pt) return a.pt > b.pt;
	return a.pp > b.pp;
}

void mono::Order(){
	// if EoM in use, 2MP == T^2
	// this variant only eliminates +- pairs on the same particle
	/*if(usingEoM){
		for(auto& p : particles){
			while(p.pm > 0 && p.pp > 0){
				coeff /= 2; // I think this is right...
				p.pm -= 1;
				p.pp -= 1;
				p.pt += 2;
			}
		}
	}*/
	// this variant eliminates + completely, turning it into Pt^2/(2Pm)
	if(usingEoM){
		for(auto& p : particles){
			if(p.pp > 0){
				coeff /= std::pow(2,p.pp); // no bit shift in case coeff < 0
				p.pm -= p.pp;
				p.pt += 2*p.pp;
				p.pp = 0;
			}
		}
	}
	std::sort(particles.begin(), particles.end(), ParticlePrecedence);
}

mono mono::OrderCopy() const{
	mono ret(*this);
	ret.Order();
	return ret;
}

bool mono::MirrorIsBetter(const mono& original){
	if(original.TotalPm() != original.TotalPp())
		return original.TotalPp() > original.TotalPm();
	if(original.MaxPm() != original.MaxPp())
		return original.MaxPp() > original.MaxPm();
	mono mirror(original);
	mirror.MirrorPM();
	std::cout << "Compare " << original << " to its mirror: " << mirror << std::endl;

	for(auto i = 0u; i < original.NParticles(); ++i){
		if(original.Pt(i) != mirror.Pt(i)) return mirror.Pt(i) > original.Pt(i);
	}

	std::cout << "Asked to compare " << original << " to its mirror but they "
		<< "were the same..." << std::endl;

	return false;
}

/*
// this should be called whenever something happens that changes the momenta in
// a monomial (new or existing) to ensure that it's always ordered. To ensure
// the original operation isn't disrupted, call it in the scope of the function
// which is iterating over the individual particles, not within a particle-scope
// 1: order by pm
// 2: identify pm nodes, order by pt
// 3: identify pm&pt nodes, order by pp
void mono::Order(){
	std::sort(particles.begin(), particles.end(),
			[](particle a, particle b){return a.pm > b.pm;});

	std::vector<int> nodes(IdentifyPmNodes());
	auto i = 0u;
	for(auto& node : nodes){
		std::sort(particles.begin() + i, particles.begin() + i + node,
				[](particle a, particle b){return a.pt > b.pt;});
		i += node;
	}

	nodes = ::IdentifyNodes([this](unsigned int i){
			return std::array<int,2>({{Pm(i), Pt(i)}}); }, NParticles() );
	i = 0u;
	for(auto& node : nodes){
		std::sort(particles.begin() + i, particles.begin() + i + node,
				[](particle a, particle b){return a.pp > b.pp;});
		i += node;
	}
}*/

std::vector<int> mono::IdentifyNodes() const{
	return ::IdentifyNodes(particles);
}

std::vector<int> mono::IdentifyPmNodes() const{
	return ::IdentifyNodes([this](unsigned int i){return Pm(i);}, NParticles());
}

std::vector<int> mono::IdentifyPtNodes() const{
	return ::IdentifyNodes([this](unsigned int i){return Pt(i);}, NParticles());
}

std::vector<int> mono::IdentifyPpNodes() const{
	return ::IdentifyNodes([this](unsigned int i){return Pp(i);}, NParticles());
}

void mono::MirrorPM(){
	int temp;
	for(auto& p : particles){
		temp = p.pm;
		p.pm = p.pp;
		p.pp = temp;
	}
	Order();
}

// NOTE! These break ordering, so you have to re-order when you're done!
mono mono::DerivPm(const unsigned int part) const{
	if(part >= particles.size()){
		std::cerr << "Error: monomial told to take a derivative of momentum Pm["
			<< part << "], but it only knows about " << NParticles() << "."
			<< std::endl;
		return *this;
	}
	mono ret(*this);
	ret.coeff *= ret.Pm(part);
	ret.Pm(part)--;
	return ret;
}

mono mono::DerivPt(const unsigned int part) const{
	if(part >= particles.size()){
		std::cerr << "Error: monomial told to take a derivative of momentum Pt["
			<< part << "], but it only knows about " << NParticles() << "."
			<< std::endl;
		return *this;
	}
	mono ret(*this);
	ret.coeff *= ret.Pt(part);
	ret.Pt(part)--;
	return ret;
}

mono mono::DerivPp(const unsigned int part) const{
	if(part >= particles.size()){
		std::cerr << "Error: monomial told to take a derivative of momentum Pp["
			<< part << "], but it only knows about " << NParticles() << "."
			<< std::endl;
		return *this;
	}
	mono ret(*this);
	ret.coeff *= ret.Pp(part);
	ret.Pp(part)--;
	return ret;
}

std::vector<mono> mono::DerivPm() const{
	std::vector<mono> ret;
	for(unsigned int particle = 0; particle < NParticles(); ++particle){
		ret.push_back(DerivPm(particle).OrderCopy());
		ret[particle].Order();
	}
	return ret;
}

std::vector<mono> mono::DerivPt() const{
	std::vector<mono> ret;
	for(unsigned int particle = 0; particle < NParticles(); ++particle){
		ret.push_back(DerivPt(particle).OrderCopy());
		ret[particle].Order();
	}
	return ret;
}

std::vector<mono> mono::DerivPp() const{
	std::vector<mono> ret;
	for(unsigned int particle = 0; particle < NParticles(); ++particle){
		ret.push_back(DerivPp(particle).OrderCopy());
	}
	return ret;
}

mono mono::MultPm(const unsigned int particle) const{
	mono ret(*this);
	++ret.Pm(particle);
	return ret;
}

mono mono::MultPt(const unsigned int particle) const{
	mono ret(*this);
	++ret.Pt(particle);
	return ret;
}

mono mono::MultPp(const unsigned int particle) const{
	mono ret(*this);
	if(usingEoM){
		ret.Pm(particle) -= 1;
		ret.Pt(particle) += 2;
	} else {
		ret.Pp(particle) += 1;
	}
	return ret;
}

std::array<mono, 4> mono::K1(const unsigned int particle,
		const coeff_class delta) const{
	std::array<mono, 4> ret({{
			2*this->DerivPp(particle).DerivPp(particle).MultPp(particle),
			2*this->DerivPt(particle).DerivPp(particle).MultPt(particle),
			2*delta*this->DerivPp(particle),
			this->DerivPt(particle).DerivPt(particle).MultPm(particle)}});
	return ret;
}

std::array<mono, 5> mono::K2(const unsigned int particle,
		const coeff_class delta) const{
	std::array<mono, 5> ret({{
			-2*this->DerivPt(particle).DerivPp(particle).MultPp(particle),
			-2*this->DerivPt(particle).DerivPm(particle).MultPm(particle),
			-this->DerivPt(particle).DerivPt(particle).MultPt(particle),
			-2*delta*this->DerivPt(particle),
			-2*this->DerivPm(particle).DerivPp(particle).MultPt(particle)}});
	return ret;
}

std::array<mono, 4> mono::K3(const unsigned int particle,
		const coeff_class delta) const{
	std::array<mono, 4> ret({{
			2*this->DerivPm(particle).DerivPm(particle).MultPm(particle),
			2*this->DerivPt(particle).DerivPm(particle).MultPt(particle),
			2*delta*this->DerivPm(particle),
			this->DerivPt(particle).DerivPt(particle).MultPp(particle)}});
	return ret;
}

std::vector<mono> mono::K1(const coeff_class delta) const{
	std::vector<mono> ret;
	for(unsigned int i = 0; i < NParticles(); ++i){
		for(auto& outMono : this->K1(i, delta)) ret.push_back(outMono.OrderCopy());
	}
	return ret;
}

std::vector<mono> mono::K2(const coeff_class delta) const{
	std::vector<mono> ret;
	for(unsigned int i = 0; i < NParticles(); ++i){
		for(auto& outMono : this->K2(i, delta)) ret.push_back(outMono.OrderCopy());
	}
	return ret;
}

std::vector<mono> mono::K3(const coeff_class delta) const{
	std::vector<mono> ret;
	for(unsigned int i = 0; i < NParticles(); ++i){
		for(auto& outMono : this->K3(i, delta)) ret.push_back(outMono.OrderCopy());
	}
	return ret;
}

// These three L are only consistent if you use the equations of motion!
std::array<mono,2> mono::L1(const unsigned int targetParticle) const{
	return std::array<mono,2>({{this->DerivPm(targetParticle).MultPt(targetParticle),
			this->DerivPt(targetParticle).MultPp(targetParticle)}});
}

std::array<mono,1> mono::L2(const unsigned int targetParticle) const{
	return std::array<mono,1>({{this->DerivPm(targetParticle).MultPm(targetParticle)}});
}

std::array<mono,1> mono::L3(const unsigned int targetParticle) const{
	return std::array<mono,1>({{this->DerivPt(targetParticle).MultPm(targetParticle)}});
}

std::vector<mono> mono::L1() const{
	std::vector<mono> ret;
	for(unsigned int i = 0; i < NParticles(); ++i){
		for(auto& outMono : this->L1(i)) ret.push_back(outMono.OrderCopy());
	}
	return ret;
}

std::vector<mono> mono::L2() const{
	std::vector<mono> ret;
	for(unsigned int i = 0; i < NParticles(); ++i){
		for(auto& outMono : this->L2(i)) ret.push_back(outMono.OrderCopy());
	}
	return ret;
}

std::vector<mono> mono::L3() const{
	std::vector<mono> ret;
	for(unsigned int i = 0; i < NParticles(); ++i){
		for(auto& outMono : this->L3(i)) ret.push_back(outMono.OrderCopy());
	}
	return ret;
}

