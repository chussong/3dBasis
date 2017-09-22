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

const int& mono::Pm(const int i) const{
	return particles[i].pm;
}

int& mono::Pm(const int i){
	return particles[i].pm;
}

const int& mono::Pt(const int i) const{
	return particles[i].pt;
}

int& mono::Pt(const int i){
	return particles[i].pt;
}

const int& mono::Pp(const int i) const{
	if(usingEoM){
		std::cerr << "Warning: someone has requested Pp(" << i << ") from a "
			<< "monomial with the EoM enabled." << std::endl;
	}
	return particles[i].pp;
}

int& mono::Pp(const int i){
	if(usingEoM){
		std::cerr << "Warning: someone has requested Pp(" << i << ") from a "
			<< "monomial with the EoM enabled." << std::endl;
	}
	return particles[i].pp;
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
	if(std::abs(out.coeff - 1) < EPSILON) return os << out.particles;
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
	// WARNING: this assumes that the sign of the coefficient will be accounted
	// for in the function calling this! Only the absolute value is attached!
	if(std::abs(Coeff() - 1) > EPSILON) os << std::abs(Coeff()) << "*{";
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
	if(std::abs(Coeff() - 1) > EPSILON) os << "}";
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

// return a vector containing one entry per distinguishable particle in *this.
// Each entry is the number of those particles contained.
std::vector<size_t> mono::CountIdentical() const{
	std::vector<size_t> ret;
	std::vector<bool> counted(NParticles(), false);
	for(auto i = 0u; i < NParticles(); ++i){
		if(counted[i]) continue;
		int count = 1;
		for(auto j = i+1; j < NParticles(); ++j){
			if(particles[i] == particles[j]){
				++count;
				counted[j] = true;
			}
		}
		ret.push_back(count);
	}
	return ret;
}

// return a sorted vector with one entry per particle, each entry being the 
// location within this monomial of the FIRST particle which is identical to it
std::vector<size_t> mono::PermutationVector() const {
	std::vector<size_t> ret;
	std::vector<bool> counted(NParticles(), false);
	for(auto i = 0u; i < NParticles(); ++i){
		if(counted[i]) continue;
		ret.push_back(i);
		for(auto j = i+1; j < NParticles(); ++j){
			if(particles[i] == particles[j]){
				counted[j] = true;
				ret.push_back(i);
			}
		}
	}
	return ret;
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
	// this variant eliminates Pp completely, turning it into Pt^2/(2Pm)
	if(usingEoM){
		for(auto& p : particles){
			if(p.pp != 0){
				//std::cout << "Reordering " << HumanReadable();
				coeff /= std::pow(2,p.pp); // no bit shift in case coeff < 0
				p.pm -= p.pp;
				p.pt += 2*p.pp;
				p.pp = 0;
				//std::cout << " to " << HumanReadable() << "." << std::endl;
			}
		}
		//NullIfIllegal(); // I'm suspicious about this; should we be doing it?
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
	mono mirror(original.MirrorPM());
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

mono mono::MirrorPM() const{
	mono ret(*this);
	for(auto i = 0u; i < NParticles(); ++i){
		ret.Pm(i) = this->Pp(i);
		ret.Pp(i) = this->Pm(i);
		/*if(ret.Pp(i) < 0){
			ret.Pm(i) -= ret.Pp(i);
			ret.Pp(i) = 0;
		}*/
	}
	ret.Order();
	return ret;
}

// returns whether or not this is a Dirichlet-compliant monomial, i.e. that each
// particle has a Pm
bool mono::IsDirichlet() const{
	for(auto& p : particles){
		if(p.pm < 1){
			/*std::cout << "Found that " << HumanReadable() << " is non-Dirichlet"
				<< "." << std::endl;*/
			return false;
		}
	}
	return true;
}

// returns whether or not this particle will have a zero norm due to being a
// p-perp descendant; for 2-particle states, this is just everything where both
// Pm are equal. For higher particle numbers, it's states where each particle is
// identical except exactly one has 1 higher Pt than the others.
bool mono::IsNull() const {
	if (NParticles() < 2 ) {
		std::cerr << "Error: tried to call IsNull() on a monomial with only " 
			<< NParticles() << " particles." << std::endl;
		return false;
	}
	if (NParticles() == 2) return (Pm(0) == Pm(1)) && (TotalPt()%2 == 1);

	for (auto i = 1u; i < NParticles(); ++i) {
		if (Pm(i) != Pm(0) || Pt(i) != Pt(0) - 1) return false;
	}
	return true;
}

/*void mono::NullIfIllegal(){
	for(particle& p : particles){
		if(2*p.pm + p.pt < 0){
			// this is not actually a legal state, just throw it out?
			coeff = 0;
		}
	}
}*/

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
		ret.Pt(particle) += 2;
		ret.Pm(particle) -= 1;
		ret.Coeff() /= 2;
	} else {
		ret.Pp(particle) += 1;
	}
	return ret;
}

std::array<mono, 4> mono::K1(const unsigned int particle,
		const coeff_class delta) const{
	return std::array<mono, 4> ({{
			2*this->DerivPp(particle).DerivPp(particle).MultPp(particle),
			2*this->DerivPt(particle).DerivPp(particle).MultPt(particle),
			2*delta*this->DerivPp(particle),
			this->DerivPt(particle).DerivPt(particle).MultPm(particle)}});
}

std::array<mono, 5> mono::K2(const unsigned int particle,
		const coeff_class delta) const{
	return std::array<mono, 5> ({{
			-2*this->DerivPt(particle).DerivPp(particle).MultPp(particle),
			-2*this->DerivPt(particle).DerivPm(particle).MultPm(particle),
			-this->DerivPt(particle).DerivPt(particle).MultPt(particle),
			-2*delta*this->DerivPt(particle),
			-2*this->DerivPm(particle).DerivPp(particle).MultPt(particle)}});
}

std::array<mono, 4> mono::K3(const unsigned int particle,
		const coeff_class delta) const{
	return std::array<mono, 4> ({{
			2*this->DerivPm(particle).DerivPm(particle).MultPm(particle),
			2*this->DerivPt(particle).DerivPm(particle).MultPt(particle),
			2*delta*this->DerivPm(particle),
			this->DerivPt(particle).DerivPt(particle).MultPp(particle)}});
}

std::array<mono, 1> mono::K1_EoM(const unsigned int particle,
		const coeff_class) const{
	return std::array<mono, 1> ({{
			this->DerivPt(particle).DerivPt(particle).MultPm(particle)}});
}

std::array<mono, 3> mono::K2_EoM(const unsigned int particle,
		const coeff_class delta) const{
	return std::array<mono, 3> ({{
			-2*this->DerivPm(particle).DerivPt(particle).MultPm(particle),
			-this->DerivPt(particle).DerivPt(particle).MultPt(particle),
			-2*delta*this->DerivPt(particle)}});
}

std::array<mono, 4> mono::K3_EoM(const unsigned int particle,
		const coeff_class delta) const{
	return std::array<mono, 4> ({{
			2*this->DerivPm(particle).DerivPm(particle).MultPm(particle),
			2*this->DerivPm(particle).DerivPt(particle).MultPt(particle),
			2*delta*this->DerivPm(particle),
			this->DerivPt(particle).DerivPt(particle).MultPp(particle)}});
}

std::vector<mono> mono::K1(const coeff_class delta) const{
	std::vector<mono> ret;
	for(unsigned int i = 0; i < NParticles(); ++i){
		if(usingEoM){
			auto oneParticleK = K1_EoM(i, delta);
			for(auto& outMono : oneParticleK) ret.push_back(outMono.OrderCopy());
		} else {
			auto oneParticleK = K1(i, delta);
			for(auto& outMono : oneParticleK) ret.push_back(outMono.OrderCopy());
		}
	}
	//std::cout << "Action of K1 on " << HumanReadable() << ": " << ret << std::endl;
	return ret;
}

std::vector<mono> mono::K2(const coeff_class delta) const{
	std::vector<mono> ret;
	for(unsigned int i = 0; i < NParticles(); ++i){
		if(usingEoM){
			auto oneParticleK = K2_EoM(i, delta);
			for(auto& outMono : oneParticleK) ret.push_back(outMono.OrderCopy());
		} else {
			auto oneParticleK = K2(i, delta);
			for(auto& outMono : oneParticleK) ret.push_back(outMono.OrderCopy());
		}
	}
	//std::cout << "Action of K2 on " << HumanReadable() << ": " << ret << std::endl;
	return ret;
}

std::vector<mono> mono::K3(const coeff_class delta) const{
	std::vector<mono> ret;
	for(unsigned int i = 0; i < NParticles(); ++i){
		if(usingEoM){
			auto oneParticleK = K3_EoM(i, delta);
			for(auto& outMono : oneParticleK) ret.push_back(outMono.OrderCopy());
		} else {
			auto oneParticleK = K3(i, delta);
			for(auto& outMono : oneParticleK) ret.push_back(outMono.OrderCopy());
		}
	}
	//std::cout << "Action of K3 on " << HumanReadable() << ": " << ret << std::endl;
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

bool mono::IPPermutationCheck(const std::vector<int>& vec){
	for(auto i = 0u; i < vec.size()-1; i += 2){
		if(vec[i] > vec[i+1]) return false;
	}
	for(auto i = 0u; i < vec.size()-2; i += 2){
		if(vec[i] > vec[i+2]) return false;
	}
	return true;
}

std::tuple<coeff_class, int, int> mono::IPPairData(const mono& A, const mono& B,
		const unsigned int index1, const unsigned int index2){
	const bool index1InB = index1 >= A.NParticles();
	const bool index2InB = index2 >= A.NParticles();
	const int a1 = index1InB ? B.Pm(index1 - A.NParticles()) : A.Pm(index1);
	const int a2 = index2InB ? B.Pm(index2 - A.NParticles()) : A.Pm(index2);
	const int b1 = index1InB ? B.Pt(index1 - A.NParticles()) : A.Pt(index1);
	const int b2 = index2InB ? B.Pt(index2 - A.NParticles()) : A.Pt(index2);

	coeff_class coeff {1};
	const int totalPm {a1+a2};
	const int totalPt {b1+b2};

	if((a1 + b2)%2 == 1) coeff = -coeff;
	coeff *= Pochhammer(0.5, totalPm);
	coeff *= Pochhammer(2*totalPm + 1, totalPt);
	coeff *= 1 << totalPm;

	return std::make_tuple(coeff, totalPm, totalPt);
}

coeff_class mono::InnerProduct(const mono& A, const mono& B, 
		const GammaCache& cache, const KVectorCache& kCache){
	return IPZuhair(A, B, cache, kCache);
}

// Return true if a full cycle has been completed, otherwise return false.
// If we were using degeneracies we could reset kVector[i] = kVector[i-1], maybe
bool mono::IPMIncrementK(std::vector<int>& kVector, const std::vector<int>& totalPt){
	if(kVector.size() != totalPt.size()){
		std::cout << "kVector: " << kVector << "; totalPt: " << totalPt 
			<< std::endl;
		throw std::runtime_error("IP: kVector & totalPt size mismatch");
	}
	++kVector.back();
	for(unsigned int i = kVector.size()-1; i < kVector.size(); --i){
		if(kVector[i] > totalPt[i]/2){ // maybe (totalPt[i]+1)/2?
			if(i == 0) return true;
			++kVector[i-1];
			kVector[i] = 0;
		}
	}
	return false;
}

// this attempts to implement Matt's formulation of the inner product
coeff_class mono::IPMatt(const mono& A, const mono& B){
	// if perp parity doesn't match, it has to be zero
	if((A.TotalPt() + B.TotalPt()) % 2 == 1) return 0;

	std::cout << "Inner product of " << A << " with " << B << ":" << std::endl;

	unsigned int totalParticles = A.NParticles() + B.NParticles();
	// generate all permutations of (1,2,3,4,...,N) where two adjacent elements
	// represent a contraction, e.g. 1 with 2, 3 with 4, etc. Avoid overcounting
	// by demanding that each pair (i,j) have i < j and that i1 < i2 < i3 < ...
	
	// This is the dumbest possible way to do that.
	std::vector<std::vector<int>> permutations;
	std::vector<int> candidate;
	for(auto i = 0u; i < totalParticles; ++i) candidate.push_back(i);
	do {
		if(IPPermutationCheck(candidate)) permutations.push_back(candidate);
	} while(std::next_permutation(candidate.begin(), candidate.end()));

	// This method does not take advantage of the numerous degeneracies.

	coeff_class total {0};
	std::vector<int> kVector;
	for(auto& p : permutations){
		coeff_class coeff = 1;
		std::vector<int> totalPm;
		std::vector<int> totalPt;
		for(auto i = 0u; i < p.size() - 1; i += 2){
			auto pairData = IPPairData(A, B, p[i], p[i+1]);
			coeff *= std::get<0>(pairData);
			totalPm.push_back(std::get<1>(pairData));
			totalPt.push_back(std::get<2>(pairData));
		}
		std::cout << "Permutation: " << p << std::endl;
		std::cout << "Total Pm: " << totalPm << std::endl;
		std::cout << "Total Pt: " << totalPt << std::endl;
		kVector.clear();
		kVector.resize(p.size()/2, 0);
		coeff_class kSum {0};
		do{
			coeff_class kCoeff {1};
			int a {0};
			int b {0};
			int c {0};
			int n {0};
			for(auto i = 0u; i < kVector.size(); ++i){
				kCoeff *= Binomial(totalPt[i], 2*kVector[i]);
				kCoeff *= Pochhammer(0.5, kVector[i]);
				kCoeff /= Pochhammer(totalPm[i] + 1, kVector[i]);
				kCoeff *= 1 << kVector[i];

				a += totalPm[i] + kVector[i];
				b += kVector[i];
				c += totalPt[i]/2 - kVector[i];
				n += 1;
			}
			std::cout << "(a, b, c, n): (" << a << ", " << b << ", " << c
				<< ", " << n << ")" << std::endl;
			kCoeff *= IPFourier(a, b, c, n);
			std::cout << "Fourier result: " << IPFourier(a, b, c, n) << std::endl;
			std::cout << "Final kCoeff: " << kCoeff << std::endl;
			kSum += kCoeff;
		} while(!IPMIncrementK(kVector, totalPt));
		std::cout << "Coefficient of this permutation: " << coeff << std::endl;
		std::cout << "Final kSum: " << kSum << std::endl;
		std::cout << "Total from this permutation: " << coeff*kSum << std::endl;
		total += coeff * kSum;
	}


	return total;
}

// this is (5.9) from Matt's Revisiting3D.pdf except for the \mu and P_- terms
coeff_class mono::IPFourier(const int a, const int b, const int c, const int n){
	coeff_class A = static_cast<coeff_class>(a);
	coeff_class B = static_cast<coeff_class>(b);
	coeff_class C = static_cast<coeff_class>(c);
	coeff_class N = (static_cast<coeff_class>(n)-1)/2;
	coeff_class ret = 1;
	if(c % 2 == 1) ret = -ret;
	ret *= 1 << (a+2);	// a.k.a. 2^(a+2)
	ret *= M_PI * M_PI;
	ret *= Pochhammer(C + 1, c);
	ret *= Pochhammer(A + B + C + N, c);

	coeff_class logGamma = 2*std::lgamma(A + B + C + N);
	logGamma -= std::lgamma(A + C + N);
	logGamma -= std::lgamma(B + C + N);
	logGamma -= std::lgamma(2*A + 2*B + 4*C + 2*N);
	ret *= std::exp(logGamma);
	return ret;
}

// this attempts to implement Zuhair's formulation of the inner product. Note:
// it requires both monomials to have the same number of particles.
coeff_class mono::IPZuhair(const mono& A, const mono& B, 
		const GammaCache& cache, const KVectorCache& kCache){

	constexpr bool debug = false;

	int totalPm = 0;
	int totalPt = 0;
	for(auto i = 0u; i < A.NParticles(); ++i){
		totalPm += A.Pm(i);
		totalPt += A.Pt(i);
	}
	for(auto i = 0u; i < B.NParticles(); ++i){
		totalPm += B.Pm(i);
		totalPt += B.Pt(i);
	}
	if(A.NParticles() != B.NParticles() || totalPt % 2 == 1) return 0;

	coeff_class multiplicity = 1;
	for(auto& count : B.CountIdentical()) multiplicity *= Factorial(count);

	if(debug){
	std::cout << "--------------------" << std::endl;
	std::cout << "Inner product " << A << " x " << B << std::endl;
	std::cout << "Prefactor: " << cache.Prefactor(totalPm, totalPt) <<std::endl;
	std::cout << "Multiplicity: " << multiplicity << std::endl;
	}

	//std::vector<size_t> perm(B.NParticles());
	//for(auto i = 0u; i < B.NParticles(); ++i) perm[i] = i;
	std::vector<size_t> perm(B.PermutationVector());

	coeff_class sum = 0;
	std::vector<char> ptVector(A.NParticles());
	std::vector<std::array<char,2>> contractions(A.NParticles());;
	for(int totalK = 0; totalK <= totalPt/2; ++totalK){
		coeff_class totalKPrefactor = cache.Middle(totalPm, totalPt, totalK);
		if(debug){
		std::cout << "Middle coefficient @ (" << totalPm << ", " << totalPt
			<< ", " << totalK << "): " << totalKPrefactor << std::endl;
		}

		// - permute monomial B
		// - sum over all kVectors whose total is totalK and where each entry
		//   kVector[i] <= totalPt[i] (noting that totalPt[i] changes with perm)

		do{ // for each permutation of B
			// for each kVector whose total is totalK
			if(debug) std::cout << "Permutation: " << perm << std::endl;
			/*for(auto i = 0u; i < ptVector.size(); ++i){
				ptVector[i] = static_cast<char>((A.Pt(i) + B.Pt(perm[i]))/2);
			}
			std::sort(ptVector.begin(), ptVector.end(), std::greater<char>());
			std::cout << "PtVector: " << ptVector << std::endl;
			const KVectorBundle& kVectors = kCache.FromPt(ptVector, totalK);
			for(size_t startPt = 0; startPt < kVectors.size(); startPt += A.NParticles()){
				if(debug){
					std::cout << "k vector: "
					<< std::vector<char>(kVectors.begin() + startPt, 
							kVectors.begin() + startPt + A.NParticles())
					<< std::endl;
				}
				coeff_class logOfProduct = 0;
				for(auto i = 0u; i < A.NParticles(); ++i){
					logOfProduct += cache.Inner(A.Pm(i) + B.Pm(perm[i]),
							A.Pt(i) + B.Pt(perm[i]), kVectors[startPt + i]);
				}
				if(debug){
				std::cout << "Contribution from this permutation and kVector: "
					<< std::exp(logOfProduct) << std::endl;
				}
				sum += totalKPrefactor*std::exp(logOfProduct);
			}*/
			for(auto i = 0u; i < contractions.size(); ++i){
				contractions[i][0] = static_cast<char>(A.Pm(i) + B.Pm(perm[i]));
				contractions[i][1] = static_cast<char>(A.Pt(i) + B.Pt(perm[i]));
			}
			std::sort(contractions.begin(), contractions.end(),
					[](std::array<char,2> a, std::array<char,2> b){
						return a[1] > b[1];
					} );
			const KVectorBundle& kVectors = kCache.FromCont(contractions, totalK);
			for(size_t startPt = 0; startPt < kVectors.size(); startPt += A.NParticles()){
				if(debug){
					std::cout << "k vector: "
					<< std::vector<char>(kVectors.begin() + startPt, 
							kVectors.begin() + startPt + A.NParticles())
					<< std::endl;
				}
				coeff_class logOfProduct = 0;
				for(auto i = 0u; i < A.NParticles(); ++i){
					logOfProduct += cache.Inner(contractions[i][0],
							contractions[i][1], kVectors[startPt + i]);
				}
				if(debug){
				std::cout << "Contribution from this permutation and kVector: "
					<< std::exp(logOfProduct) << std::endl;
				}
				sum += totalKPrefactor*std::exp(logOfProduct);
			}
			/*for(const auto& kVector : VectorsAtK(totalK, perm, A, B, 0)){
				if(debug) std::cout << "k vector: " << kVector << std::endl;
				coeff_class logOfProduct = 0;
				for(auto i = 0u; i < A.NParticles(); ++i){
					logOfProduct += cache.Inner(A.Pm(i) + B.Pm(perm[i]),
							A.Pt(i) + B.Pt(perm[i]), kVector[i]);
				}
				if(debug){
				std::cout << "Contribution from this permutation and kVector: "
					<< std::exp(logOfProduct) << std::endl;
				}
				sum += totalKPrefactor*std::exp(logOfProduct);
			}*/
		}while(std::next_permutation(perm.begin(), perm.end()));
	}
	return A.Coeff()*B.Coeff()*cache.Prefactor(totalPm, totalPt)*multiplicity*sum;
}

// return a vector of all the vectors of length A.size() whose total is totalK,
// subject to the constraint that kVector[i] <= (A.Pt(i) + B.Pt(perm[i]))/2.
// This version operates by recursively calling itself on the remaining part
// of the vector until it reaches the end.
std::vector<std::vector<char>> mono::VectorsAtK(const char totalK, 
		const std::vector<size_t>& perm, const mono& A, const mono& B,
		const size_t start){
	//std::cout << "VectorsAtK(" << totalK << ", " << perm << ", " << A << ", "
		//<< B << ", " << start << "):" << std::endl;
	std::vector<std::vector<char>> ret;

	if(start >= A.NParticles()) return ret;
	const char maxK = (A.Pt(start) + B.Pt(perm[start]))/2;
	if(start == A.NParticles()-1 && totalK > maxK) return ret;
	if(start == A.NParticles()-1 && totalK <= maxK) return {{totalK}};

	for(int thisK = 0; thisK <= std::min(maxK, totalK); ++thisK){
		for(const auto& kv : VectorsAtK(totalK - thisK, perm, A, B, start+1)){
			ret.push_back(kv);
			ret.back().insert(ret.back().begin(), thisK);
		}
	}
	//for(auto& r : ret) std::cout << r << std::endl;
	return ret;
}

// return a vector of all the vectors of length A.size() whose total is totalK,
// subject to the constraint that kVector[i] <= (A.Pt(i) + B.Pt(perm[i]))/2.
// This version operates by starting at the beginning and adding new elements
// until it runs out of K.
// WARNING: THIS DOESN'T ACTUALLY WORK RIGHT NOW. I'M NOT PLANNING TO FINISH IT
// BECAUSE IT SEEMS TO BE NO FASTER THAN THE RECURSIVE ONE, WHICH WORKS ALREADY.
/*std::vector<std::vector<int>> mono::VectorsAtK(const int totalK,
		const std::vector<size_t>& perm, const mono& A, const mono& B,
		const size_t){
	std::vector<std::vector<int>> kVectors(1,std::vector<int>(1, totalK));
	std::vector<int> remainingK;
	int runningCount = 0;
	for(auto pos = A.NParticles()-1; pos < -1u; --pos){
		runningCount += (A.Pt(pos) + B.Pt(perm[pos]))/2;
		remainingK.insert(remainingK.begin(), runningCount);
	}
	for(auto pos = 0u; pos < A.NParticles() - 1; ++pos){
		std::vector<std::vector<int>> newKVectors;
		int maxK = (A.Pt(pos) + B.Pt(perm[pos]))/2;
		for(auto& partialVector : kVectors){
			for(int newK = 0; newK <= partialVector.back() && newK <= maxK
					&& partialVector.back() - newK < remainingK[pos]; ++newK){
				newKVectors.push_back(partialVector);
				newKVectors.back().back() = newK;
				newKVectors.back().push_back(partialVector.back() - newK);
			}
		}
		std::swap(newKVectors, kVectors);
	}
	//std::cout << kVectors.size() << " k vectors." << std::endl;
	return kVectors;
}*/

// just like an inner product, but including the intermediate operator "op" in
// its contractions. as before, all fields inside objects must be contracted 
// with a field inside a different object.
/*coeff_class mono::MatrixElement(const mono& A, const mono& B, const mono& op,
		const GammaCache& cache, const KVectorCache& kCache){

	constexpr bool debug = false;

	// generate (a series of) vector(s) of the Pm and Pt to be found in each
	// contraction

	int totalPm = 0;
	int totalPt = 0;
	for(auto i = 0u; i < A.NParticles(); ++i){
		totalPm += A.Pm(i);
		totalPt += A.Pt(i);
	}
	for(auto i = 0u; i < B.NParticles(); ++i){
		totalPm += B.Pm(i);
		totalPt += B.Pt(i);
	}
	if(A.NParticles() != B.NParticles() || totalPt % 2 == 1) return 0;

	coeff_class multiplicity = 1;
	for(auto& count : B.CountIdentical()) multiplicity *= Factorial(count);

	if(debug){
	std::cout << "--------------------" << std::endl;
	std::cout << "Inner product " << A << " x " << B << std::endl;
	std::cout << "Prefactor: " << cache.Prefactor(totalPm, totalPt) <<std::endl;
	std::cout << "Multiplicity: " << multiplicity << std::endl;
	}

	std::vector<size_t> perm(B.PermutationVector());

	coeff_class sum = 0;
	std::vector<char> ptVector(A.NParticles());
	std::vector<std::array<char,2>> contractions(A.NParticles());;
	for(int totalK = 0; totalK <= totalPt/2; ++totalK){
		coeff_class totalKPrefactor = cache.Middle(totalPm, totalPt, totalK);
		if(debug){
		std::cout << "Middle coefficient @ (" << totalPm << ", " << totalPt
			<< ", " << totalK << "): " << totalKPrefactor << std::endl;
		}

		// - permute monomial B
		// - sum over all kVectors whose total is totalK and where each entry
		//   kVector[i] <= totalPt[i] (noting that totalPt[i] changes with perm)

		do{ // for each permutation of B
			// for each kVector whose total is totalK
			if(debug) std::cout << "Permutation: " << perm << std::endl;
			for(auto i = 0u; i < contractions.size(); ++i){
				contractions[i][0] = static_cast<char>(A.Pm(i) + B.Pm(perm[i]));
				contractions[i][1] = static_cast<char>(A.Pt(i) + B.Pt(perm[i]));
			}
			std::sort(contractions.begin(), contractions.end(),
					[](std::array<char,2> a, std::array<char,2> b){
						return a[1] > b[1];
					} );
			const KVectorBundle& kVectors = kCache.FromCont(contractions, totalK);
			for(size_t startPt = 0; startPt < kVectors.size(); startPt += A.NParticles()){
				if(debug){
					std::cout << "k vector: "
					<< std::vector<char>(kVectors.begin() + startPt, 
							kVectors.begin() + startPt + A.NParticles())
					<< std::endl;
				}
				coeff_class logOfProduct = 0;
				for(auto i = 0u; i < A.NParticles(); ++i){
					logOfProduct += cache.Inner(contractions[i][0],
							contractions[i][1], kVectors[startPt + i]);
				}
				if(debug){
				std::cout << "Contribution from this permutation and kVector: "
					<< std::exp(logOfProduct) << std::endl;
				}
				sum += totalKPrefactor*std::exp(logOfProduct);
			}
		}while(std::next_permutation(perm.begin(), perm.end()));
	}
	return cache.Prefactor(totalPm, totalPt)*multiplicity*sum;
}*/
