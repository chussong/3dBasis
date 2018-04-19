#include "mono.hpp"

Mono::Mono(const std::vector<particle>& particles, const coeff_class& coeff): 
    coeff(coeff), particles(particles) {
    Order();
}

Mono::Mono(const std::vector<int>& pm, const std::vector<int>& pt,
		const coeff_class& coeff): coeff(coeff) {
    if(pm.size() != pt.size()){
        std::cerr << "Error: attempted to construct a monomial out of particle "
            << "data with different sizes: {" << pm.size() << "," << pt.size()
            << "}. It will be blank instead." << std::endl;
        return;
    }
    particles.resize(pm.size());
    for(auto i = 0u; i < pm.size(); ++i){
        particles[i].pm = pm[i];
        particles[i].pt = pt[i];
    }
    Order();
}

const char& Mono::Pm(const int i) const {
    return particles[i].pm;
}

// can't have non-const accessors since they can make a mono fail to be ordered
/*char& Mono::Pm(const int i){
    return particles[i].pm;
}*/

void Mono::ChangePm(const int i, const char newValue) {
    particles[i].pm = newValue;
    Order();
}

const char& Mono::Pt(const int i) const {
    return particles[i].pt;
}

// can't have non-const accessors since they can make a mono fail to be ordered
/*char& Mono::Pt(const int i){
    return particles[i].pt;
}*/

void Mono::ChangePt(const int i, const char newValue) {
    particles[i].pt = newValue;
    Order();
}

// note: like all operations on completed Monos, this assumes both are ordered
bool Mono::operator==(const Mono& other) const {
    if(particles.size() != other.particles.size()){
        std::cerr << "Error: asked to compare two monomials with different "
            << "numbers of particles." << std::endl;
        return false;
    }
    for(unsigned int i = 0; i < particles.size(); ++i){
        if(particles[i].pm != other.particles[i].pm
                    || particles[i].pt != other.particles[i].pt) return false;
    }
    return true;
}

std::ostream& operator<<(std::ostream& os, const Mono& out) {
    return os << out.HumanReadable();
}

std::string MathematicaOutput(const Mono& out) {
    std::stringstream ss;
    if (std::abs<builtin_class>(out.coeff - 1) < EPSILON) {
        ss << out.particles;
    } else {
        ss << out.coeff << " * " << out.particles;
    }
    return ss.str();
}

std::string Mono::HumanReadable() const {
    std::ostringstream os;
    // WARNING: this assumes that the sign of the coefficient will be accounted
    // for in the function calling this! Only the absolute value is attached!
    if(std::abs<builtin_class>(Coeff() - 1) > EPSILON) {
        os << std::abs<builtin_class>(Coeff()) << "*{";
    }
    for(auto& p : particles) {
        if(p.pm != 0){
            os << "M";
            if(p.pm != 1) os << "^" << std::to_string(p.pm);
        }
        if(p.pt != 0){
            os << "T";
            if(p.pt != 1) os << "^" << std::to_string(p.pt);
        }
        os << "O"; // FIXME?? it would be nice if this were Φ
    }
    if(std::abs<builtin_class>(Coeff() - 1) > EPSILON) os << "}";
    return os.str();
}

Mono Mono::operator-() const {
    Mono ret(*this);
    ret.coeff *= -1;
    return ret;
}

int Mono::TotalPm() const {
    int total = 0;
    for(auto& p : particles) total += p.pm;
    return total;
}

int Mono::TotalPt() const {
    int total = 0;
    for(auto& p : particles) total += p.pt;
    return total;
}

// if the Mono is ordered, particles[0] is guaranteed to have the highest Pm
int Mono::MaxPm() const {
    if(particles.size() < 1) return -1;
    return particles[0].pm;
}

// for this, we have to actually check
int Mono::MaxPt() const {
    int max = -1;
    for(auto& p : particles) if(p.pt > max) max = p.pt;
    return max;
}

// return a vector containing one entry per distinguishable particle in *this.
// Each entry is the number of those particles contained.
std::vector<size_t> Mono::CountIdentical() const {
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
std::vector<size_t> Mono::PermutationVector() const {
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

bool Mono::ParticlePrecedence(const particle& a, const particle& b) {
    if(a.pm != b.pm) return a.pm > b.pm;
    return a.pt > b.pt;
}

void Mono::Order() {
    std::sort(particles.begin(), particles.end(), ParticlePrecedence);
}

Mono Mono::OrderCopy() const {
    Mono ret(*this);
    ret.Order();
    return ret;
}

std::vector<int> Mono::IdentifyNodes() const {
    return ::IdentifyNodes(particles);
}

std::vector<int> Mono::IdentifyPmNodes() const {
    return ::IdentifyNodes([this](unsigned int i){return Pm(i);}, NParticles());
}

std::vector<int> Mono::IdentifyPtNodes() const {
    return ::IdentifyNodes([this](unsigned int i){return Pt(i);}, NParticles());
}

// returns whether or not this particle will have a zero norm due to being a
// p-perp descendant; for 2-particle states, this is just everything where both
// Pm are equal. For higher particle numbers, it's states where each particle is
// identical except exactly one has 1 higher Pt than the others.
bool Mono::IsNull() const {
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

// NOTE! These break ordering, so you have to re-order when you're done!
Mono Mono::DerivPm(const unsigned int part) const {
    if(part >= particles.size()){
        std::cerr << "Error: monomial told to take a derivative of momentum Pm["
            << part << "], but it only knows about " << NParticles() << "."
            << std::endl;
        return *this;
    }
    Mono ret(*this);
    ret.coeff *= ret.Pm(part);
    ret.ChangePm(part, ret.Pm(part)-1);
    return ret;
}

Mono Mono::DerivPt(const unsigned int part) const {
    if(part >= particles.size()){
        std::cerr << "Error: monomial told to take a derivative of momentum Pt["
            << part << "], but it only knows about " << NParticles() << "."
            << std::endl;
        return *this;
    }
    Mono ret(*this);
    ret.coeff *= ret.Pt(part);
    ret.ChangePt(part, ret.Pt(part)-1);
    return ret;
}

std::vector<Mono> Mono::DerivPm() const {
    std::vector<Mono> ret;
    for(unsigned int particle = 0; particle < NParticles(); ++particle){
        ret.push_back(DerivPm(particle).OrderCopy());
        ret[particle].Order();
    }
    return ret;
}

std::vector<Mono> Mono::DerivPt() const {
    std::vector<Mono> ret;
    for(unsigned int particle = 0; particle < NParticles(); ++particle){
        ret.push_back(DerivPt(particle).OrderCopy());
        ret[particle].Order();
    }
    return ret;
}

Mono Mono::MultPm(const unsigned int particle) const {
    Mono ret(*this);
    ret.ChangePm(particle, ret.Pm(particle)+1);
    return ret;
}

Mono Mono::MultPt(const unsigned int particle) const {
    Mono ret(*this);
    ret.ChangePt(particle, ret.Pt(particle)+1);
    return ret;
}

Mono Mono::MultPp(const unsigned int particle) const {
    Mono ret(*this);
    ret.ChangePt(particle, ret.Pt(particle)+2);
    ret.ChangePm(particle, ret.Pm(particle)-1);
    ret.Coeff() /= 2;
    return ret;
}
