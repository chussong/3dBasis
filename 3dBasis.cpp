#include "3dBasis.hpp"

int main(int argc, char* argv[]) {
	if(argc < 3 || argc > 3){
		std::cerr << "Error: provide number of particles and degree when "
			<< "calling." << std::endl;
		return EXIT_FAILURE;
	}
	std::string numP(argv[1]);
	std::string deg (argv[2]);

	splitBasis startingBasis(std::stoi(numP), std::stoi(deg));
	std::cout << "Constructed a basis with " << startingBasis.EvenBasis().size()
		<< " even elements and " << startingBasis.OddBasis().size()
		<< " odd elements." << std::endl;
	basis targetBasis(std::stoi(numP), std::stoi(deg)-1);
	std::cout << "Constructed a target basis with " << targetBasis.size()
		<< " elements." << std::endl;
	// we could split the target basis by perp parity as well but it would
	// require constructing four matrices, since K_\perp takes (even)<->(odd)

	// - create matrix of K acting on each element of startingBasis
	Matrix evenKActions = KMatrix(startingBasis.EvenBasis(), targetBasis);
	Matrix oddKActions = KMatrix(startingBasis.OddBasis(), targetBasis);
	//std::cout << KActions << std::endl;

	// - find kernel of above matrix and output
	std::vector<poly> evenKernel = Kernel(evenKActions, startingBasis.EvenBasis());
	std::vector<poly> oddKernel = Kernel(oddKActions, startingBasis.OddBasis());

	/*std::cout << "Found the following " << kernel.size() << "-dimensional kernel:"
		<< std::endl;
	for(auto& kernelVec : kernel) std::cout << kernelVec << std::endl;*/
	
	return EXIT_SUCCESS;
}

mono::mono(const std::vector<int>& pm, const std::vector<int>& pt,
		const std::vector<int>& pp,	const coeff_class& coeff): coeff(coeff){
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
	return os << out.coeff << " * " << out.particles;
}

std::ostream& operator<<(std::ostream& os, const poly& out){
	for(auto& component : out) os << component << " + ";
	return out.size() > 0 ? os << "\b\b \b\b" : os;
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

/* The following two functions might actually work just as well as Order():
bool ParticlePrecedence(const particle& a, const particle& b){
	if(a.pm != b.pm) return a.pm > b.pm;
	if(a.pt != b.pt) return a.pt > b.pt;
	return a.pp > b.pp;
}

void mono::Order(){
	std::sort(particles.begin(), particles.end(), ParticlePrecedence(a,b));
}*/

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
}

mono mono::OrderCopy() const{
	mono ret(*this);
	ret.Order();
	return ret;
}

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
	++ret.Pp(particle);
	return ret;
}

std::array<mono, 4> mono::K1(const unsigned int particle) const{
	std::array<mono, 4> ret({{
			2*this->DerivPp(particle).DerivPp(particle).MultPp(particle),
			2*this->DerivPt(particle).DerivPp(particle).MultPt(particle),
			/*(technically 2*delta*)*/this->DerivPp(particle),
			this->DerivPt(particle).DerivPt(particle).MultPm(particle)}});
	return ret;
}

std::array<mono, 5> mono::K2(const unsigned int particle) const{
	std::array<mono, 5> ret({{
			-2*this->DerivPt(particle).DerivPp(particle).MultPp(particle),
			-2*this->DerivPt(particle).DerivPm(particle).MultPm(particle),
			-this->DerivPt(particle).DerivPt(particle).MultPt(particle),
			-/*(technically 2*delta*)*/this->DerivPt(particle),
			-2*this->DerivPm(particle).DerivPp(particle).MultPt(particle)}});
	return ret;
}

std::array<mono, 4> mono::K3(const unsigned int particle) const{
	std::array<mono, 4> ret({{
			2*this->DerivPm(particle).DerivPm(particle).MultPm(particle),
			2*this->DerivPt(particle).DerivPm(particle).MultPt(particle),
			/*(technically 2*delta*)*/this->DerivPm(particle),
			this->DerivPt(particle).DerivPt(particle).MultPp(particle)}});
	return ret;
}

std::vector<mono> mono::K1() const{
	std::vector<mono> ret;
	for(unsigned int i = 0; i < NParticles(); ++i){
		for(auto& outMono : this->K1(i)) ret.push_back(outMono.OrderCopy());
	}
	return ret;
}

std::vector<mono> mono::K2() const{
	std::vector<mono> ret;
	for(unsigned int i = 0; i < NParticles(); ++i){
		for(auto& outMono : this->K2(i)) ret.push_back(outMono.OrderCopy());
	}
	return ret;
}

std::vector<mono> mono::K3() const{
	std::vector<mono> ret;
	for(unsigned int i = 0; i < NParticles(); ++i){
		for(auto& outMono : this->K3(i)) ret.push_back(outMono.OrderCopy());
	}
	return ret;
}

poly::poly(std::vector<mono> terms){
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

poly poly::K1() const{
	poly ret;
	for(auto& term : terms){
		for(auto& newTerm : term.K1()) ret += newTerm;
	}
	return ret;
}

poly poly::K2() const{
	poly ret;
	for(auto& term : terms){
		for(auto& newTerm : term.K2()) ret += newTerm;
	}
	return ret;
}

poly poly::K3() const{
	poly ret;
	for(auto& term : terms){
		for(auto& newTerm : term.K3()) ret += newTerm;
	}
	return ret;
}

basis::basis(const int numP, const int degree) {
	// 1: generate all possibilities for P_-
	// 2: identify nodes in P_-, generate possible distributions of P_\perp to
	// 		the nodes and then P_\perp within each node, adding each at top level
	// 3: identify nodes in (P_- and P_\perp), repeating step 2 for P_+
	// 4: take list from step 3 and create a mono from each entry, then store
	// 		the list of these as basisMonos
	// NOTE: this does not make any use of the symmetries between the components
	// 		in constructing the basis; for instance, the states where every
	// 		pm = 0 are followed by a copy of the earlier parts of the basis
	std::vector<std::vector<int>> minus = GetStatesUpToDegree(numP, degree);
	std::vector<std::vector<particle>> particleCfgs;
	for(auto& minusCfg : minus){
		std::vector<particle> newCfg(minusCfg.size());
		for(auto i = 0u; i < newCfg.size(); ++i) newCfg[i].pm = minusCfg[i];
		particleCfgs.push_back(newCfg);
	}

	std::vector<int> nodes;
	std::vector<std::vector<particle>> newCfgs;
	for(auto& configuration : particleCfgs){
		nodes = IdentifyNodes(configuration);
		int remainingEnergy = degree;
		for(auto& part : configuration) remainingEnergy -= part.pm;
		std::vector<std::vector<int>> perp(CfgsFromNodes(remainingEnergy, nodes,
															false));
		for(auto& newCfg : CombinedCfgs(configuration, perp, 2)){
			newCfgs.push_back(newCfg);
		}
	}

	particleCfgs.clear();
	for(auto& cfg : newCfgs){
		nodes = IdentifyNodes(cfg);
		int remainingEnergy = degree;
		for(auto& part : cfg) remainingEnergy -= part.pm + part.pt;
		std::vector<std::vector<int>> plus(CfgsFromNodes(remainingEnergy, nodes,
															true));
		for(auto& newCfg : CombinedCfgs(cfg, plus, 3)){
			particleCfgs.push_back(newCfg);
		}
	}

	for(auto& cfg : particleCfgs){
		basisMonos.emplace_back(cfg);
		//std::cout << cfg << std::endl;
	}
}

std::vector<std::vector<particle>> basis::CombinedCfgs(
		const std::vector<particle>& baseCfg,
		const std::vector<std::vector<int>>& newCfgs, const int componentToChange){
	std::vector<std::vector<particle>> ret;
	if(componentToChange < 1 || componentToChange > 3){
		std::cerr << "Error: asked to change invalid component number "
			<< componentToChange << "; this number must be in [1,3]." << std::endl;
		return ret;
	}
	std::vector<particle> toAdd;
	for(auto& newCfg : newCfgs){
		if(newCfg.size() != baseCfg.size()){
			std::cerr << "Error: attempted to combine base cfg of size "
				<< baseCfg.size() << " with node data of size " << newCfg.size()
				<< "; these must be the same." << std::endl;
			throw(std::runtime_error("CombinedCfgs"));
		}
		toAdd = baseCfg;
		for(auto i = 0u; i < toAdd.size(); ++i){
			if(componentToChange == 3){
				toAdd[i].pp = newCfg[i];
				continue;
			}
			if(componentToChange == 2){
				toAdd[i].pt = newCfg[i];
				continue;
			}
			if(componentToChange == 1){
				toAdd[i].pm = newCfg[i];
				continue;
			}
		}
		ret.push_back(toAdd);
	}
	return ret;
}

std::vector<std::vector<int>> basis::CfgsFromNodes(const int remainingEnergy,
		const std::vector<int>& nodes, const bool exact){
	std::vector<std::vector<int>> ret;
	std::vector<std::vector<int>> nodeEnergies(GetStatesByDegree(nodes.size(),
				remainingEnergy, exact, 0));
	nodeEnergies = Permute(nodeEnergies);
	std::vector<std::vector<int>> perpCfgs(CfgsFromNodePartition(nodes, 
				nodeEnergies));
	return CfgsFromNodePartition(nodes, nodeEnergies);
}

std::vector<std::vector<int>> basis::CfgsFromNodePartition(const std::vector<int> nodes,
		std::vector<std::vector<int>> nodeEnergies){
	if(nodeEnergies.size() == 0 || nodes.size() != nodeEnergies[0].size()
			|| nodes.size() == 0){
		std::cerr << "Error: asked to turn node configurations into regular "
			<< "cfgs but the arguments were malformed:" << std::endl;
		std::cerr << "nodes is the following: " << nodes << std::endl;
		std::cerr << "nodeEnergies has size " << nodeEnergies.size() << std::endl;
		if(nodeEnergies.size() > 0){
			std::cerr << " and its first entry is: ";
			for(auto& node: nodeEnergies[0]) std::cerr << node << ",";
			std::cerr << "\b " << std::endl;
		}
		throw(std::runtime_error("CfgsFromNodePartition"));
	}
	std::vector<std::vector<int>> finalCfgs;
	int nPart = 0;
	for(auto& nodeSize : nodes) nPart += nodeSize;
	//std::cout << "NODES:          " << nodes << std::endl;
	for(auto& possibility : nodeEnergies){
		// nodePoss layers: nodes->possible configurations->particle energies
		//std::cout << "NODE PARTITION: " << possibility << std::endl;
		std::vector<std::vector<std::vector<int>>> nodePoss;
		nodePoss.resize(nodes.size());
		for(unsigned int i = 0; i < nodes.size(); ++i){
			nodePoss[i] = GetStatesAtDegree(nodes[i], possibility[i]);
		}
		std::vector<unsigned int> cfgSpec;
		cfgSpec.resize(nodes.size(), 0);
		std::vector<int> cfg;
		cfg.resize(nPart);
		unsigned int p;
		bool done = false;
		while(!done){
			p = 0;
			for(auto n = 0u; n < nodePoss.size(); ++n){
				for(auto& partEnergy : nodePoss[n][cfgSpec[n]]){
					cfg[p] = partEnergy;
					++p;
				}
			}
			//std::cout << "Perp cfg constructed: " << cfg << std::endl;
			finalCfgs.push_back(cfg);
			for(unsigned int n = cfgSpec.size()-1; n < nodePoss.size(); --n){
				if(++cfgSpec[n] < nodePoss[n].size()){
					break;
				} else {
					cfgSpec[n] = 0;
					if(n == 0) done = true;
				}
			}
		}
	}
	return finalCfgs;
}

unsigned int basis::FindInBasis(const std::vector<int>& pm, const std::vector<int>& pt,
		const std::vector<int>& pp){
	return FindInBasis(mono(pm, pt, pp));
}

unsigned int basis::FindInBasis(const mono& wildMono){
	for(auto i = 0u; i < basisMonos.size(); ++i){
		if(basisMonos[i] == wildMono) return i;
	}
	std::cout << "Warning! Failed to find the following mono in our basis: "
		<< wildMono << std::endl;
	return -1u;
}

/*std::array<basis,2> basis::ParitySplit() const{
	std::array<basis,2> ret = {{*this, *this}};
	for(int i = 0; i < 2; ++i){
		ret[i].basisMonos.erase(std::remove_if(ret[i].begin(), ret[i].end(), 
				[i](mono& m){ return m.TotalPt()%2 == i; } ));
	}
	return ret;
}*/

void basis::DeleteOdd(){
	basisMonos.erase(std::remove_if(basisMonos.begin(), basisMonos.end(), 
				splitBasis::IsOdd), basisMonos.end());
}

void basis::DeleteEven(){
	basisMonos.erase(std::remove_if(basisMonos.begin(), basisMonos.end(),
				splitBasis::IsEven), basisMonos.end());
}

Triplet basis::ExpressMono(const mono& toExpress, const int column,
		const int rowOffset) const{
	for(auto i = 0u; i < basisMonos.size(); ++i){
		if(toExpress == basisMonos[i])
			return Triplet(rowOffset+i, column, toExpress.Coeff());
	}
	std::cerr << "Error: tried to express the monomial " << toExpress
	<< " on the given basis but was not able to identify it." << std::endl;
	return Triplet(-1, -1, toExpress.Coeff());
}

std::list<Triplet> basis::ExpressPoly(const poly& toExpress, 
		const int column, const int rowOffset) const{
	std::list<Triplet> ret;
	if(toExpress.size() == 0) return ret;
	unsigned int hits = 0u;
	for(auto i = 0u; i < basisMonos.size(); ++i){
		for(auto& term : toExpress){
			if(term == basisMonos[i]){
				ret.emplace_front(rowOffset+i, column, term.Coeff());
				++hits;
				break;
			}
		}
		if(hits == toExpress.size()) break;
	}
	if(hits < toExpress.size()){
		std::cerr << "Error: tried to express the polynomial " << toExpress
		<< " on the given basis but was only able to identify " << hits
		<< " of " << toExpress.size() << " terms." << std::endl;
	}
	//std::cout << "Expressed " << toExpress << " as the following triplets:" << std::endl;
	//for(auto& trip : ret) std::cout << trip << std::endl;
	return ret;
}

bool splitBasis::IsOdd(const mono& toTest){
	return toTest.TotalPt()%2 == 1;
}

bool splitBasis::IsEven(const mono& toTest){
	return !IsOdd(toTest);
}

// this could definitely be done more intelligently if speed were important
splitBasis::splitBasis(const int numP, const int degree): 
	evenBasis(numP, degree), oddBasis(numP, degree){
	evenBasis.DeleteOdd();
	oddBasis.DeleteEven();
}

std::list<Triplet> splitBasis::ExpressPoly(const poly& toExpress, 
		const int column, const int rowOffset) const{
	std::list<Triplet> ret;
	for(auto& term : toExpress){
		if(IsOdd(term)){
			ret.push_front(oddBasis.ExpressMono(term, column, rowOffset));
		} else {
			ret.push_front(evenBasis.ExpressMono(term, column, rowOffset));
		}
	}
	return ret;
}

// note: triplets displayed (row, column, value) despite matrix's storage type
std::ostream& operator<<(std::ostream& os, const Triplet& out){
	return os << "(" << out.row() << "," << out.col() << "," << out.value()
		<< ")";
}

/*
// note: I recommend starting a new line before this, but I'm not your parents
std::ostream& operator<<(std::ostream& os, const Matrix& out){
	out.isVector() ? os << "( " : os << "/ ";
	for(auto col = 0; col < out.cols(); ++col){
		os << out.coeff(0, col) << ", ";
	}
	out.isVector() ? return os << "\b\b )" : os << "\b\b \\\n";

	os << "| ";
	for(auto row = 1; row < out.rows()-1; ++row){
		for(auto col = 0; col < out.cols(); ++col){
			os << out.coeff(row, col) << ", ";
		}
	}
	os << "\b\b |\n";

	os << "\\ ";
	for(auto col = 0; col < out.cols(); ++col){
		of << out.coeff(out.rows()-1, col) << ", ";
	}
	return os << "\b\b /";
}*/

std::vector<std::vector<int>> Permute(const std::vector<std::vector<int>> ordered){
	std::vector<std::vector<int>> ret;
	for(auto& base : ordered){
		std::vector<int> current(base);
		std::sort(current.begin(), current.end(), std::greater<int>());
		do{
			ret.push_back(current);
		} while(std::prev_permutation(current.begin(), current.end()));
	}
	return ret;
}

std::vector<std::vector<int>> GetStatesByDegree(const int numP, 
		const int deg, const bool exact, const int min){
	std::vector<std::vector<int>> ret;
	if(numP == 0 || deg == 0){
		ret.resize(1);
		ret[0].resize(numP, 0);
		return ret;
	}
	int lowerBound;
	if(exact){
		lowerBound = deg/numP;
		if(deg % numP == 0) lowerBound--;
	} else {
		lowerBound = -1;
	}
	lowerBound = std::max(lowerBound, min/numP + (min%numP!=0) - 1 - (min<0));
	for(int i = deg; i > lowerBound; --i){
		std::vector<std::vector<int>> candidates = GetStatesByDegree(numP-1,
				deg-i, exact, min-i);
		for(auto& cand : candidates){
			if(cand.size() == 0 || cand[0] <= i){
				cand.insert(cand.begin(), i);
				ret.push_back(cand);
			}
		}
	}
	return ret;
}

std::vector<std::vector<int>> GetStatesUpToDegree(const int numP,
		const int deg){
	return GetStatesByDegree(numP, deg, false, 0);
}

std::vector<std::vector<int>> GetStatesAtDegree(const int numP,
		const int deg){
	return GetStatesByDegree(numP, deg, true, 0);
}

Matrix KMatrix(const basis& startingBasis, const basis& targetBasis){
	//std::cout << "Constructing polynomials of K actions..." << std::endl;
	std::vector<poly> K1Actions, K2Actions, K3Actions;
	for(auto& basisMono : startingBasis){
		K1Actions.emplace_back(basisMono.K1());
		K2Actions.emplace_back(basisMono.K2());
		K3Actions.emplace_back(basisMono.K3());
	}

	//std::cout << "Converting K actions to triplets..." << std::endl;
	std::list<Triplet> entries = ConvertToRows(K1Actions, targetBasis, 0);
	entries.splice(entries.end(), ConvertToRows(K2Actions, targetBasis, 1));
	entries.splice(entries.end(), ConvertToRows(K3Actions, targetBasis, 2));

	//std::cout << "List of triplets done, matrixifying..." << std::endl;
	Matrix ret(3*targetBasis.size(), startingBasis.size());
	// if this is slow, can do ret.reserve(3*numP) to speed it up
	ret.setFromTriplets(entries.begin(), entries.end());
	// matrix must be compressed here but setFromTriplets does it automatically
	return ret;
}

std::list<Triplet> ConvertToRows(const std::vector<poly>& polyForms, 
		const basis& targetBasis, const Eigen::Index rowOffset){
	std::list<Triplet> ret = targetBasis.ExpressPoly(polyForms[0], 0,
			rowOffset*targetBasis.size());
	for(auto i = 1u; i < polyForms.size(); ++i){
		ret.splice(ret.end(), targetBasis.ExpressPoly(polyForms[i], i,
				rowOffset*targetBasis.size()));
	}
	return ret;
}

// takes QR decomposition of the matrix and returns the polynomial forms of its
// rightmost N columns, which are the N orthonormal basis vectors of the kernel
std::vector<poly> Kernel(const Matrix& KActions, const basis& startBasis){
	std::cout << "Computing kernel from K matrix..." << std::endl;
	QRSolver solver;
	solver.compute(KActions.transpose());
	std::cout << "Solved. Found rank " << solver.rank() << ", i.e. "
		<< startBasis.size() - solver.rank() << " primaries." << std::endl;

	std::cout << "Converting the kernel to polynomials..." << std::endl;

	/*DMatrix projector(startBasis.size(), startBasis.size() - solver.rank());
	for(auto row = 0; row < projector.rows(); ++row){
		for(auto col = 0; col < projector.cols(); ++col){
			if(row == col + solver.rank()){
				projector(row, col) = 1;
			} else {
				projector(row, col) = 0;
			}
		}
	}
	DMatrix kernelMatrix = solver.matrixQ()*projector;

	std::vector<poly> ret;
	std::cout << "Solving done: kernel matrix is the following " 
		<< kernelMatrix.rows() << "x" << kernelMatrix.cols() << " matrix, equal"
		<< " to the rank " << solver.rank() << ":\n" << kernelMatrix << std::endl;
	for(auto col = 0; col < kernelMatrix.cols(); ++col){
		ret.push_back(ColumnToPoly(kernelMatrix, col, startBasis));
	}*/

	DVector projector = Eigen::VectorXd::Zero(startBasis.size());
	DVector kernelVector(startBasis.size());
	std::vector<poly> ret;
	ret.resize(startBasis.size() - solver.rank());
	for(auto col = 0u; col < startBasis.size() - solver.rank(); ++col){
		projector(solver.rank() + col-1) = 0;
		projector(solver.rank() + col) = 1;
		kernelVector = solver.matrixQ()*projector;
		//std::cout << kernelVector << "\n----------" << std::endl;
		ret[col] = VectorToPoly(kernelVector, startBasis);
	}

	return ret;
}

/*poly VectorToPoly(const Vector& kernelVector, const basis& startBasis){
	poly ret;
	if(static_cast<size_t>(kernelMatrix.rows()) != startBasis.size()){
		std::cerr << "Error: the given Q matrix has " << kernelMatrix.rows()
			<< " rows, " << "but the given basis has " << startBasis.size() 
			<< " monomials. These must be the same." << std::endl;
		return ret;
	}
	for(auto row = 0; row < kernelVector.rows(); ++row){
		if(kernelVector.coeff(row, col) == 0) continue;
		ret += kernelMatrix.coeff(row, col)*startBasis[row];
	}

	if(ret.size() == 0) return ret;
	coeff_class smallestCoeff = std::abs(ret[0].Coeff());
	for(auto& term : ret) smallestCoeff = std::min(std::abs(term.Coeff()), smallestCoeff);
	for(auto& term : ret) term /= smallestCoeff;
	return ret;
}*/

poly VectorToPoly(const DVector& kernelVector, const basis& startBasis){
	poly ret;
	if(static_cast<size_t>(kernelVector.rows()) != startBasis.size()){
		std::cerr << "Error: the given Q column has " << kernelVector.rows()
			<< " rows, " << "but the given basis has " << startBasis.size() 
			<< " monomials. These must be the same." << std::endl;
		return ret;
	}
	for(auto row = 0; row < kernelVector.rows(); ++row){
		if(kernelVector.coeff(row) == 0) continue;
		ret += kernelVector.coeff(row)*startBasis[row];
	}

	if(ret.size() == 0) return ret;
	coeff_class smallestCoeff = std::abs(ret[0].Coeff());
	for(auto& term : ret) smallestCoeff = std::min(std::abs(term.Coeff()), smallestCoeff);
	for(auto& term : ret) term /= smallestCoeff;
	return ret;
}

poly ColumnToPoly(const Matrix& kernelMatrix, const Eigen::Index col, 
		const basis& startBasis){
	poly ret;
	if(static_cast<size_t>(kernelMatrix.rows()) != startBasis.size()){
		std::cerr << "Error: the given Q matrix has " << kernelMatrix.rows()
			<< " rows, " << "but the given basis has " << startBasis.size() 
			<< " monomials. These must be the same." << std::endl;
		return ret;
	}
	for(Eigen::Index row = 0; row < kernelMatrix.rows(); ++row){
		if(kernelMatrix.coeff(row, col) == 0) continue;
		ret += kernelMatrix.coeff(row, col)*startBasis[row];
	}

	if(ret.size() == 0) return ret;
	coeff_class smallestCoeff = std::abs(ret[0].Coeff());
	for(auto& term : ret) smallestCoeff = std::min(std::abs(term.Coeff()), smallestCoeff);
	for(auto& term : ret) term /= smallestCoeff;
	return ret;
}

poly ColumnToPoly(const DMatrix& kernelMatrix, const Eigen::Index col, 
		const basis& startBasis){
	poly ret;
	if(static_cast<size_t>(kernelMatrix.rows()) != startBasis.size()){
		std::cerr << "Error: the given Q matrix has " << kernelMatrix.rows()
			<< " rows, " << "but the given basis has " << startBasis.size() 
			<< " monomials. These must be the same." << std::endl;
		return ret;
	}
	for(Eigen::Index row = 0; row < kernelMatrix.rows(); ++row){
		if(kernelMatrix.coeff(row, col) == 0) continue;
		ret += kernelMatrix.coeff(row, col)*startBasis[row];
	}

	if(ret.size() == 0) return ret;
	coeff_class smallestCoeff = std::abs(ret[0].Coeff());
	for(auto& term : ret) smallestCoeff = std::min(std::abs(term.Coeff()), smallestCoeff);
	for(auto& term : ret) term /= smallestCoeff;
	return ret;
}

