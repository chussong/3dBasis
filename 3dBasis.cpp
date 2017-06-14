#include "3dBasis.hpp"

int main(int argc, char* argv[]) {
	if(argc < 3 || argc > 4){
		std::cerr << "Error: provide number of particles and degree when "
			<< "calling." << std::endl;
		return EXIT_FAILURE;
	}
	std::string numP(argv[1]);
	std::string deg (argv[2]);

	if(argc == 3){
		basis startingBasis(std::stoi(numP), std::stoi(deg));
	}
	if(argc == 4){
		std::string M(argv[3]);
		basis startingBasis(std::stoi(numP), std::stoi(deg), std::stoi(M));
	}
	return EXIT_SUCCESS;
}

bool mono::operator==(const mono& other) const{
	if(pm.size() != other.pm.size() || pp.size() != other.pp.size()){
		std::cerr << "Error: attempted to compare two monomials with different "
			<< "numbers of particles." << std::endl;
		return false;
	}
	for(unsigned int i = 0; i < pm.size(); ++i){
		if(pm[i] != other.pm[i] || pp[i] != other.pp[i]) return false;
	}
	return true;
}

std::ostream& operator<<(std::ostream& os, const mono& out){
	os << "{" << out.pm << "," << out.pp << "}";
	return os;
}

mono mono::DerivPm(const unsigned int particle){
	if(particle >= pm.size()){
		std::cerr << "Error: monomial told to take a derivative of momentum Pm["
			<< particle << "], but it only knows about " << pm.size() << "."
			<< std::endl;
		return *this;
	}
	mono ret(*this);
	ret.coeff *= ret.pm[particle];
	ret.pm[particle]--;
	return ret;
}

mono mono::DerivPp(const unsigned int particle){
	if(particle >= pp.size()){
		std::cerr << "Error: monomial told to take a derivative of momentum Pp["
			<< particle << "], but it only knows about " << pp.size() << "."
			<< std::endl;
		return *this;
	}
	mono ret(*this);
	ret.coeff *= ret.pp[particle];
	ret.pp[particle]--;
	return ret;
}

std::vector<mono> mono::DerivPm(){
	std::vector<mono> ret;
	for(unsigned int particle = 0; particle < pm.size(); ++particle){
		ret.push_back(DerivPm(particle));
	}
	return ret;
}

std::vector<mono> mono::DerivPp(){
	std::vector<mono> ret;
	for(unsigned int particle = 0; particle < pp.size(); ++particle){
		ret.push_back(DerivPp(particle));
	}
	return ret;
}

basis::basis(const int numP, const int degree): degree(degree), numP(numP) {
	std::vector<std::vector<int>> minus = GetStatesUpToDegree(numP, degree);
	basisMonos = AddPlusUpToDegree(minus, degree);
}

basis::basis(const int numP, const int degree, const int M): degree(degree),
		numP(numP) {
	/*std::vector<std::vector<int>> minus = GetStatesUpToDegree(degree, numP, M);
	std::vector<mono> basisMonos = AddPlusUpToDegree(minus, degree, M);*/

	if(std::abs(M) > degree){
		std::cerr << "Error: there are no states with |M| > degree!" << std::endl;
		throw(std::runtime_error("basis constructor"));
	}
	std::vector<std::vector<int>> minus = GetStatesUpToDegree(numP, (degree+M)/2, M);
	basisMonos = AddPlusUpToDegree(minus, degree, M);
}

std::ostream& operator<<(std::ostream& os, const std::vector<int>& out){
	os << "{";
	for(auto& i : out){
		if(i >= 0) os << " ";
		os << i << ",";
	}
	os << "\b }";
	return os;
}

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
		const int deg, const bool strict, const int min){
	std::vector<std::vector<int>> ret;
	if(numP == 0 || deg == 0){
		ret.resize(1);
		ret[0].resize(numP, 0);
		return ret;
	}
	int lowerBound;
	if(strict){
		lowerBound = deg/numP;
		if(deg % numP == 0) lowerBound--;
	} else {
		lowerBound = -1;
	}
	lowerBound = std::max(lowerBound, min/numP + (min%numP!=0) - 1 - (min<0));
	for(int i = deg; i > lowerBound; --i){
		std::vector<std::vector<int>> candidates = GetStatesByDegree(numP-1,
				deg-i, strict, min-i);
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

std::vector<std::vector<int>> GetStatesUpToDegree(const int numP,
		const int deg, const int M){
	return GetStatesByDegree(numP, deg, false, M);
}

std::vector<std::vector<int>> GetStatesAtDegree(const int numP,
		const int deg){
	return GetStatesByDegree(numP, deg, true, 0);
}

std::vector<std::vector<int>> GetStatesAtDegree(const int numP,
		const int deg, const int M){
	return GetStatesByDegree(numP, deg, true, M);
}

// we only put P_+ onto particles with 0 P_- in order to avoid overcounting
std::vector<std::shared_ptr<mono>> basis::AddPlusUpToDegree(
		const std::vector<std::vector<int>>& minus, const int degree){
	std::vector<std::shared_ptr<mono>> ret;
	for(auto& entry : minus){
		int startingDeg = 0;
		int zeroCount = 0;
		for(auto& p : entry){
			startingDeg += p;
			if(p == 0) ++zeroCount;
		}
		std::vector<std::vector<int>> perp(GetStatesUpToDegree(zeroCount,
					degree - startingDeg));
		for(auto& entryPerp : perp){
			entryPerp.insert(entryPerp.begin(), entry.size()-zeroCount, 0);
			ret.push_back(std::make_shared<mono>(entry, entryPerp));
		}
	}
	ret = mono::FinishPerpAtDegree(ret, degree); // they are not valid monos without this!
	return ret;
}

std::vector<std::shared_ptr<mono>> basis::AddPlusUpToDegree(
		const std::vector<std::vector<int>>& minus, const int degree, const int M){
	std::vector<std::shared_ptr<mono>> ret;
	for(auto& entry : minus){
		int startingDeg = 0;
		int zeroCount = 0;
		for(auto& p : entry){
			startingDeg += p;
			if(p == 0) ++zeroCount;
		}
		std::vector<std::vector<int>> perp(GetStatesAtDegree(zeroCount,
					std::min(degree - startingDeg, startingDeg - M)));
		for(auto& entryPerp : perp){
			entryPerp.insert(entryPerp.begin(), entry.size()-zeroCount, 0);
			ret.push_back(std::make_shared<mono>(entry, entryPerp));
		}
	}
	ret = mono::FinishPerpAtDegree(ret, degree); // they are not valid monos without this!
	return ret;
}

// this function takes fake monos whose Pp field has been filled with P_+ and 
// converts them to proper ones with only Pm and Pp using the equation of motion
std::vector<std::shared_ptr<mono>> mono::FinishPerpAtDegree(
		const std::vector<std::shared_ptr<mono>>& combined,
		const int degree){
	std::cout << "Asked to finish the following monos:" << std::endl;
	std::vector<std::shared_ptr<mono>> ret;
	for(auto& mn : combined){
		std::cout << *mn << " -> " << std::endl;
		mono baseMono(*mn);
		int remainingEnergy = degree;
		for(unsigned int i = 0; i < baseMono.NParticles(); ++i){
			remainingEnergy -= mn->pm[i] + mn->pp[i];
			baseMono.pm[i] -= mn->pp[i];
			baseMono.pp[i] *= 2;
		}
		std::vector<int> nodes(baseMono.IdentifyNodes());
		std::vector<std::vector<int>> nodeEnergies(GetStatesAtDegree(nodes.size(),
					remainingEnergy));
		nodeEnergies = Permute(nodeEnergies);
		std::vector<std::vector<int>> perpCfgs(CfgsFromNodes(nodes, 
					nodeEnergies));
		std::shared_ptr<mono> newMono;
		for(auto& cfg : perpCfgs){
			newMono = std::make_shared<mono>(baseMono);
			std::cout << *newMono << " + " << cfg;
			for(auto i = 0; i < newMono->pp.size(); ++i){
				newMono->pp[i] += cfg[i];
			}
			std::cout << " = " << *newMono << std::endl;
			ret.push_back(newMono);
		}
		std::cout << "--------------------" << std::endl;
	}
	return ret;
}

std::vector<int> mono::IdentifyNodes() const{
	std::vector<int> nodes;
	int newNode;
	newNode = 1;
	for(unsigned int i = 1; i < NParticles(); ++i){
		if(pm[i] != pm[i-1] || pp[i] != pp[i-1]){
			nodes.push_back(newNode);
			newNode = 1;
		} else {
			++newNode;
		}
	}
	nodes.push_back(newNode);
	return nodes;
}

std::vector<std::vector<int>> mono::CfgsFromNodes(const std::vector<int> nodes,
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
		throw(std::runtime_error("CfgsFromNodes"));
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
			for(auto n = 0; n < nodePoss.size(); ++n){
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
