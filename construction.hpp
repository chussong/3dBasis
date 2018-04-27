#ifndef CONSTRUCTION_HPP
#define CONSTRUCTION_HPP

#include <array>

#include "io.hpp"

// T can be any type or class with an == operator; value indexes T by uint
template<typename Accessor>
inline std::vector<int> IdentifyNodes(Accessor A, const size_t size){
    std::vector<int> nodes;
    int newNode;
    newNode = 1;
    for(unsigned int i = 1; i < size; ++i){
        if(A(i) != A(i-1)){
            nodes.push_back(newNode);
            newNode = 1;
        } else {
            ++newNode;
        }
    }
    nodes.push_back(newNode);
    return nodes;
}

// should work for any containerish thing T with operator[]
template<typename T>
inline std::vector<int> IdentifyNodes(const T& container){
    return IdentifyNodes([container](unsigned int i){return container[i];},
                    container.size());
}

// default behavior for a vector of particles; used by Mono::IdentifyNodes()
template<>
inline std::vector<int> IdentifyNodes(const std::vector<particle>& particles){
    return IdentifyNodes([particles](unsigned int i){return 
                std::array<int, 2>({{particles[i].pm, particles[i].pt}}); },
                particles.size());
}

inline std::vector<std::vector<int>> Permute(const std::vector<std::vector<int>> ordered){
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

template<typename T>
inline std::vector<std::vector<T>> GetStatesByDegree(const T numP, 
		const T deg, const bool exact, const T min) {
    std::vector<std::vector<T>> ret;
    if(numP == 0 || deg == 0){
        ret.resize(1);
        ret[0].resize(numP, 0);
        return ret;
    }
    T lowerBound;
    if(exact){
        lowerBound = deg/numP;
        if(deg % numP == 0) lowerBound--;
    } else {
        lowerBound = -1;
    }
    lowerBound = std::max<T>(lowerBound, min/numP + (min%numP!=0) - 1 - (min<0));
    for (T i = deg; i > lowerBound; --i) {
        std::vector<std::vector<T>> candidates = GetStatesByDegree<T>(numP-1,
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

template<typename T>
inline std::vector<std::vector<T>> GetStatesUpToDegree(const T numP,
		const T deg, const T M = 0) {
    return GetStatesByDegree<T>(numP, deg, false, M);
}

template<typename T>
inline std::vector<std::vector<T>> GetStatesAtDegree(const T numP,
		const T deg, const T M = 0) {
    return GetStatesByDegree<T>(numP, deg, true, M);
}

inline std::vector<std::vector<particle>> CombinedCfgs(
    const std::vector<particle>& baseCfg,
    const std::vector<std::vector<int>>& newCfgs, const int componentToChange) {
    std::vector<std::vector<particle>> ret;
    if(componentToChange < 1 || componentToChange > 2){
        std::cerr << "Error: asked to change invalid component number "
            << componentToChange << "; this number must be in [1,2]." << std::endl;
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

inline std::vector<std::vector<int>> CfgsFromNodePartition(
		const std::vector<int> nodes,
		std::vector<std::vector<int>> nodeEnergies){
    if(nodeEnergies.size() == 0 || nodes.size() != nodeEnergies[0].size()
                    || nodes.size() == 0){
        std::cerr << "Error: asked to turn node configurations into regular "
            << "cfgs but the arguments were malformed:" << std::endl;
        // std::cerr << "nodes is the following: " << nodes << std::endl;
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
            nodePoss[i] = GetStatesAtDegree<int>(nodes[i], possibility[i]);
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

inline std::vector<std::vector<int>> CfgsFromNodes(const int remainingEnergy,
		const std::vector<int>& nodes, const bool exact) {
    std::vector<std::vector<int>> ret;
    std::vector<std::vector<int>> nodeEnergies(GetStatesByDegree<int>(nodes.size(),
                            remainingEnergy, exact, 0));
    nodeEnergies = Permute(nodeEnergies);
    std::vector<std::vector<int>> perpCfgs(CfgsFromNodePartition(nodes, 
                            nodeEnergies));
    return CfgsFromNodePartition(nodes, nodeEnergies);
}

#endif
