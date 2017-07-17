#ifndef IO_HPP
#define IO_HPP

#include "constants.hpp"

// stream output operator template for vectors
template<typename T>
inline std::ostream& operator<<(std::ostream& os, const std::vector<T>& out){
	os << "{";
	for(auto& element : out){
		if(element >= 0) os << " ";
		os << element << ",";
	}
	os << "\b }";
	return os;
}

// specialization of above template which "transposes" particle vectors
template<>
inline std::ostream& operator<<(std::ostream& os, const std::vector<particle>& out){
	os << "{";
	for(auto& p : out){
		if(p.pm >= 0) os << " ";
		os << p.pm << ",";
	}
	os << "\b }{";
	for(auto& p : out){
		if(p.pt >= 0) os << " ";
		os << p.pt << ",";
	}
	os << "\b }{";
	for(auto& p : out){
		if(p.pp >= 0) os << " ";
		os << p.pp << ",";
	}
	os << "\b }";
	return os;
}

#endif
