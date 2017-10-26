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
	return os << "\b }";
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
	os << "\b }";
	return os;
}

template<>
inline std::ostream& operator<<(std::ostream& os, const std::vector<char>& out){
	os << "{";
	for(auto& element : out){
		if(element >= 0) os << " ";
		os << static_cast<int>(element) << ",";
	}
	return os << "\b }";
}

template<typename T, int N>
inline std::ostream& operator<<(std::ostream& os, const std::array<T,N>& out){
	os << "{";
	for(auto& element : out){
		if(element >= 0) os << " ";
		os << element << ",";
	}
	return os << "\b }";
}

#endif
