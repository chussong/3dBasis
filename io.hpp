#ifndef IO_HPP
#define IO_HPP

#include <iostream>
#include "constants.hpp"

// stream output operator template for vectors
template<typename T>
inline std::ostream& operator<<(std::ostream& os, const std::vector<T>& out){
	os << "{";
	if (out.size() == 0) return os << " }";
	for (std::size_t i = 0; i < out.size()-1; ++i) {
		if (out[i] >= 0) os << " ";
		os << out[i] << ",";
	}
	if (out.back() >= 0) os << " ";
	return os << out.back() << " }";
}

// specialization of above template which "transposes" particle vectors
template<>
inline std::ostream& operator<<(std::ostream& os, const std::vector<particle>& out){
	if (out.empty()) return os << "{{ }{ }}";

	os << "{{";
	for (std::size_t i = 0; i < out.size()-1; ++i) {
		if (out[i].pm >= 0) os << " ";
		os << std::to_string(out[i].pm) << ",";
	}
	if (out.back().pm >= 0) os << " ";
	os << out.back().pm << " }{";

	for (std::size_t i = 0; i < out.size()-1; ++i) {
		if (out[i].pt >= 0) os << " ";
		os << std::to_string(out[i].pt) << ",";
	}
	if (out.back().pt >= 0) os << " ";
	os << out.back().pt << " }}";
	return os;
}

// specialization for vectors of chars which displays them as numbers
template<>
inline std::ostream& operator<<(std::ostream& os, const std::vector<char>& out){
	os << "{";
	if (out.size() == 0) return os << " }";
	for (std::size_t i = 0; i < out.size()-1; ++i) {
		if (out[i] >= 0) os << " ";
		os << static_cast<int>(out[i]) << ",";
	}
	if (out.back() >= 0) os << " ";
	return os << out.back() << " }";
}

template<typename T, int N>
inline std::ostream& operator<<(std::ostream& os, const std::array<T,N>& out){
	os << "{";
	for (std::size_t i = 0; i < out.size()-1; ++i) {
		if (out[i] >= 0) os << " ";
		os << out[i] << ",";
	}
	if (out.back() >= 0) os << " ";
	return os << out.back() << " }";
}

/*
// note: triplets displayed (row, column, value) despite matrix's storage type
inline std::ostream& operator<<(std::ostream& os, const Triplet& out) {
	return os << "(" << out.row() << "," << out.col() << "," << out.value()
		<< ")";
}
*/

#endif
