#ifndef IO_HPP
#define IO_HPP

#include <iostream>
#include <sstream> // for mathematica output formatting
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

inline std::string MathematicaOutput(const DMatrix& out) {
	std::stringstream stream;
	stream << "{";
	for (Eigen::Index row = 0; row < out.rows() - 1; ++row) {
		stream << "{";
		for (Eigen::Index col = 0; col < out.cols() - 1; ++col) {
			stream << out(row, col) << ", ";
		}
		stream << out(row, out.cols()-1) << "},\n";
	}

	stream << "{";
	for (Eigen::Index col = 0; col < out.cols() - 1; ++col) {
		stream << out(out.rows()-1, col) << ", ";
	}
	stream << out(out.rows()-1, out.cols()-1) << "}}";
	stream.flush();
	return stream.str();
}

#endif
