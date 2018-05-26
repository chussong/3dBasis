#ifndef IO_HPP
#define IO_HPP

#include <iostream>
#include <sstream> // for mathematica output formatting
#include <iomanip>
#include <limits>
#ifndef NO_GUI
#include <QtCore/QTextStream>
#include <QtCore/QString>
#endif
#include "constants.hpp"

// stream output for particles
inline std::ostream& operator<<(std::ostream& os, const particle& out) {
    return os << "{ " << static_cast<int>(out.pm) << ", "
        << static_cast<int>(out.pt) << " }";
}

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

// specialization of above template which does not "transpose" particle vectors
template<>
inline std::ostream& operator<<(std::ostream& os, 
    const std::vector<particle>& out) {
    os << "{";
    if (out.size() == 0) return os << " }";
    for (std::size_t i = 0; i < out.size()-1; ++i) {
        os << out[i] << ", ";
    }
    return os << out.back() << "}";
}

// specialization of above template which does "transpose" particle vectors
/*template<>
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
}*/

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
    return os << static_cast<int>(out.back()) << " }";
}

// stream output for arrays
template<typename T, std::size_t N>
inline std::ostream& operator<<(std::ostream& os, const std::array<T,N>& out) {
    os << "{";
    for (std::size_t i = 0; i < out.size()-1; ++i) {
        if (out[i] >= 0) os << " ";
        os << out[i] << ",";
    }
    if (out.back() >= 0) os << " ";
    return os << out.back() << " }";
}

// partially specialized array output for characters
template<std::size_t N>
inline std::ostream& operator<<(std::ostream& os, const std::array<char,N>& out) {
    os << "{";
    for (std::size_t i = 0; i < out.size()-1; ++i) {
        if (out[i] >= 0) os << " ";
        os << static_cast<int>(out[i]) << ",";
    }
    if (out.back() >= 0) os << " ";
    return os << static_cast<int>(out.back()) << " }";
}

// note: triplets displayed (row, column, value) despite matrix's storage type
inline std::ostream& operator<<(std::ostream& os, const Triplet& out) {
	return os << "(" << out.row() << "," << out.col() << "," << out.value()
		<< ")";
}

namespace {
    constexpr char hexMap[16] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', 
        '9', 'A', 'B', 'C', 'D', 'E', 'F'};
} // anonymous namespace

// reinterpret a string as a sequence of 8-bit numbers
inline std::string MVectorOut(std::string mVector) {
	for (auto i = 0u; i < mVector.size(); ++i) {
		if (mVector[i] < 0 || mVector[i] > 15) {
			throw(std::out_of_range("MVectorOut range error"));
		}
		mVector[i] = hexMap[static_cast<int>(mVector[i])];
	}
	return mVector;
}

// this should output an exact decimal form of "out"; in particular, it should
// not use the "e" notation
inline std::string MathematicaOutput(const coeff_class out) {
    std::stringstream ss;
    ss.precision(std::numeric_limits<builtin_class>::max_digits10);
    ss << out;
    std::string stringForm = ss.str();

    std::size_t e = stringForm.find('e');
    if (e != std::string::npos) {
        stringForm.replace(e, 1, "*10^(");
        stringForm.append(")");
    }

    return stringForm;
}

inline std::string MathematicaOutput(const DMatrix& out) {
    if (out.rows() == 0 || out.cols() == 0) return "{ }";
    std::stringstream stream;
    stream << "{";
    for (Eigen::Index row = 0; row < out.rows() - 1; ++row) {
        stream << "{";
        for (Eigen::Index col = 0; col < out.cols() - 1; ++col) {
            stream << MathematicaOutput(out(row, col)) << ", ";
        }
        stream << MathematicaOutput(out(row, out.cols()-1)) << "},\n";
    }

    stream << "{";
    for (Eigen::Index col = 0; col < out.cols() - 1; ++col) {
        stream << MathematicaOutput(out(out.rows()-1, col)) << ", ";
    }
    stream << MathematicaOutput(out(out.rows()-1, out.cols()-1)) << "}}";
    stream.flush();
    return stream.str();
}

inline std::string MathematicaOutput(const SMatrix& out) {
    std::stringstream stream;
    stream << "SparseArray[{";
    for (Eigen::Index i = 0; i < out.outerSize(); ++i) {
        for (SMatrix::InnerIterator it(out, i); it; ++it) {
            stream << '{' << it.row()+1 << ',' << it.col()+1 << "} -> " 
                << it.value() << ", ";
        }
    }
    stream.flush();
    std::string output = stream.str();
    output.replace(output.size()-2, 2, "}]");
    return output;
}

#ifndef NO_GUI
template<typename T>
inline QTextStream& operator<<(QTextStream& os, T out) {
    std::stringstream ss;
    ss << out;
    // return os << QString(ss.str().c_str());
    // std::string str(ss.str());
    // return os << str.c_str();
    return os << ss.str().c_str();
}
// 
// template<>
// inline QTextStream& operator<<(QTextStream& os, const std::string& out) {
    // return os << QString(out.c_str());
// }
#endif

#endif
