#ifndef MAIN_HPP
#define MAIN_HPP

#include <iostream>
#include <exception>
#include <string>
#include <vector>
#include <array>
#include <list>
#include <utility>		// std::pair
#include <algorithm>	// std::remove_if

constexpr char VERSION[] = "0.4.3";
constexpr char RELEASE_DATE[] = __DATE__;

#include "constants.hpp"
#include "construction.hpp"
#include "mono.hpp"
#include "poly.hpp"
#include "basis.hpp"
#include "io.hpp"

std::ostream& operator<<(std::ostream& os, const Triplet& out);

// startup and input parsing --------------------------------------------------

arguments ParseArguments(int argc, char* argv[]);
int ParseOptions(std::vector<std::string> options);

int FindPrimaries(const arguments& args);
int FindPrimariesParityOnly(const arguments& args);
int FindPrimariesBruteForce(const arguments& args);
int FindPrimariesByM(const arguments& args);
void AddPrimariesAtL(const mBasis& startBasis, const mBasis& targetBasis,
		const unsigned int L, std::vector<poly>& primaries, 
		const coeff_class delta);

// functions interfacing with Eigen ------------------------------------------

Matrix KMatrix(const Basis<mono>& startingBasis, const Basis<mono>& targetBasis,
		const coeff_class delta);
Matrix K13Matrix(const Basis<mono>& startingBasis, const Basis<mono>& targetBasis,
		const coeff_class delta);
Matrix K2Matrix(const Basis<mono>& startingBasis, const Basis<mono>& targetBasis,
		const coeff_class delta);
std::array<Matrix,4> KMatrices(const splitBasis<mono>& startingBasis,
		const splitBasis<mono>& targetBasis, const coeff_class delta);
std::list<Triplet> ConvertToRows(const std::vector<poly>& polyForms, 
		const Basis<mono>& targetBasis, const Eigen::Index rowOffset);
std::vector<poly> Kernel(const Matrix& KActions, const Basis<mono>& startBasis,
		const bool outputKernel);
std::vector<poly> CombineKernels(const std::vector<poly>& kernel1,
		const std::vector<poly>& kernel2);
poly VectorToPoly(const Vector& kernelVector, const Basis<mono>& startBasis);
poly VectorToPoly(const DVector& kernelVector, const Basis<mono>& startBasis);
poly ColumnToPoly(const Matrix& kernelMatrix, const Eigen::Index col, 
		const Basis<mono>& startBasis);
poly ColumnToPoly(const DMatrix& kernelMatrix, const Eigen::Index col, 
		const Basis<mono>& startBasis);


// templates -----------------------------------------------------------------

template<typename ParseTo>
ParseTo ReadArg(const std::string& arg){
	std::cerr << "Error: attempted to parse the argument " << arg << " to a "
		<< "type with no known parsing function. Please specialize the ReadArg "
		<< "template to your coeff_class." << std::endl;
	return ParseTo();
}

template<>
int ReadArg<int>(const std::string& arg){
	int ret;
	try{ret = std::stoi(arg);}
	catch(const std::invalid_argument &e){
		std::cerr << "Error: this non-option argument could not be "
			<< "converted to an integer: " << arg << std::endl;
		throw;
	}
	catch(const std::out_of_range &e){
		std::cerr << "Error: specification of N or degree is too "
			<< "large to store. This computation would never finish"
			<< " anyway..." << std::endl;
		throw;
	}
	return ret;
}

template<>
double ReadArg<double>(const std::string& arg){
	double ret;
	try{ret = std::stod(arg);}
	catch(const std::invalid_argument &e){
		std::cerr << "Error: this non-option argument could not be "
			<< "converted to an integer: " << arg << std::endl;
		throw;
	}
	catch(const std::out_of_range &e){
		std::cerr << "Error: specification of Delta is too large to store."
			<< std::endl;
		throw;
	}
	return ret;
}

#endif
