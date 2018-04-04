#ifndef MAIN_HPP
#define MAIN_HPP

#include <iostream>
#include <fstream>
#include <exception>
#include <string>
#include <vector>
#include <array>
#include <list>
#include <utility>		// std::pair
#include <algorithm>	// std::remove_if
#include <type_traits>	// std::is_same
#include <gsl/gsl_errno.h>  // handling for GSL errors

constexpr char VERSION[] = "0.9.5";
constexpr char RELEASE_DATE[] = __DATE__;

#include "constants.hpp"
#include "construction.hpp"
#include "mono.hpp"
#include "poly.hpp"
#include "basis.hpp"
#include "io.hpp"
#include "timer.hpp"
#include "gram-schmidt.hpp"
#include "matrix.hpp"
#include "multinomial.hpp" // the coefficients are initialized in main()
#include "discretization.hpp"
#include "test.hpp"

// startup and input parsing --------------------------------------------------

Arguments ParseArguments(int argc, char* argv[]);
int ParseOptions(std::vector<std::string> options);
void GSLErrorHandler(const char* reason, const char* file, int line, int err);

// actual computations --------------------------------------------------------

std::vector<Poly> ComputeBasisStates(const Arguments& args);
std::vector<Poly> ComputeBasisStates_SameParity(
        const std::vector<Basis<Mono>>& inputBases, const Arguments& args);
DMatrix PolysOnMinBasis(const Basis<Mono>& minimalBasis,
        const std::vector<Poly> orthogonalized, std::ostream& outStream);
DMatrix ComputeHamiltonian(const Arguments& args);
DMatrix ComputeHamiltonian_SameParity(const std::vector<Basis<Mono>>& inputBases,
                                      const Arguments& args);

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
