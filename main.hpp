#ifndef MAIN_HPP
#define MAIN_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

constexpr char VERSION[] = "0.9.7";
constexpr char RELEASE_DATE[] = __DATE__;

#include "constants.hpp"
#include "gui/main_window.hpp"
#include "calculation.hpp"
#include "test.hpp"

// startup and input parsing --------------------------------------------------

Arguments ParseArguments(int argc, char* argv[]);
int ParseOptions(std::vector<std::string> options);

// templates for reading inputs -----------------------------------------------

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
