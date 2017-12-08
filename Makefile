#-------------------------------------------------------------------------------
# Central Makefile for 3dBasis
#-------------------------------------------------------------------------------

BASEDIR = $(CURDIR)

CXX = clang++

CXXFLAGS = -IEigen -Wall -Wextra -pedantic -O3 -g -c -std=c++14
LDFLAGS = -lgsl -lblas -lpthread

EXECUTABLE = 3dBasis

SOURCES = mono.cpp poly.cpp 3dBasis.cpp multinomial.cpp matrix.cpp \
		  gram-schmidt.cpp
OBJECTS = mono.o poly.o 3dBasis.o multinomial.o matrix.o gram-schmidt.o

default: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $+ $(LDFLAGS) -o $@

3dBasis.o: 3dBasis.cpp 3dBasis.hpp mono.hpp poly.hpp basis.hpp \
	timer.hpp matrix.hpp multinomial.hpp gram-schmidt.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

mono.o: mono.cpp mono.hpp io.hpp constants.hpp construction.hpp 
	$(CXX) $(CXXFLAGS) $< -o $@

poly.o: poly.cpp poly.hpp mono.hpp io.hpp constants.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

multinomial.o: multinomial.cpp multinomial.hpp constants.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

gram-schmidt.o: gram-schmidt.cpp constants.hpp timer.hpp basis.hpp mono.hpp \
	poly.hpp matrix.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

matrix.o: matrix.cpp matrix.hpp multinomial.hpp mono.hpp basis.hpp io.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f *.o && rm -f $(EXECUTABLE)
