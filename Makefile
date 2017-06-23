CXX = clang++

CXXFLAGS = -Wall -Wextra -pedantic -O3 -g -c -std=c++14
LDFLAGS = -lmpfr -lgmp -lspqr -lcholmod
EXECUTABLE = 3dBasis

SOURCES = 3dBasis.cpp
OBJECTS = 3dBasis.o

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $+ $(LDFLAGS) -o $@

3dBasis.o: 3dBasis.cpp 3dBasis.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f *.o && rm -f $(EXECUTABLE)

# eventually we want make to:
# 0? compile LAPACK or BLAS for local system
# 1: compile suitesparse for local system (linking to LAPACK/BLAS if applicable)
# 2: compile GMP for local system
# 3: compile mpfr for local system
# 4: compile 3dBasis
# So all of the above should be included in subdirectories.
