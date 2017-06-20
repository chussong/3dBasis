CXX = clang++

CXXFLAGS = -Wall -Wextra -pedantic -g -c -std=c++14
LDFLAGS = -lgmp
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
