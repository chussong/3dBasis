#-------------------------------------------------------------------------------
# Central Makefile for 3dBasis
#-------------------------------------------------------------------------------

BASEDIR = $(CURDIR)

CXX = clang++

CXXFLAGS = -IEigen -ISuiteSparse/include -Wall -Wextra -pedantic -O3 -g -c -std=c++14
LDFLAGS = -L$(BASEDIR)/lib -Wl,-rpath=$(BASEDIR)/lib,-rpath=$(BASEDIR)/SuiteSparse/lib \
		  -lgsl -lspqr -lcholmod -lmetis -lopenblas -lpthread
EXECUTABLE = 3dBasis

SOURCES = mono.cpp poly.cpp 3dBasis.cpp cache.cpp binomial.cpp multinomial.cpp \
		  matrix.cpp
OBJECTS = mono.o poly.o 3dBasis.o cache.o binomial.o multinomial.o matrix.o

default: $(EXECUTABLE)

all:  | lib
	( cd LAPACK && $(MAKE) lib )
	( ln -s $(BASEDIR)/LAPACK/liblapack.a $(BASEDIR)/lib/liblapack.a )
	( cd OpenBLAS && $(MAKE) libs )
	( ln -s $(BASEDIR)/OpenBLAS/libopenblas.a $(BASEDIR)/lib/libopenblas.a )
	( cd SuiteSparse && $(MAKE) library \
		LDFLAGS="-L$(BASEDIR)/SuiteSparse/lib -L$(BASEDIR)/lib" \
		BLAS="-lopenblas -lgfortran -lpthread" )
	( ln -s $(BASEDIR)/SuiteSparse/lib/libspqr.so $(BASEDIR)/lib/libspqr.so )
	( ln -s $(BASEDIR)/SuiteSparse/lib/libcholmod.so $(BASEDIR)/lib/libcholmod.so )
	( ln -s $(BASEDIR)/SuiteSparse/lib/libmetis.so $(BASEDIR)/lib/libmetis.so )
	( $(MAKE) $(EXECUTABLE) )

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $+ $(LDFLAGS) -o $@

3dBasis.o: 3dBasis.cpp 3dBasis.hpp mono.hpp poly.hpp basis.hpp lib/libspqr.so \
	timer.hpp cache.hpp matrix.hpp multinomial.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

mono.o: mono.cpp mono.hpp io.hpp constants.hpp construction.hpp cache.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

poly.o: poly.cpp poly.hpp mono.hpp io.hpp constants.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

cache.o: cache.cpp cache.hpp constants.hpp io.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

binomial.o: binomial.cpp binomial.hpp constants.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

multinomial.o: multinomial.cpp multinomial.hpp constants.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

matrix.o: matrix.cpp matrix.hpp binomial.hpp multinomial.hpp mono.hpp basis.hpp \
	io.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

lib/libspqr.so: lib/liblapack.a lib/libopenblas.a
	( cd SuiteSparse && $(MAKE) library \
		LDFLAGS="-L$(BASEDIR)/SuiteSparse/lib -L$(BASEDIR)/lib" \
		BLAS="-lopenblas -lgfortran -lpthread" )
	( ln -s $(BASEDIR)/SuiteSparse/lib/libspqr.so $(BASEDIR)/lib/libspqr.so )
	( ln -s $(BASEDIR)/SuiteSparse/lib/libcholmod.so $(BASEDIR)/lib/libcholmod.so )
	( ln -s $(BASEDIR)/SuiteSparse/lib/libmetis.so $(BASEDIR)/lib/libmetis.so )

lib/liblapack.a: | lib
	( cd LAPACK && $(MAKE) lib )
	( ln -s $(BASEDIR)/LAPACK/liblapack.a $(BASEDIR)/lib/liblapack.a )

lib/libopenblas.a: | lib
	( cd OpenBLAS && $(MAKE) libs )
	( ln -s $(BASEDIR)/OpenBLAS/libopenblas.a $(BASEDIR)/lib/libopenblas.a )

lib:
	mkdir -p lib

clean:
	rm -f *.o && rm -f $(EXECUTABLE)
