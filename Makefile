#-------------------------------------------------------------------------------
# Central Makefile for 3dBasis
#-------------------------------------------------------------------------------

BASEDIR = $(CURDIR)

# if you're not using clang, change this (e.g. to g++)
CXX = clang++

# these should be the locations of your qt installation
QTINC = /usr/include/x86_64-linux-gnu/qt5
QTLIB = /usr/lib/x86_64-linux-gnu/
INC = -I$(QTINC) -I$(BASEDIR)

# depending on your configuration of Qt, you may have to remove -fPIC. It should
# tell you something about position-independent code if you need to do this
CXXFLAGS = -IEigen -Wall -Wextra -pedantic -Wno-c++1z-extensions -fPIC -O3 -g \
	  -c -std=c++14 $(INC)
LDFLAGS = -lgsl -lblas -L$(QTLIB) -lQt5Widgets -lQt5Gui -lQt5Core -lpthread 

EXECUTABLE = 3dBasis

SOURCES = main.cpp calculation.cpp mono.cpp poly.cpp multinomial.cpp \
	  matrix.cpp gram-schmidt.cpp discretization.cpp test.cpp \
	  gui/main_window.cpp gui/moc_main_window.cpp resources.qrc

OBJECTS = $(SOURCES:.cpp=.o) resources.rcc

default: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $+ $(LDFLAGS) -o $@

main.o: main.cpp main.hpp constants.hpp calculation.hpp test.hpp \
	gui/main_window.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

calculation.o: calculation.cpp calculation.hpp constants.hpp construction.hpp \
	mono.hpp poly.hpp basis.hpp io.hpp timer.hpp gram-schmidt.hpp \
	matrix.hpp multinomial.hpp discretization.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

mono.o: mono.cpp mono.hpp io.hpp constants.hpp construction.hpp 
	$(CXX) $(CXXFLAGS) $< -o $@

poly.o: poly.cpp poly.hpp mono.hpp io.hpp constants.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

multinomial.o: multinomial.cpp multinomial.hpp constants.hpp io.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

gram-schmidt.o: gram-schmidt.cpp constants.hpp timer.hpp basis.hpp mono.hpp \
	poly.hpp matrix.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

matrix.o: matrix.cpp matrix.hpp multinomial.hpp mono.hpp basis.hpp io.hpp \
    	discretization.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

discretization.o: discretization.cpp discretization.hpp constants.hpp mono.hpp \
	basis.hpp hypergeo.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

test.o: test.cpp test.hpp io.hpp discretization.hpp matrix.hpp gram-schmidt.hpp\
    	hypergeo.hpp constants.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

gui/main_window.o: gui/main_window.cpp gui/main_window.hpp constants.hpp \
    	calculation.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

gui/moc_main_window.o: gui/moc_main_window.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

gui/moc_main_window.cpp: gui/main_window.hpp
	moc $< -o $@

gui/resources.rcc: gui/resources.qrc
	rcc $< -o $@

clean:
	rm -f *.o && rm -f $(EXECUTABLE)
