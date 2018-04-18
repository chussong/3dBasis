#-------------------------------------------------------------------------------
# Central Makefile for 3dBasis
#-------------------------------------------------------------------------------

BASEDIR := $(CURDIR)

# if you're not using clang, change this (e.g. to g++)
CXX := clang++

# depending on your configuration of Qt, you may have to remove -fPIC. It should
# tell you something about position-independent code if you need to do this
CXXFLAGS := -IEigen -Wall -Wextra -pedantic -fPIC -O3 -g -c -I$(BASEDIR)

# these should be the locations of your qt installation
QTINC := /usr/include/x86_64-linux-gnu/qt5
QTLIB := /usr/lib/x86_64-linux-gnu/

CXXFLAGS_CORE := -$(CXXFLAGS) -std=c++14
CXXFLAGS_QT := $(CXXFLAGS) -I$(QTINC) -std=c++17

LDFLAGS := -lgsl -lblas -lpthread 

LDFLAGS_CORE := $(LDFLAGS)
LDFLAGS_QT := -L$(QTLIB) -lQt5Widgets -lQt5Gui -lQt5Core $(LDFLAGS)

EXECUTABLE := 3dBasis

SOURCES_CORE := main.cpp calculation.cpp mono.cpp poly.cpp multinomial.cpp \
		matrix.cpp gram-schmidt.cpp discretization.cpp test.cpp \
SOURCES_QT := gui/main_window.cpp gui/moc_main_window.cpp gui/calc_widget.cpp \
	  gui/moc_calc_widget.cpp gui/file_widget.cpp gui/moc_file_widget.cpp \
	  gui/console_widget.cpp gui/moc_console_widget.cpp

#RESOURCES := gui/resources.rcc

OBJECTS_CORE := $(SOURCES_CORE:.cpp:=.o)
OBJECTS_QT := $(SOURCES_QT:.cpp:=.o)

default: $(EXECUTABLE)_GUI

nogui: $(EXECUTABLE)_NOGUI

$(EXECUTABLE)_GUI: $(OBJECTS_CORE) $(OBJECTS_QT)
	$(CXX) $(OBJECTS_CORE) $(OBJECTS_QT) $(LDFLAGS_QT) -o $(EXECUTABLE)

$(EXECUTABLE)_NOGUI: $(OBJECTS_CORE)
	$(CXX) $(OBJECTS_CORE) $(LDFLAGS) -o $(EXECUTABLE)

main.o: main.cpp main.hpp constants.hpp calculation.hpp gui/main_window.hpp
	$(CXX) $(CXXFLAGS_QT) $< -o $@

calculation.o: calculation.cpp calculation.hpp constants.hpp construction.hpp \
	mono.hpp poly.hpp basis.hpp io.hpp timer.hpp gram-schmidt.hpp \
	matrix.hpp multinomial.hpp discretization.hpp test.hpp
	$(CXX) $(CXXFLAGS_CORE) $< -o $@

mono.o: mono.cpp mono.hpp io.hpp constants.hpp construction.hpp 
	$(CXX) $(CXXFLAGS_CORE) $< -o $@

poly.o: poly.cpp poly.hpp mono.hpp io.hpp constants.hpp
	$(CXX) $(CXXFLAGS_CORE) $< -o $@

multinomial.o: multinomial.cpp multinomial.hpp constants.hpp io.hpp
	$(CXX) $(CXXFLAGS_CORE) $< -o $@

gram-schmidt.o: gram-schmidt.cpp constants.hpp timer.hpp basis.hpp mono.hpp \
	poly.hpp matrix.hpp
	$(CXX) $(CXXFLAGS_CORE) $< -o $@

matrix.o: matrix.cpp matrix.hpp multinomial.hpp mono.hpp basis.hpp io.hpp \
    	discretization.hpp constants.hpp
	$(CXX) $(CXXFLAGS_CORE) $< -o $@

discretization.o: discretization.cpp discretization.hpp constants.hpp mono.hpp \
	basis.hpp hypergeo.hpp
	$(CXX) $(CXXFLAGS_CORE) $< -o $@

test.o: test.cpp test.hpp io.hpp discretization.hpp matrix.hpp gram-schmidt.hpp\
    	hypergeo.hpp constants.hpp
	$(CXX) $(CXXFLAGS_CORE) $< -o $@

gui/main_window.o: gui/main_window.cpp gui/main_window.hpp constants.hpp \
    	gui/calc_widget.hpp gui/file_widget.hpp gui/console_widget.hpp
	$(CXX) $(CXXFLAGS_QT) $< -o $@

gui/moc_main_window.o: gui/moc_main_window.cpp
	$(CXX) $(CXXFLAGS_QT) $< -o $@

gui/moc_main_window.cpp: gui/main_window.hpp
	moc $< -o $@

gui/calc_widget.o: gui/calc_widget.cpp gui/calc_widget.hpp constants.hpp \
    	calculation.hpp
	$(CXX) $(CXXFLAGS_QT) $< -o $@

gui/moc_calc_widget.o: gui/moc_calc_widget.cpp
	$(CXX) $(CXXFLAGS_QT) $< -o $@

gui/moc_calc_widget.cpp: gui/calc_widget.hpp
	moc $< -o $@

gui/file_widget.o: gui/file_widget.cpp gui/file_widget.hpp
	$(CXX) $(CXXFLAGS_QT) $< -o $@

gui/moc_file_widget.o: gui/moc_file_widget.cpp
	$(CXX) $(CXXFLAGS_QT) $< -o $@

gui/moc_file_widget.cpp: gui/file_widget.hpp
	moc $< -o $@

gui/console_widget.o: gui/console_widget.cpp gui/console_widget.hpp
	$(CXX) $(CXXFLAGS_QT) $< -o $@

gui/moc_console_widget.o: gui/moc_console_widget.cpp
	$(CXX) $(CXXFLAGS_QT) $< -o $@

gui/moc_console_widget.cpp: gui/console_widget.hpp
	moc $< -o $@

#gui/resources.rcc: gui/resources.qrc
	#rcc $< -o $@

clean:
	rm -f $(OBJECTS_CORE) $(OBJECTS_QT) && rm -f $(EXECUTABLE)
