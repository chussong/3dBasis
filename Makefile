#-------------------------------------------------------------------------------
# Central Makefile for 3dBasis
#-------------------------------------------------------------------------------

.PHONY: default nogui static clean

# these should be the locations of your Qt installation and Qt's "moc" tool; if
# necessary, change the one for your OS, or just add yours at the end
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    QTINC := -I/usr/include/x86_64-linux-gnu/qt5
    QTLIB := -L/usr/lib/x86_64-linux-gnu -lQt5Widgets -lQt5Gui -lQt5Core
    MOC := moc
endif
ifeq ($(UNAME_S),Darwin)
    QTINC := -I/usr/local/opt/qt/include
    QTLIB := -F/usr/local/opt/qt/lib -framework QtWidgets -framework QtGui \
	-framework QtCore
    MOC := /usr/local/opt/qt/bin/moc
endif
#FIXME: SEND PREPROCESSOR USING_GCC FLAG


ifndef SECOND_PASS
#-------------------------------------------------------------------------------
# set up environment on first pass
#-------------------------------------------------------------------------------

export SECOND_PASS := true

default static: export CXXFLAGS_EXTRA := $(QTINC)
default static: export GUI_FROM_MAIN := gui/main_window.hpp
default:
	@$(MAKE) default
static:
	@$(MAKE) static
nogui: export CXXFLAGS_EXTRA := -DNO_GUI
nogui: export GUI_FROM_MAIN :=
nogui:
	@$(MAKE) nogui

clean:
	@$(MAKE) clean

# FIXME: GNU requires 'install' targets as well

else
#-------------------------------------------------------------------------------
# do the build on the second pass
#-------------------------------------------------------------------------------

BASEDIR := $(CURDIR)

# default values for GNU standard arguments the user could override
# if you're not using clang, change CXX (e.g. to g++)
CXX := clang++
CFLAGS :=
CXXFLAGS := $(CFLAGS) -O3 -g
LDFLAGS :=

# depending on your configuration of Qt, you may have to remove -fPIC. It should
# tell you something about position-independent code if you need to do this
CXXFLAGS_GLOBAL := -IEigen -Wall -Wextra -pedantic -fPIC -c -I$(BASEDIR) \
    		   -std=c++14 -Wno-c++1z-extensions $(CXXFLAGS)

CXXFLAGS_CORE := $(CXXFLAGS_GLOBAL) $(CXXFLAGS_EXTRA)
CXXFLAGS_QT := $(CXXFLAGS_GLOBAL) $(QTINC)

LDFLAGS_GLOBAL := -lgsl -lblas -lpthread -lboost_filesystem -lboost_system $(LDFLAGS)

LDFLAGS_CORE := $(LDFLAGS_GLOBAL)
LDFLAGS_QT := $(QTLIB) $(LDFLAGS_GLOBAL)

EXECUTABLE := 3dBasis

SOURCES_CORE := main.cpp calculation.cpp mono.cpp poly.cpp multinomial.cpp \
		matrix.cpp gram-schmidt.cpp discretization.cpp test.cpp \
		hypergeo.cpp
SOURCES_QT := gui/main_window.cpp gui/moc_main_window.cpp gui/calc_widget.cpp \
	  gui/moc_calc_widget.cpp gui/file_widget.cpp gui/moc_file_widget.cpp \
	  gui/console_widget.cpp gui/moc_console_widget.cpp

#RESOURCES := gui/resources.rcc

OBJECTS_CORE := $(SOURCES_CORE:.cpp=.o)
OBJECTS_QT := $(SOURCES_QT:.cpp=.o)
OBJECTS_ALL := $(OBJECTS_QT) $(OBJECTS_CORE)

#-------------------------------------------------------------------------------
# meta targets
#-------------------------------------------------------------------------------

default: $(EXECUTABLE)_GUI

nogui: $(EXECUTABLE)_NOGUI

$(EXECUTABLE)_GUI: $(OBJECTS_ALL)
	$(CXX) $(OBJECTS_ALL) $(LDFLAGS_QT) -o $(EXECUTABLE)

static: $(OBJECTS_ALL)
	$(CXX) -static $(OBJECTS_ALL) $(LDFLAGS_QT) -o $(EXECUTABLE)

$(EXECUTABLE)_NOGUI: $(OBJECTS_CORE)
	$(CXX) $(OBJECTS_CORE) $(LDFLAGS_CORE) -o $(EXECUTABLE)

clean:
	rm -f $(OBJECTS_ALL) $(EXECUTABLE)

#-------------------------------------------------------------------------------
# core targets
#-------------------------------------------------------------------------------

main.o: main.cpp main.hpp constants.hpp calculation.hpp $(GUI_FROM_MAIN)
	$(CXX) $(CXXFLAGS_CORE) $< -o $@

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

discretization.o: discretization.cpp discretization.hpp constants.hpp \
	hypergeo.hpp multinomial.hpp
	$(CXX) $(CXXFLAGS_CORE) $< -o $@

hypergeo.o: hypergeo.cpp hypergeo.hpp constants.hpp
	$(CXX) $(CXXFLAGS_CORE) $< -o $@

test.o: test.cpp test.hpp io.hpp discretization.hpp matrix.hpp gram-schmidt.hpp\
    	hypergeo.hpp constants.hpp
	$(CXX) $(CXXFLAGS_CORE) $< -o $@

#-------------------------------------------------------------------------------
# gui targets
#-------------------------------------------------------------------------------

gui/main_window.o: gui/main_window.cpp gui/main_window.hpp constants.hpp \
    	gui/calc_widget.hpp gui/file_widget.hpp gui/console_widget.hpp
	$(CXX) $(CXXFLAGS_QT) $< -o $@

gui/moc_main_window.o: gui/moc_main_window.cpp
	$(CXX) $(CXXFLAGS_QT) $< -o $@

gui/moc_main_window.cpp: gui/main_window.hpp
	$(MOC) $< -o $@

gui/calc_widget.o: gui/calc_widget.cpp gui/calc_widget.hpp constants.hpp \
    	calculation.hpp
	$(CXX) $(CXXFLAGS_QT) $< -o $@

gui/moc_calc_widget.o: gui/moc_calc_widget.cpp
	$(CXX) $(CXXFLAGS_QT) $< -o $@

gui/moc_calc_widget.cpp: gui/calc_widget.hpp
	$(MOC) $< -o $@

gui/file_widget.o: gui/file_widget.cpp gui/file_widget.hpp
	$(CXX) $(CXXFLAGS_QT) $< -o $@

gui/moc_file_widget.o: gui/moc_file_widget.cpp
	$(CXX) $(CXXFLAGS_QT) $< -o $@

gui/moc_file_widget.cpp: gui/file_widget.hpp
	$(MOC) $< -o $@

gui/console_widget.o: gui/console_widget.cpp gui/console_widget.hpp
	$(CXX) $(CXXFLAGS_QT) $< -o $@

gui/moc_console_widget.o: gui/moc_console_widget.cpp
	$(CXX) $(CXXFLAGS_QT) $< -o $@

gui/moc_console_widget.cpp: gui/console_widget.hpp
	$(MOC) $< -o $@

#gui/resources.rcc: gui/resources.qrc
#	rcc $< -o $@
endif #SECOND_PASS
