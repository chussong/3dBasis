# 3dBasis
A program for numerical determination of conformal primary operators in three 
dimensions.  

Copyright 2017-2018 Charles Hussong  

Project homepage:	https://github.com/chussong/3dBasis  
Contact email:		chussong@jhu.edu  

## License

3dBasis uses [Eigen](http://eigen.tuxfamily.org/) to perform its basic matrix 
operations, which is included. You will also need a Basic Linear Algebra
System (BLAS), which is not included; many systems come with one pre-installed,
but if yours doesn't, try [OpenBLAS](http://www.openblas.net/).  

The work-in-progress GUI uses [Qt](https://www.qt.io), an LGPL cross-platform 
GUI library; this project will be released under the GPL to comply with their 
licensing requirements, but I haven't actually done this yet.  

## Installation

To compile, get a terminal in this directory and type "make". If you have Clang
available everything will ideally be handled automatically; if you'd rather use
another compiler, edit the top of the Makefile.  

Mac users can get a C++ compiler by running xcode-select --install to get the 
Command Line Tools. This will give you Clang, so you won't have to change CXX in
the Makefile.  

To build with the GUI (the default), you'll need to install 
[Qt](https://www.qt.io/). On Ubuntu, you can just do 
"sudo apt install qt5-default" and then "make" should just work; if you're not 
on Ubuntu, you may have to change QTINC and QTLIB at the top of the Makefile to 
the appropriate directories for your installation. QTINC is the "include" 
directory containing subdirectories QtCore, QtGui, and QtWidgets; QTLIB is the 
"lib" directory containing libQt5Widgets.so, libQt5Gui.so, and libQt5Core.so.  

If you don't want to bother with Qt, you can also build without the GUI; do this
with "make nogui" instead of the usual "make". This build should successfully
complete without Qt installed, and makes exactly the same calculations but
parameters must be entered from the command line.  

Note that if you switch between the default and nogui versions, you will likely
have to rebuild everything. Use "make clean" to delete all the existing compiled
files, after which you can "make" or "make nogui" as usual.  

Finally, it's possible to make the GUI version with static linking so that it 
will run on a system without Qt; doing this requires building Qt with the static
versions of its libraries, which is not the default, so I haven't really tested
it yet. The command to do this is "make static".  

## Usage

There are two ways to invoke the program: if called without any arguments, as a
simple "./3dBasis", it will open up a GUI that controls all of the maintained
options. Of particular importance are "n", the number of particles, "l", the
number of non-Dirichlet derivatives, and "p", the number of mu partitions.  

These arguments can also be supplied in the initial invocation, in which case
the GUI will not be launched at all. The syntax for this is "./3dBasis n l p",
and additional options can be placed anywhere preceded by a minus sign "-". The
full list of options is given below, but the most important are "-o" and "-O",
which cause the results of the computation to be output as Mathematica code 
instead of written to the terminal. "-o file.txt" will cause output to be 
written to file.txt, appending to the file if it exists already, and 
"-O file.txt" will do the same but overwrite the file if it exists already.  

The actual output of the computation is a basis of orthogonal states found via
the Gram-Schmidt process, followed by a succession of matrices representing 
terms of the finite-dimensional Hamiltonian on this basis. These are given first
as matrices between the "minimal basis" of monomials used to express the basis
states, and then as matrices between the basis states themselves.  

Gram-Schmidt is done using a custom algorithm which makes it easy to compute
matrix elements but it's unfortunately not very stable. I'm looking into doing
this with Householder reflections, but the basis obtained from this is somewhat
lower quality when rounding errors are not an issue. Details of this are in
gram-schmidt.cpp.  

The current version of 3dBasis is configured to use quadruple precision floating
point numbers, with fallback to double precision for many math routines. These 
can be changed in constants.hpp by changing the lines with the typedefs for 
coeff\_class and builtin\_class; coeff\_class is the actual type of the 
coefficients used all over the program, while builtin\_class is a type to which
coeff\_class is convertible and for which the \<cmath\> functions are defined. 
If you need the higher precision in all cases, you can define both of these 
types to be the same and add overloads of the \<cmath\> functions for them. This
would certainly cause a severe performance degradation, so it would probably be
best to try to figure out where exactly you need the greater precision and only
introduce it there.  

## Options

There are several options available: 

| Option | Description |
| ------ | ----------- |
| -d | debug mode, producing some extra output (currently always on) |
| -i | include interaction terms in the Hamiltonian (the default is a free theory) |
| -m | perform a test of the multinomial module, then exit |
| -M | use only all-minus states with no transverse momentum |
| -o \<filename\> | write non-error output to \<filename\> instead of the terminal. If this file exists, it will be APPENDED TO |
| -O \<filename\> | write non-error output to \<filename\> instead of the terminal. If this file exists, it will be OVERWRITTEN |
| -s | only orthogonal basis states, output them, then exit |
| -t | perform all automated unit tests, then exit |
| -v | instead of running, print the version and date of release, then exit |
