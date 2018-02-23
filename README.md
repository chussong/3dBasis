# 3dBasis
A program for numerical determination of conformal primary operators in three 
dimensions.  

Copyright 2017 Charles Hussong  

Project homepage:	https://github.com/chussong/3dBasis  
Contact email:		chussong@jhu.edu  

## License

3dBasis uses [Eigen](http://eigen.tuxfamily.org/) to perform its basic matrix 
operations, which is included. You will also need a Basic Linear Algebra
System (BLAS), which is not included; many systems come with one pre-installed,
but if yours doesn't, try [OpenBLAS](http://www.openblas.net/).

## Installation

To compile, get a terminal in this directory and type "make". If you have clang
available everything will ideally be handled automatically; if you'd rather use
another compiler, edit the top of the Makefile.  

Mac users can get a C++ compiler by running xcode-select --install to get the 
Command Line Tools.  

## Usage

For now, the program only does one thing: when invoked with ./3dBasis N D, where
N and D are nonnegative integers, it performs a Gram-Schmidt reduction on all
Dirichlet operators with N particles and D derivatives to determine an 
orthogonal basis.  

Gram-Schmidt is done using a custom algorithm which makes it easy to compute
matrix elements but it's unfortunately not very stable. I'm looking into doing
this with Householder reflections, but the basis obtained from this is somewhat
lower quality when rounding errors are not an issue. Details of this are in
gram-schmidt.cpp.

The current version of 3dBasis is configured to use double precision floating
point numbers. This can be changed in 3dBasis.hpp by changing the line with the
typedef for coeff\_class. Hopefully we will not need higher precision, but if
necessary it should be possible to increase it arbitrarily. This would certainly
cause a severe performance degradation.  

## Options

There are several options available: 

| Option | Description |
| ------ | ----------- |
| -d | debug mode, producing lots of extra output (currently always on) |
| -m | perform a test of the multinomial module, then exit |
| -M | use only all-minus states with no transverse momentum |
| -o \<filename\> | write non-error output to \<filename\> instead of the terminal. If this file exists, it will be APPENDED TO |
| -O \<filename\> | write non-error output to \<filename\> instead of the terminal. If this file exists, it will be OVERWRITTEN |
| -t | perform all automated unit tests, then exit |
| -v | instead of running, print the version and date of release, then exit |
