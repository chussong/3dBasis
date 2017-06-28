# 3dBasis
A program for numerical determination of conformal primary operators in three 
dimensions.  

Copyright 2017 Charles Hussong  

Project homepage:	https://github.com/chussong/3dBasis  
Contact email:		chussong@jhu.edu  

## License

3dBasis uses [Eigen](http://eigen.tuxfamily.org/) to perform its basic matrix 
operations, and [SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html)
to perform the QR decomposition at the core of the computation. When building
SuiteSparse for 3dBasis, we also use the [LAPACK](http://www.netlib.org/lapack/)
and [OpenBLAS](http://www.openblas.net/) low-level matrix manipulation libraries.  

## Installation

To compile, get a terminal in this directory and type "make". If you have C++
and Fortran compilers available everything will ideally be handled automatically.  

#Note: the dependencies currently do not build correctly inside a directory with
a space in its path. For now, make sure you don't have one (this means you can't
build it in the shared Dropbox folder).#  

Mac users can get a C++ compiler by running xcode-select --install to get the 
Command Line Tools. Unfortunately, they do not contain a Fortran compiler, so 
you will also need to get gfortran; I suggest downloading it from [the GCC wiki]
(http://gcc.gnu.org/wiki/GFortranBinariesMacOS). You can also use a different
Fortran compiler by changing "gfortran" to something else in the LAPACK Makefile.  

## Usage

For now, the program only does one thing: when invoked with ./3dBasis N D, where
N and D are nonnegative integers, it outputs the number of primary operators
with N particles and D derivatives. It also constructs an orthonormal basis for
the space of these primaries, but is currently not set to output it; if you want
this information, there is a commented-out line at the end of main() in 
3dBasis.cpp. The primary vectors will be different from those output by
Mathematica, but they should both span the same space.  

The current version of 3dBasis is configured to use double precision floating
point numbers. This can be changed in 3dBasis.hpp by changing the line with the
typedef for coeff\_class. Hopefully we will not need higher precision, but if
necessary it's possible to increase it arbitrarily. This would certainly cause a
severe performance degradation.  

## Options

There are several options available: 

| Option | Description |
| ------ | ----------- |
| -b | perform pure brute force calculation without using subspaces |
| -e | require that monomials obey the equations of motion (experimental) |
| -v | instead of running, print the version and date of the release and exit |

## Notes on Efficiency

The underlying libraries are configured to make use of multiple cores and to run
computations for the QR decomposition on the GPU, so this program will generally
perform best on a computer with several cores and a graphics card.  

The slow part of the computation is the QR decomposition; I was not able to find
a discussion of the decomposition's complexity, but it is likely order 
(size of the basis)^3^, so the bet way to improve efficiency is to split the 
computation into parts which have factorized bases. I already do this for 
odd/even perp parity, but there is also a symmetry under exchange of the minus
and plus components which allows a reduction in basis size. The target basis
can also be split using the distinct perp-parity properties of the different
boost generators. I have not yet implemented either of these.  
