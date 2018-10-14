Copyright German Sinuco 2018
Distributed under the Boost Software License, Version 1.0.
(See a copy at http://www.boost.org/LICENSE_1_0.txt)


This library implements the multimode floquet approach to calculate the
time-evolution operator of a time-dependent quantum system with discrete 
spectrum, where the time-dependence is given by a superposition of harmonic
coupligs.

HOW TO INSTALL

A simple Makefile is included in the distribution. As a requirement,
your system should have LAPACK. The library also builds a sparse representation
of the Hamiltonian, which can be used with sparse dedicated software such as the intel MKL.
Subroutines using functions of the intel-MKL are also implements, and to use them you would 
also need to have this library install in your systems.

With this, to install the software you mus specify the path to the required libraries and 
conrresponding include folders (LAPACK and MKL (optional)).

The full library is build with the command:

make lib  

the LAPACK dependent componets with

make lib_lapack

Both commands build all needed object files, put toghether the file
lib/libmultimodefloquet.a and copy the fortran modules in the directory 
include/.


Examples are in the folders examples/FORTRAN and examples/CPP. In each 
directory, a Makefile takes care of creating executables.  



