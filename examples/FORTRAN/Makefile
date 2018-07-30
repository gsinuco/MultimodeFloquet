#export LD_LIBRARY_PATH="/opt/intel/compilers_and_libraries_2017/linux/mkl/lib/intel64"
# SET FORTRAN AND CPP COMPILERS
CPP = g++
CC  = gcc
GF  = gfortran
AR  = ar 
RANLIB = ranlib

# SET REQUIRED FLAGS
GFFLAGS    =  -lgsl -lgslcblas  -llapack -lblas -g
GFFLAGS_SP =  -m64  -w -fno-second-underscore -x f77-cpp-input  -lpthread -lm -ldl -lgsl -lgslcblas -O3 -llapack -lblas -g
MKLFLAGS   =  -lmkl_gf_lp64 -lmkl_sequential -lmkl_core
CFLAGS     =  -lgsl -lgslcblas   -static

#SET MKL-intel LIBRARY PATH
MKLLIBS = /opt/intel/compilers_and_libraries/linux/mkl/lib/intel64
#SET MKL-intel INCLUDE PATH
MKLINC = /opt/intel/compilers_and_libraries/linux/mkl/include	

#SET MULTIMODEFLOQUET LIBRARY PATH
MMFLIB = ../../lib/
#SET MULTIMODEFLOQUET INCLUDE PATH
MMFINC = ../../include/

###################################
# MAKE CPP EXAMPLES
###################################

all: Example_lib_f Example_lib_f_sp



Example_lib_f: main_qubit.f90 main_DressedQubit.f90
	$(GF) -o qubit         main_qubit.f90        -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) 
	$(GF) -o qubit_driver  main_qubit_DRIVER.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) 
	$(GF) -o spin_driver   main_spin_DRIVER.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) 
	$(GF) -o dressedqubit  main_DressedQubit.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) 
	$(GF) -o dressedqubit_DRIVER  main_DressedQubit_DRIVER.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) 
	$(GF) -o dressedRb87_DRIVER   main_87Rb_RFMW_DRIVER.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) 

Example_lib_f_sp: main_qubit_SP.f90 main_DressedQubit_SP.f90
	$(GF) -o  qubit_sp          main_qubit_SP.f90        -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(GF) -o  dressedqubit_sp   main_DressedQubit_SP.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)


############################
# CLEAN
############################

clean:
	rm qubit dressedqubit qubit_sp dressedqubit_sp dressedRb87_DRIVER spin_driver dressedqubit_DRIVER qubit_driver

