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
MMFLIB = ../../../lib/
#SET MULTIMODEFLOQUET INCLUDE PATH
MMFINC = ../../../include/

###################################
# MAKE CPP EXAMPLES
###################################

all: Example_lib_f 



Example_lib_f: main_87Rb_RFMW_SPECTRUM.f90 main_87Rb_RFMW_SPECTRUM_SP.f90 main_87Rb_RFx2_SPECTRUM.f90 main_87Rb_RFx2_SPECTRUM_SP.f90 main_87Rb_RFx2_SPECTRUM_B.f90 main_87Rb_RFx1_SPECTRUM.f90 main_87Rb_RFx2_SPECTRUM_C.f90
	$(GF) -o dressedRb87_SPECTRUM    main_87Rb_RFMW_SPECTRUM.f90    -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) 
	$(GF) -o dressedRb87_SPECTRUM_SP main_87Rb_RFMW_SPECTRUM_SP.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(GF) -o dressedRb87RFx2_SPECTRUM    main_87Rb_RFx2_SPECTRUM.f90    -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) 
	$(GF) -o dressedRb87RFx2_SPECTRUM_B    main_87Rb_RFx2_SPECTRUM_B.f90    -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) 
	$(GF) -o dressedRb87RFx2_SPECTRUM_C    main_87Rb_RFx2_SPECTRUM_C.f90    -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) 
	$(GF) -o dressedRb87RFx2_SPECTRUM_SP main_87Rb_RFx2_SPECTRUM_SP.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(GF) -o dressedRb87RFx1_SPECTRUM    main_87Rb_RFx1_SPECTRUM.f90    -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) 

############################
# CLEAN
############################

clean:
	rm dressedRb87_SPECTRUM dressedRb87_SPECTRUM_SP dressedRb87RFx2_SPECTRUM dressedRb87RFx2_SPECTRUM_B dressedRb87RFx2_SPECTRUM_C dressedRb87RFx2_SPECTRUM_SP dressedRb87RFx1_SPECTRUM *.*~

