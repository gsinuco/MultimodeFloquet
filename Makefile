#export LD_LIBRARY_PATH="/opt/intel/compilers_and_libraries_2017/linux/mkl/lib/intel64"
# SET FORTRAN AND CPP COMPILERS
CPP = g++
CC  = gcc
GF  = gfortran
AR  = ar 
RANLIB = ranlib

# SET REQUIRED FLAGS
GFFLAGS    =  -llapack -lblas -g
GFFLAGS_SP =  -m64  -w -fno-second-underscore -x f77-cpp-input  -lpthread -lm -ldl  -llapack -lblas -g
MKLFLAGS   =  -lmkl_gf_lp64 -lmkl_sequential -lmkl_core
CFLAGS     =  -static

#SET MKL-intel LIBRARY PATH
MKLLIBS = /opt/intel/compilers_and_libraries/linux/mkl/lib/intel64
#SET MKL-intel INCLUDE PATH
MKLINC = /opt/intel/compilers_and_libraries/linux/mkl/include	

###################################
# MAKE LIBRARY AND ALL EXECUTABLES
###################################

all: all_examples lib lib_lapack

all_examples: Example_lib Example_lib_sp Example_lib_c Example_lib_c_sp


lib:build/modes.o build/Modules.o build/Modules_release.o build/delta_kr.o build/Floquet.o \
 build/I_and_J_representations.o build/F_representation.o build/LapackEigenValues.o \
 build/MKLSparseEigenValues.o build/util.o build/quick-sort-index-table.o build/VarCRSPacking.o \
 build/sparse_utils.o build/MultimodeHamiltonian_SP.o build/MultimodeHamiltonian.o \
 build/MultimodeFloquetTE.o build/MultimodeFloquetTE_DRIVER.o build/MultimodeMicroMotion.o \
 build/MultimodeTransitionAVG.o build/MultimodeDressedBasis.o build/MultimodeDressedBasis_SP.o \
 build/util_c.o build/modes_C.o build/Floquet_init_C.o build/MultimodeHamiltonian_SP_C.o  \
 build/MultimodeHamiltonian_C.o build/LapackEigenValues_C.o build/MultimodeTransitionAVG_C.o \
 build/MultimodeMicroMotion_C.o build/MultimodeFloquetTE_DRIVER_C.o build/MultimodeFloquetTE_C.o \
 build/MultimodeDressedBasis_C.o build/MultimodeDressedBasis_SP_C.o build/MKLSparseEigenValues_C.o 
	$(AR) -urv lib/libmultimodefloquet.a build/*.o
	$(RANLIB) lib/libmultimodefloquet.a
	cp *.mod ./include/

lib_lapack :build/modes.o build/Modules.o build/Modules_release.o build/delta_kr.o build/Floquet.o \
 build/I_and_J_representations.o build/F_representation.o build/LapackEigenValues.o \
 build/util.o build/quick-sort-index-table.o build/VarCRSPacking.o \
 build/sparse_utils.o build/MultimodeHamiltonian.o \
 build/MultimodeFloquetTE.o build/MultimodeFloquetTE_DRIVER.o build/MultimodeMicroMotion.o \
 build/MultimodeTransitionAVG.o build/MultimodeDressedBasis.o \
 build/util_c.o build/modes_C.o build/Floquet_init_C.o \
 build/MultimodeHamiltonian_C.o build/LapackEigenValues_C.o build/MultimodeTransitionAVG_C.o \
 build/MultimodeMicroMotion_C.o build/MultimodeFloquetTE_DRIVER_C.o build/MultimodeFloquetTE_C.o \
 build/MultimodeDressedBasis_C.o  build/MKLSparseEigenValues_C.o 
	$(AR) -urv lib/libmultimodefloquet.a build/*.o
	$(RANLIB) lib/libmultimodefloquet.a
	cp *.mod ./include/


Example_lib: ./examples/FORTRAN/main_qubit.f90  ./examples/FORTRAN/main_DressedQubit.f90
	$(GF) -o ./examples/FORTRAN/qubit  ./examples/FORTRAN/main_qubit.f90 -I./include/ -L./lib/ -lmultimodefloquet $(GFFLAGS)
	$(GF) -o ./examples/FORTRAN/dressedqubit  ./examples/FORTRAN/main_DressedQubit.f90 -I./include/ -L./lib/ -lmultimodefloquet $(GFFLAGS)

Example_lib_sp: ./examples/FORTRAN/main_qubit_SP.f90 ./examples/FORTRAN/main_DressedQubit_SP.f90
	$(GF) -o ./examples/FORTRAN/qubit_sp  ./examples/FORTRAN/main_qubit_SP.f90 -I./include/ -L./lib/ -lmultimodefloquet -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(GF) -o ./examples/FORTRAN/dressedqubit_sp  ./examples/FORTRAN/main_DressedQubit_SP.f90 -I./include/ -L./lib/ -lmultimodefloquet -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)

Example_lib_c: ./examples/CPP/main_qubit.cpp  ./examples/CPP/main_DressedQubit.cpp
	$(CPP) -o ./examples/CPP/qubit  ./examples/CPP/main_qubit.cpp -I./include/ -L./lib/ -lmultimodefloquet -lgfortran $(GFFLAGS)
	$(CPP) -o ./examples/CPP/dressedqubit  ./examples/CPP/main_DressedQubit.cpp -I./include/ -L./lib/ -lmultimodefloquet -lgfortran $(GFFLAGS)

Example_lib_c_sp: ./examples/CPP/main_qubit_sp.cpp ./examples/CPP/main_DressedQubit_SP.cpp
	$(CPP) -o  ./examples/CPP/qubit_sp         ./examples/CPP/main_qubit_sp.cpp        -I./include/ -L./lib/ -lmultimodefloquet -lgfortran -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)         
	$(CPP) -o  ./examples/CPP/dressedqubit_sp  ./examples/CPP/main_DressedQubit_SP.cpp -I./include/ -L./lib/ -lmultimodefloquet -lgfortran -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)


####################################
# BUILD OBJECT FILES FOR CPP WRAPPER
#####################################
build/util_c.o: build/util.o build/modes.o src/util_c.f90
	$(GF) -c -o $@ build/modes.o build/util.o src/util_c.f90 -g

build/modes_C.o: build/modes.o src/modes_C.f90
	$(GF) -c -o $@ build/modes.o src/modes_C.f90  -g

build/Floquet_init_C.o: build/modes.o build/modes_C.o src/Floquet_init_C.f90
	$(GF) -c -o $@ build/modes.o build/modes_C.o src/Floquet_init_C.f90  -g

build/MultimodeHamiltonian_SP_C.o: build/MultimodeHamiltonian_SP.o src/MultimodeHamiltonian_SP_C.f90
	$(GF) -c -o $@  src/MultimodeHamiltonian_SP_C.f90 -g

build/MultimodeHamiltonian_C.o: build/MultimodeHamiltonian.o src/MultimodeHamiltonian_C.f90
	$(GF) -c -o $@  build/MultimodeHamiltonian.o src/MultimodeHamiltonian_C.f90  -g

build/LapackEigenValues_C.o: build/Modules.o build/LapackEigenValues.o src/LapackEigenValues_C.f90
	$(GF) -c -o $@ build/Modules.o build/LapackEigenValues.o src/LapackEigenValues_C.f90  -g

build/MultimodeTransitionAVG_C.o: build/modes_C.o build/MultimodeTransitionAVG.o src/MultimodeTransitionAVG_C.f90
	$(GF) -c -o $@  build/modes_C.o build/MultimodeTransitionAVG.o src/MultimodeTransitionAVG_C.f90  -g

build/MultimodeMicroMotion_C.o: build/modes_C.o build/MultimodeMicroMotion.o src/MultimodeMicroMotion_C.f90
	$(GF) -c -o $@  build/modes_C.o build/MultimodeMicroMotion.o src/MultimodeMicroMotion_C.f90  -g

build/MultimodeFloquetTE_DRIVER_C.o: build/modes_C.o build/MultimodeFloquetTE_DRIVER.o src/MultimodeFloquetTE_DRIVER_C.f90
	$(GF) -c -o $@  build/modes_C.o build/MultimodeFloquetTE_DRIVER.o src/MultimodeFloquetTE_DRIVER_C.f90  -g

build/MultimodeFloquetTE_C.o: build/modes_C.o build/MultimodeFloquetTE.o src/MultimodeFloquetTE_C.f90
	$(GF) -c -o $@  build/modes_C.o build/MultimodeFloquetTE.o src/MultimodeFloquetTE_C.f90  -g

build/MultimodeDressedBasis_C.o: build/modes_C.o  build/MultimodeDressedBasis.o src/MultimodeDressedBasis_C.f90
	$(GF) -o $@ -c build/modes_C.o  build/MultimodeDressedBasis.o src/MultimodeDressedBasis_C.f90  -g

build/MultimodeDressedBasis_SP_C.o: build/modes_C.o  build/MultimodeDressedBasis_SP.o src/MultimodeDressedBasis_SP_C.f90
	$(GF) -o $@ -c build/modes_C.o  build/MultimodeDressedBasis_SP.o src/MultimodeDressedBasis_SP_C.f90  -g

build/MKLSparseEigenValues_C.o: build/MKLSparseEigenValues.o src/MKLSparseEigenvalues_C.f90
	$(GF) -c -o $@ src/MKLSparseEigenvalues_C.f90  -g


#build/util_c.o
#build/modes_C.o
#build/Floquet_init_C.o
#build/MultimodeHamiltonian_SP_C.o
#build/MultimodeHamiltonian_C.o
#build/LapackEigenValues_C.o
#build/MultimodeTransitionAVG_C.o
#build/MultimodeMicroMotion_C.o
#build/MultimodeFloquetTE_DRIVER_C.o
#build/MultimodeFloquetTE_C.o
#build/MultimodeDressedBasis_C.o
#build/MultimodeDressedBasis_SP_C.o
#build/MKLSparseEigenValues_C.o


############################
# BUILD FORTRAN OBJECT FILES
############################

build/modes.o: src/modes.f90
	$(GF) -c -o $@ src/modes.f90  

build/Modules.o: src/Modules.f90
	$(GF) -c -o $@ src/Modules.f90  -g

build/Modules_release.o: build/Modules.o src/Modules_release.f90
	$(GF) -c -o $@ src/Modules_release.f90  -g

build/delta_kr.o: src/delta_kr.f90
	$(GF) -c -o $@ src/delta_kr.f90  -g

build/Floquet.o: build/Modules.o build/modes.o src/Floquet_init.f90
	$(GF) -c -o $@ src/Floquet_init.f90  -g

build/I_and_J_representations.o: src/I_and_J_representations.f90
	$(GF) -c  -o $@ src/I_and_J_representations.f90  -g

build/F_representation.o: src/F_representation.f90
	$(GF) -c  -o $@ src/F_representation.f90 

build/LapackEigenValues.o:src/LapackEigenValues.f90
	$(GF) -c -o $@ src/LapackEigenValues.f90  -g

build/MKLSparseEigenValues.o:src/MKLSparseEigenvalues.f90
	$(GF) -c -o $@ src/MKLSparseEigenvalues.f90  -g

build/util.o: src/util.f90
	$(GF) -c -o $@ src/util.f90  -g

build/quick-sort-index-table.o: src/quick-sort-index-table.f90
	$(GF) -o $@ -c src/quick-sort-index-table.f90  -g

build/VarCRSPacking.o: src/VarCRSPacking.f90
	$(GF) -o $@ -c src/VarCRSPacking.f90  -g

build/sparse_utils.o:src/sparse_utils.f90
	$(GF) -o $@ -c src/sparse_utils.f90  -g

build/MultimodeHamiltonian_SP.o:src/MultimodeHamiltonian_SP.f90
	$(GF) -o $@ -c src/MultimodeHamiltonian_SP.f90  -g

build/MultimodeHamiltonian.o:src/MultimodeHamiltonian.f90
	$(GF) -o $@ -c src/MultimodeHamiltonian.f90 -g

build/MultimodeFloquetTE.o: build/MultimodeHamiltonian.o src/MultimodeFloquetTE.f90
	$(GF) -o $@ -c build/MultimodeHamiltonian.o src/MultimodeFloquetTE.f90 -g

build/MultimodeFloquetTE_DRIVER.o: build/MultimodeHamiltonian.o src/MultimodeFloquetTE_DRIVER.f90
	$(GF) -o $@ -c build/MultimodeHamiltonian.o src/MultimodeFloquetTE_DRIVER.f90  -g

build/MultimodeMicroMotion.o:src/MultimodeMicroMotion.f90
	$(GF) -o $@ -c src/MultimodeMicroMotion.f90  -g

build/MultimodeTransitionAVG.o:src/MultimodeTransitionAVG.f90
	$(GF) -o $@ -c src/MultimodeTransitionAVG.f90  -g

build/MultimodeDressedBasis.o:src/MultimodeDressedBasis.f90
	$(GF) -o $@ -c src/MultimodeDressedBasis.f90 -g

build/MultimodeDressedBasis_SP.o:src/MultimodeDressedBasis_SP.f90
	$(GF) -o $@ -c src/MultimodeDressedBasis_SP.f90  -g

build/MultimodeFloquet.o:src/MultimodeFloquet.f90
	$(GF) -o $@ -c src/MultimodeFloquet.f90  -g



#build/modes.o
#build/Modules.o
#build/Modules_release.o
#build/delta_kr.o
#build/Floquet.o
#build/I_and_J_representations.o
#build/F_representation.o
#build/LapackEigenValues.o
#build/MKLSparseEigenValues.o
#build/util.o
#build/quick-sort-index-table.o
#build/VarCRSPacking.o
#build/sparse_utils.o
#build/MultimodeHamiltonian_SP.o
#build/MultimodeHamiltonian.o
#build/MultimodeFloquetTE.o
#build/MultimodeFloquetTE_DRIVER.o
#build/MultimodeMicroMotion.o
#build/MultimodeTransitionAVG.o
#build/MultimodeDressedBasis.o
#build/MultimodeDressedBasis_SP.o
#build/MultimodeFloquet.o

############################
# CLEAN
############################

clean:
	rm build/*.o ./*mod build/MultimodeFloquet* include/* lib/*
