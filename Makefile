export LD_LIBRARY_PATH="/opt/intel/compilers_and_libraries_2017/linux/mkl/lib/intel64"

CPP = g++
CC = gcc
GF  = gfortran
AR  = ar 
RANLIB = ranlib
# Use in zeeman
GFFLAGS    =  -lgsl -lgslcblas -O3 -llapack -lblas -g
GFFLAGS_SP =  -m64  -w -fno-second-underscore -x f77-cpp-input  -lpthread -lm -ldl -lgsl -lgslcblas -O3 -llapack -lblas -g
MKLFLAGS   =  -lmkl_gf_lp64 -lmkl_sequential -lmkl_core
CFLAGS     =  -lgsl -lgslcblas -O3  -static

MKLLIBS = /opt/intel/compilers_and_libraries/linux/mkl/lib/intel64
MKLINC = /opt/intel/compilers_and_libraries/linux/mkl/include	

# Use in apollo:
#GFFLAGS = -L/usr/lib64/atlas/lib $(GSL_INCL) $(GSL_LIBS) -llapack -lblas -g -lgsl -lgslcblas -O3  
#CFLAGS =   $(GSL_INCL) $(GSL_LIBS) -lgsl -lgslcblas -O3


all: build/MultimodeFloquetDressedQubit_C build/MultimodeFloquetDressedQubit build/MultimodeFloquetQubit build/MultimodeFloquetQubit_C build/MultimodeFloquetQubit_SP build/MultimodeFloquetQubit_SP_C build/MultimodeFloquetDressedQubit_SP_C build/MultimodeFloquetDressedQubit_SP
	
#build/MultimodeFloquetQubit_SP Example_lib_sp build/MultimodeFloquetQubit_new Example_lib  build/MultimodeFloquetDressedQubit build/MultimodeFloquetDressedQubit_SP build/MultimodeFloquetCinterface

lib: build/modes.o build/Modules.o build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o build/VarCRSPacking.o build/quick-sort-index-table.o build/sparse_utils.o  build/LapackEigenValues.o build/MKLSparseEigenValues.o build/MultimodeFloquetTE.o  build/MultimodeMicroMotion.o build/MultimodeTransitionAVG.o build/MultimodeHamiltonian_SP.o build/MultimodeHamiltonian.o 
	$(AR) -urv lib/libmultimodefloquet.a build/*.o
	$(RANLIB) lib/libmultimodefloquet.a

Example_lib: ./examples/FORTRAN/main_DressedQubit.f90 
	$(GF) -o ./examples/FORTRAN/$@  ./examples/FORTRAN/main_DressedQubit.f90 -L./lib/ -lmultimodefloquet $(GFFLAGS)

#Example_lib_sp: ./examples/FORTRAN/main_DressedQubit_SP.f90 
#$(GF) -o ./examples/FORTRAN/$@  ./examples/FORTRAN/main_DressedQubit_SP.f90 -L./lib/ -lmultimodefloquet -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)

Example_lib_sp: ./examples/FORTRAN/main_qubit_SP.f90 
	$(GF) -o ./examples/FORTRAN/$@  ./examples/FORTRAN/main_qubit_SP.f90 -L./lib/ -lmultimodefloquet -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)




build/MultimodeFloquetQubit_SP_C:        build/modes.o build/modes_C.o build/Modules.o build/Modules_release.o build/util.o build/util_c.o build/sparse_utils.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/LapackEigenValues.o build/VarCRSPacking.o build/quick-sort-index-table.o build/Floquet.o build/Floquet_init_C.o build/MultimodeHamiltonian_SP.o build/MultimodeHamiltonian_SP_C.o build/MKLSparseEigenValues.o build/MKLSparseEigenValues_C.o build/MultimodeTransitionAVG.o build/MultimodeTransitionAVG_C.o  build/MultimodeMicroMotion.o build/MultimodeMicroMotion_C.o build/MultimodeFloquetTE.o build/MultimodeFloquetTE_C.o src/main_qubit_sp.cpp
	$(CPP) -o $@                     build/modes.o build/modes_C.o build/Modules.o build/Modules_release.o build/util.o build/util_c.o build/sparse_utils.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/LapackEigenValues.o build/VarCRSPacking.o build/quick-sort-index-table.o build/Floquet.o build/Floquet_init_C.o build/MultimodeHamiltonian_SP.o build/MultimodeHamiltonian_SP_C.o build/MKLSparseEigenValues.o build/MKLSparseEigenValues_C.o build/MultimodeTransitionAVG.o build/MultimodeTransitionAVG_C.o  build/MultimodeMicroMotion.o build/MultimodeMicroMotion_C.o build/MultimodeFloquetTE.o build/MultimodeFloquetTE_C.o src/main_qubit_sp.cpp -lgfortran -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)

build/MultimodeFloquetDressedQubit_SP_C: build/modes.o build/modes_C.o build/Modules.o build/Modules_release.o build/util.o build/util_c.o build/sparse_utils.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/LapackEigenValues.o build/VarCRSPacking.o build/quick-sort-index-table.o build/Floquet.o build/Floquet_init_C.o build/MultimodeHamiltonian_SP.o build/MultimodeHamiltonian_SP_C.o build/MKLSparseEigenValues.o build/MKLSparseEigenValues_C.o build/MultimodeTransitionAVG.o build/MultimodeTransitionAVG_C.o  build/MultimodeMicroMotion.o build/MultimodeMicroMotion_C.o build/MultimodeFloquetTE.o build/MultimodeFloquetTE_C.o build/MultimodeDressedBasis_SP.o build/MultimodeDressedBasis_SP_C.o src/main_DressedQubit_SP.cpp
	$(CPP) -o $@                     build/modes.o build/modes_C.o build/Modules.o build/Modules_release.o build/util.o build/util_c.o build/sparse_utils.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/LapackEigenValues.o build/VarCRSPacking.o build/quick-sort-index-table.o build/Floquet.o build/Floquet_init_C.o build/MultimodeHamiltonian_SP.o build/MultimodeHamiltonian_SP_C.o build/MKLSparseEigenValues.o build/MKLSparseEigenValues_C.o build/MultimodeTransitionAVG.o build/MultimodeTransitionAVG_C.o  build/MultimodeMicroMotion.o build/MultimodeMicroMotion_C.o build/MultimodeFloquetTE.o build/MultimodeFloquetTE_C.o build/MultimodeDressedBasis_SP.o build/MultimodeDressedBasis_SP_C.o src/main_DressedQubit_SP.cpp -lgfortran  -lgfortran -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)

build/MultimodeFloquetQubit_C: build/modes.o build/modes_C.o build/Modules.o build/Modules_release.o build/util.o build/util_c.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/LapackEigenValues.o build/Floquet.o build/Floquet_init_C.o build/MultimodeHamiltonian.o build/MultimodeHamiltonian_C.o build/LapackEigenValues_C.o  build/MultimodeTransitionAVG.o build/MultimodeTransitionAVG_C.o  build/MultimodeMicroMotion.o build/MultimodeMicroMotion_C.o build/MultimodeFloquetTE.o build/MultimodeFloquetTE_C.o src/main_qubit.cpp
	$(CPP) -o $@           build/modes.o build/modes_C.o build/Modules.o build/Modules_release.o build/util.o build/util_c.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/LapackEigenValues.o build/Floquet.o build/Floquet_init_C.o build/MultimodeHamiltonian.o build/MultimodeHamiltonian_C.o build/LapackEigenValues_C.o build/MultimodeTransitionAVG.o build/MultimodeTransitionAVG_C.o build/MultimodeMicroMotion.o build/MultimodeMicroMotion_C.o build/MultimodeFloquetTE.o build/MultimodeFloquetTE_C.o src/main_qubit.cpp -lgfortran   $(GFFLAGS) 

build/MultimodeFloquetDressedQubit_C: build/modes.o build/modes_C.o build/Modules.o build/Modules_release.o build/util.o build/util_c.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/LapackEigenValues.o build/Floquet.o build/Floquet_init_C.o build/MultimodeHamiltonian.o build/MultimodeHamiltonian_C.o build/LapackEigenValues_C.o  build/MultimodeTransitionAVG.o build/MultimodeTransitionAVG_C.o  build/MultimodeMicroMotion.o build/MultimodeMicroMotion_C.o build/MultimodeFloquetTE.o build/MultimodeFloquetTE_C.o build/MultimodeDressedBasis.o build/MultimodeDressedBasis_C.o src/main_DressedQubit.cpp
	$(CPP) -o $@                  build/modes.o build/modes_C.o build/Modules.o build/Modules_release.o build/util.o build/util_c.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/LapackEigenValues.o build/Floquet.o build/Floquet_init_C.o build/MultimodeHamiltonian.o build/MultimodeHamiltonian_C.o build/LapackEigenValues_C.o  build/MultimodeTransitionAVG.o build/MultimodeTransitionAVG_C.o  build/MultimodeMicroMotion.o build/MultimodeMicroMotion_C.o build/MultimodeFloquetTE.o build/MultimodeFloquetTE_C.o build/MultimodeDressedBasis.o build/MultimodeDressedBasis_C.o src/main_DressedQubit.cpp -lgfortran   $(GFFLAGS) 




build/util_c.o: build/util.o build/modes.o src/util_c.f90
	$(GF) -c -o $@ build/modes.o build/util.o src/util_c.f90 -g
	
build/modes_C.o: build/modes.o src/modes_C.f90
	$(GF) -c -o $@ build/modes.o src/modes_C.f90  -g

build/Floquet_init_C.o: build/modes.o build/modes_C.o src/Floquet_init_C.f90
	$(GF) -c -o $@ build/modes.o build/modes_C.o src/Floquet_init_C.f90 -O3 -g

build/MultimodeHamiltonian_SP_C.o: build/MultimodeHamiltonian_SP.o src/MultimodeHamiltonian_SP_C.f90
	$(GF) -c -o $@  src/MultimodeHamiltonian_SP_C.f90 -O3 -g

build/MultimodeHamiltonian_C.o: build/MultimodeHamiltonian.o src/MultimodeHamiltonian_C.f90
	$(GF) -c -o $@  build/MultimodeHamiltonian.o src/MultimodeHamiltonian_C.f90 -O3 -g

build/LapackEigenValues_C.o: build/Modules.o build/LapackEigenValues.o src/LapackEigenValues_C.f90
	$(GF) -c -o $@ build/Modules.o build/LapackEigenValues.o src/LapackEigenValues_C.f90 -O3 -g

build/MultimodeTransitionAVG_C.o: build/modes_C.o build/MultimodeTransitionAVG.o src/MultimodeTransitionAVG_C.f90
	$(GF) -c -o $@  build/modes_C.o build/MultimodeTransitionAVG.o src/MultimodeTransitionAVG_C.f90 -O3 -g

build/MultimodeMicroMotion_C.o: build/modes_C.o build/MultimodeMicroMotion.o src/MultimodeMicroMotion_C.f90
	$(GF) -c -o $@  build/modes_C.o build/MultimodeMicroMotion.o src/MultimodeMicroMotion_C.f90 -O3 -g

build/MultimodeFloquetTE_DRIVER_C.o: build/modes_C.o build/MultimodeFloquetTE_DRIVER.o src/MultimodeFloquetTE_DRIVER_C.f90
	$(GF) -c -o $@  build/modes_C.o build/MultimodeFloquetTE_DRIVER_.o src/MultimodeFloquetTE_DRIVER_C.f90 -O3 -g

build/MultimodeFloquetTE_C.o: build/modes_C.o build/MultimodeFloquetTE.o src/MultimodeFloquetTE_C.f90
	$(GF) -c -o $@  build/modes_C.o build/MultimodeFloquetTE.o src/MultimodeFloquetTE_C.f90 -O3 -g

build/MultimodeDressedBasis_C.o: build/modes_C.o  build/MultimodeDressedBasis.o src/MultimodeDressedBasis_C.f90
	$(GF) -o $@ -c build/modes_C.o  build/MultimodeDressedBasis.o src/MultimodeDressedBasis_C.f90 -O3 -g

build/MultimodeDressedBasis_SP_C.o: build/modes_C.o  build/MultimodeDressedBasis_SP.o src/MultimodeDressedBasis_SP_C.f90
	$(GF) -o $@ -c build/modes_C.o  build/MultimodeDressedBasis_SP.o src/MultimodeDressedBasis_SP_C.f90 -O3 -g

build/MKLSparseEigenValues_C.o: build/MKLSparseEigenValues.o src/MKLSparseEigenvalues_C.f90
	$(GF) -c -o $@ src/MKLSparseEigenvalues_C.f90 -O3 -g






build/MultimodeFloquetDressedQubit_SP:  build/modes.o  build/Modules.o build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o build/VarCRSPacking.o build/quick-sort-index-table.o build/LapackEigenValues.o build/sparse_utils.o build/MKLSparseEigenValues.o build/MultimodeFloquetTE.o build/MultimodeMicroMotion.o build/MultimodeTransitionAVG.o  build/MultimodeHamiltonian_SP.o build/MultimodeDressedBasis_SP.o src/main_DressedQubit_SP.f90 
	$(GF) -o $@                     build/modes.o  build/Modules.o build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o build/VarCRSPacking.o build/quick-sort-index-table.o build/LapackEigenValues.o  build/sparse_utils.o build/MKLSparseEigenValues.o build/MultimodeFloquetTE.o build/MultimodeMicroMotion.o build/MultimodeTransitionAVG.o  build/MultimodeHamiltonian_SP.o build/MultimodeDressedBasis_SP.o build/MultimodeHamiltonian.o src/main_DressedQubit_SP.f90 -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)

build/MultimodeFloquetDressedQubit:  build/modes.o build/Modules.o build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o  build/VarCRSPacking.o build/quick-sort-index-table.o build/LapackEigenValues.o build/MultimodeFloquetTE.o  build/MultimodeMicroMotion.o build/MultimodeTransitionAVG.o build/MultimodeHamiltonian.o src/main_DressedQubit.f90 
	$(GF) -o $@ build/modes.o build/Modules.o  build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o build/VarCRSPacking.o build/quick-sort-index-table.o build/LapackEigenValues.o  build/MultimodeFloquetTE.o  build/MultimodeMicroMotion.o build/MultimodeTransitionAVG.o build/MultimodeHamiltonian.o src/main_DressedQubit.f90  $(GFFLAGS) 

build/MultimodeFloquetQubit_SP:  build/modes.o build/Modules.o build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o build/VarCRSPacking.o build/quick-sort-index-table.o build/sparse_utils.o  build/LapackEigenValues.o build/MKLSparseEigenValues.o build/MultimodeFloquetTE.o  build/MultimodeMicroMotion.o build/MultimodeTransitionAVG.o build/MultimodeHamiltonian_SP.o src/main_qubit_SP.f90 
	$(GF) -o $@              build/modes.o build/Modules.o build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o  build/VarCRSPacking.o build/quick-sort-index-table.o build/sparse_utils.o  build/LapackEigenValues.o  build/MKLSparseEigenValues.o build/MultimodeFloquetTE.o build/MultimodeMicroMotion.o build/MultimodeTransitionAVG.o  build/MultimodeHamiltonian_SP.o src/main_qubit_SP.f90  -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)

build/MultimodeFloquetQubit:  build/modes.o build/Modules.o build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o build/VarCRSPacking.o build/quick-sort-index-table.o build/sparse_utils.o  build/LapackEigenValues.o build/MultimodeFloquetTE.o build/MultimodeMicroMotion.o build/MultimodeTransitionAVG.o  build/MultimodeHamiltonian.o src/main_qubit.f90 
	$(GF) -o $@ build/modes.o build/Modules.o  build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o build/VarCRSPacking.o build/quick-sort-index-table.o build/sparse_utils.o  build/LapackEigenValues.o build/MultimodeFloquetTE.o  build/MultimodeMicroMotion.o build/MultimodeTransitionAVG.o build/MultimodeHamiltonian.o src/main_qubit.f90  $(GFFLAGS) 

#build/MultimodeFloquetQubit_old:  build/modes.o build/Modules.o build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o build/VarCRSPacking.o build/quick-sort-index-table.o build/sparse_utils.o  build/LapackEigenValues.o build/MultimodeFloquet.o src/main_qubit.f90 
#	$(GF) -o $@ build/modes.o build/Modules.o  build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o build/VarCRSPacking.o build/quick-sort-index-table.o build/sparse_utils.o  build/LapackEigenValues.o build/MultimodeFloquet.o src/main_qubit.f90  $(GFFLAGS) 

#build/MultimodeFloquetQubit:  build/modes.o build/Modules.o build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o build/VarCRSPacking.o build/quick-sort-index-table.o build/sparse_utils.o  build/LapackEigenValues.o build/MultimodeFloquetTE.o build/MultimodeHamiltonian.o src/main_qubit.f90 
#	$(GF) -o $@ build/modes.o build/Modules.o  build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o build/VarCRSPacking.o build/quick-sort-index-table.o build/sparse_utils.o  build/LapackEigenValues.o build/MultimodeFloquetTE.o build/MultimodeHamiltonian.o src/main_qubit.f90  $(GFFLAGS) 


build/MultimodeFloquetRelease:  build/modes.o build/Modules.o build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o build/VarCRSPacking.o build/quick-sort-index-table.o build/sparse_utils.o  build/LapackEigenValues.o build/MKLSparseEigenValues.o build/MultimodeFloquetTE.o build/MultimodeMicroMotion.o build/MultimodeTransitionAVG.o build/MultimodeHamiltonian_SP.o build/MultimodeHamiltonian.o src/main_Release.f90 
	$(GF) -o $@ build/modes.o build/Modules.o build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o   build/VarCRSPacking.o build/quick-sort-index-table.o build/sparse_utils.o  build/LapackEigenValues.o  build/MKLSparseEigenValues.o build/MultimodeFloquetTE.o  build/MultimodeMicroMotion.o build/MultimodeTransitionAVG.o build/MultimodeHamiltonian_SP.o build/MultimodeHamiltonian.o src/main_Release.f90  -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)




build/modes.o: src/modes.f90
	$(GF) -c -o $@ src/modes.f90  

build/Modules.o: src/Modules.f90
	$(GF) -c -o $@ src/Modules.f90 -O3 -g

build/Modules_release.o: build/Modules.o src/Modules_release.f90
	$(GF) -c -o $@ src/Modules_release.f90 -O3 -g

build/delta_kr.o: src/delta_kr.f90
	$(GF) -c -o $@ src/delta_kr.f90 -O3 -g

build/Floquet.o: build/Modules.o build/modes.o src/Floquet_init.f90
	$(GF) -c -o $@ src/Floquet_init.f90 -O3 -g

build/I_and_J_representations.o: src/I_and_J_representations.f90
	$(GF) -c  -o $@ src/I_and_J_representations.f90 -O3 -g

build/F_representation.o: src/F_representation.f90
	$(GF) -c  -o $@ src/F_representation.f90 -O3

build/LapackEigenValues.o:src/LapackEigenValues.f90
	$(GF) -c -o $@ src/LapackEigenValues.f90 -O3 -g

build/MKLSparseEigenValues.o:src/MKLSparseEigenvalues.f90
	$(GF) -c -o $@ src/MKLSparseEigenvalues.f90 -O3 -g

build/util.o: src/util.f90
	$(GF) -c -o $@ src/util.f90 -O3 -g

build/quick-sort-index-table.o: src/quick-sort-index-table.f90
	$(GF) -o $@ -c src/quick-sort-index-table.f90 -O3 -g

build/VarCRSPacking.o: src/VarCRSPacking.f90
	$(GF) -o $@ -c src/VarCRSPacking.f90 -O3 -g

build/sparse_utils.o:src/sparse_utils.f90
	$(GF) -o $@ -c src/sparse_utils.f90 -O3 -g

build/MultimodeHamiltonian_SP.o:src/MultimodeHamiltonian_SP.f90
	$(GF) -o $@ -c src/MultimodeHamiltonian_SP.f90 -O3 -g

build/MultimodeHamiltonian.o:src/MultimodeHamiltonian.f90
	$(GF) -o $@ -c src/MultimodeHamiltonian.f90 -O3 -g

build/MultimodeFloquetTE.o: build/MultimodeHamiltonian.o src/MultimodeFloquetTE.f90
	$(GF) -o $@ -c build/MultimodeHamiltonian.o src/MultimodeFloquetTE.f90 -O3 -g

build/MultimodeFloquetTE_DRIVER.o: build/MultimodeHamiltonian.o src/MultimodeFloquetTE_DRIVER.f90
	$(GF) -o $@ -c build/MultimodeHamiltonian.o src/MultimodeFloquetTE_DRIVER.f90 -O3 -g

build/MultimodeMicroMotion.o:src/MultimodeMicroMotion.f90
	$(GF) -o $@ -c src/MultimodeMicroMotion.f90 -O3 -g

build/MultimodeTransitionAVG.o:src/MultimodeTransitionAVG.f90
	$(GF) -o $@ -c src/MultimodeTransitionAVG.f90 -O3 -g

build/MultimodeDressedBasis.o:src/MultimodeDressedBasis.f90
	$(GF) -o $@ -c src/MultimodeDressedBasis.f90 -O3 -g

build/MultimodeDressedBasis_SP.o:src/MultimodeDressedBasis_SP.f90
	$(GF) -o $@ -c src/MultimodeDressedBasis_SP.f90 -O3 -g
	
build/MultimodeFloquet.o:src/MultimodeFloquet.f90
	$(GF) -o $@ -c src/MultimodeFloquet.f90 -O3 -g

clean:
	rm build/*.o ./*mod build/Multi*
