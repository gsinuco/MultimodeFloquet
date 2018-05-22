aCC = gcc
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


all: build/MultimodeFloquetDressedQubit_SP build/MultimodeFloquetQubit build/MultimodeFloquetDressedQubit build/MultimodeFloquetRelease build/MultimodeFloquetQubit_SP	
lib: build/modes.o build/Modules.o build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o build/VarCRSPacking.o build/quick-sort-index-table.o build/sparse_utils.o  build/LapackEigenValues.o build/MKLSparseEigenValues.o build/MultimodeFloquetTE.o build/MultimodeHamiltonian_SP.o build/MultimodeHamiltonian.o 
	$(AR) -urv lib/libmultimodefloquet.a build/*.o
	$(RANLIB) lib/libmultimodefloquet.a

Example_lib: ./examples/FORTRAN/main_DressedQubit.f90 
	$(GF) -o ./examples/FORTRAN/$@  ./examples/FORTRAN/main_DressedQubit.f90 -L./lib/ -lmultimodefloquet $(GFFLAGS)

Example_lib_sp: ./examples/FORTRAN/main_DressedQubit_SP.f90 
	$(GF) -o ./examples/FORTRAN/$@  ./examples/FORTRAN/main_DressedQubit_SP.f90 -L./lib/ -lmultimodefloquet -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)

build/MultimodeFloquetDressedQubit_SP:  build/modes.o  build/Modules.o build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o build/VarCRSPacking.o build/quick-sort-index-table.o build/LapackEigenValues.o build/sparse_utils.o build/MKLSparseEigenValues.o build/MultimodeFloquetTE.o build/MultimodeHamiltonian.o build/MultimodeHamiltonian_SP.o src/main_DressedQubit_SP.f90 
	$(GF) -o $@ build/modes.o  build/Modules.o  build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o build/VarCRSPacking.o build/quick-sort-index-table.o build/LapackEigenValues.o  build/sparse_utils.o build/MKLSparseEigenValues.o build/MultimodeFloquetTE.o build/MultimodeHamiltonian.o build/MultimodeHamiltonian_SP.o src/main_DressedQubit_SP.f90 -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)

build/MultimodeFloquetDressedQubit:  build/modes.o build/Modules.o build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o  build/VarCRSPacking.o build/quick-sort-index-table.o build/LapackEigenValues.o build/MultimodeFloquetTE.o build/MultimodeHamiltonian.o src/main_DressedQubit.f90 
	$(GF) -o $@ build/modes.o build/Modules.o  build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o build/VarCRSPacking.o build/quick-sort-index-table.o build/LapackEigenValues.o  build/MultimodeFloquetTE.o build/MultimodeHamiltonian.o src/main_DressedQubit.f90  $(GFFLAGS) 

build/MultimodeFloquetQubit_SP:  build/modes.o build/Modules.o build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o build/VarCRSPacking.o build/quick-sort-index-table.o build/sparse_utils.o  build/LapackEigenValues.o build/MKLSparseEigenValues.o build/MultimodeFloquetTE.o build/MultimodeHamiltonian_SP.o src/main_qubit_SP.f90 
	$(GF) -o $@ build/modes.o  build/Modules.o  build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o  build/VarCRSPacking.o build/quick-sort-index-table.o build/sparse_utils.o  build/LapackEigenValues.o  build/MKLSparseEigenValues.o build/MultimodeFloquetTE.o build/MultimodeHamiltonian_SP.o src/main_qubit_SP.f90  -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)

build/MultimodeFloquetQubit:  build/modes.o build/Modules.o build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o build/VarCRSPacking.o build/quick-sort-index-table.o build/sparse_utils.o  build/LapackEigenValues.o build/MultimodeFloquetTE.o build/MultimodeHamiltonian.o src/main_qubit.f90 
	$(GF) -o $@ build/modes.o build/Modules.o  build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o build/VarCRSPacking.o build/quick-sort-index-table.o build/sparse_utils.o  build/LapackEigenValues.o build/MultimodeFloquetTE.o build/MultimodeHamiltonian.o src/main_qubit.f90  $(GFFLAGS) 

build/MultimodeFloquetRelease:  build/modes.o build/Modules.o build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o build/VarCRSPacking.o build/quick-sort-index-table.o build/sparse_utils.o  build/LapackEigenValues.o build/MKLSparseEigenValues.o build/MultimodeFloquetTE.o build/MultimodeHamiltonian_SP.o build/MultimodeHamiltonian.o src/main_Release.f90 
	$(GF) -o $@ build/modes.o build/Modules.o build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o   build/VarCRSPacking.o build/quick-sort-index-table.o build/sparse_utils.o  build/LapackEigenValues.o  build/MKLSparseEigenValues.o build/MultimodeFloquetTE.o build/MultimodeHamiltonian_SP.o build/MultimodeHamiltonian.o src/main_Release.f90  -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)

build/modes.o: src/modes.f90
	$(GF) -c -o $@ src/modes.f90 -O3

build/Modules.o: src/Modules.f90
	$(GF) -c -o $@ src/Modules.f90 -O3

build/Modules_release.o: build/Modules.o src/Modules_release.f90
	$(GF) -c -o $@ src/Modules_release.f90 -O3

build/delta_kr.o: src/delta_kr.f90
	$(GF) -c -o $@ src/delta_kr.f90 -O3

build/Floquet.o: build/Modules.o build/modes.o src/Floquet_init.f90
	$(GF) -c -o $@ src/Floquet_init.f90 -O3 

build/I_and_J_representations.o: src/I_and_J_representations.f90
	$(GF) -c  -o $@ src/I_and_J_representations.f90 -O3

build/F_representation.o: src/F_representation.f90
	$(GF) -c  -o $@ src/F_representation.f90 -O3

build/LapackEigenValues.o:src/LapackEigenValues.f90
	$(GF) -c -o $@ src/LapackEigenValues.f90 -O3

build/MKLSparseEigenValues.o:src/MKLSparseEigenvalues.f90
	$(GF) -c -o $@ src/MKLSparseEigenvalues.f90 -O3

build/util.o: src/util.f90
	$(GF) -c -o $@ src/util.f90 -O3

build/quick-sort-index-table.o: src/quick-sort-index-table.f90
	$(GF) -o $@ -c src/quick-sort-index-table.f90 -O3

build/VarCRSPacking.o: src/VarCRSPacking.f90
	$(GF) -o $@ -c src/VarCRSPacking.f90

build/sparse_utils.o:src/sparse_utils.f90
	$(GF) -o $@ -c src/sparse_utils.f90 -O3

build/MultimodeHamiltonian_SP.o:src/MultimodeHamiltonian_SP.f90
	$(GF) -o $@ -c src/MultimodeHamiltonian_SP.f90 -O3

build/MultimodeHamiltonian.o:src/MultimodeHamiltonian.f90
	$(GF) -o $@ -c src/MultimodeHamiltonian.f90 -O3

build/MultimodeFloquetTE.o:src/MultimodeFloquetTE.f90
	$(GF) -o $@ -c src/MultimodeFloquetTE.f90 -O3

clean:
	rm build/*.o ./*mod build/Multi*
