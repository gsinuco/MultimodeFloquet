aCC = gcc
GF = gfortran
# Use in zeeman
GFFLAGS =  -m64  -w -fno-second-underscore -x f77-cpp-input  -lpthread -lm -ldl -lgsl -lgslcblas -O3 -llapack -lblas -g
MKLFLAGS = -lmkl_gf_lp64 -lmkl_sequential -lmkl_core
CFLAGS =  -lgsl -lgslcblas -O3  -static

MKLLIBS = /opt/intel/compilers_and_libraries/linux/mkl/lib/intel64
MKLINC = /opt/intel/compilers_and_libraries/linux/mkl/include	

# Use in apollo:
#GFFLAGS = -L/usr/lib64/atlas/lib $(GSL_INCL) $(GSL_LIBS) -llapack -lblas -g -lgsl -lgslcblas -O3  
#CFLAGS =   $(GSL_INCL) $(GSL_LIBS) -lgsl -lgslcblas -O3


#all: build/MultimodeCoupling_IJ
#all: build/MultimodeFloquetBare
#all: build/MultimodeFloquetDressedQubit
all: build/MultimodeFloquetDressedQubit_SP build/MultimodeFloquetQubit build/MultimodeFloquetDressedQubit build/MultimodeFloquetRelease build/MultimodeFloquetQubit_SP
#Release


build/MultimodeFloquetDressedQubit_SP:  build/modes.o build/modes_release.o build/Modules.o build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o build/SetParameters.o build/VarCRSPacking.o build/quick-sort-index-table.o build/sparse_utils.o build/LapackEigenValues.o build/MultimodeFloquet.o build/MultimodeHamiltonian_SP.o src/main_DressedQubit_SP.f90 
	$(GF) -o $@ build/modes.o build/modes_release.o build/Modules.o  build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o build/SetParameters.o build/VarCRSPacking.o build/quick-sort-index-table.o build/LapackEigenValues.o build/sparse_utils.o build/MultimodeFloquet.o build/MultimodeHamiltonian_SP.o src/main_DressedQubit_SP.f90 -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS) $(MKLFLAGS)

build/MultimodeFloquetDressedQubit:  build/modes.o build/modes_release.o build/Modules.o build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o build/SetParameters.o build/VarCRSPacking.o build/quick-sort-index-table.o build/sparse_utils.o build/LapackEigenValues.o build/MultimodeFloquet.o src/main_DressedQubit.f90 
	$(GF) -o $@ build/modes.o build/modes_release.o build/Modules.o  build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o build/SetParameters.o build/VarCRSPacking.o build/quick-sort-index-table.o build/LapackEigenValues.o build/sparse_utils.o build/MultimodeFloquet.o src/main_DressedQubit.f90 -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS) $(MKLFLAGS)

build/MultimodeFloquetQubit_SP:  build/modes.o build/modes_release.o build/Modules.o build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o build/SetParameters.o build/VarCRSPacking.o build/quick-sort-index-table.o build/sparse_utils.o  build/LapackEigenValues.o build/MultimodeFloquet.o src/main_qubit_SP.f90 
	$(GF) -o $@ build/modes.o build/modes_release.o build/Modules.o  build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o build/SetParameters.o build/VarCRSPacking.o build/quick-sort-index-table.o build/sparse_utils.o  build/LapackEigenValues.o build/MultimodeFloquet.o  build/MultimodeHamiltonian_SP.o src/main_qubit_SP.f90  -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS) $(MKLFLAGS)


build/MultimodeFloquetQubit:  build/modes.o build/modes_release.o build/Modules.o build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o build/SetParameters.o build/VarCRSPacking.o build/quick-sort-index-table.o build/sparse_utils.o  build/LapackEigenValues.o build/MultimodeFloquet.o src/main_qubit.f90 
	$(GF) -o $@ build/modes.o build/modes_release.o build/Modules.o  build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o build/SetParameters.o build/VarCRSPacking.o build/quick-sort-index-table.o build/sparse_utils.o  build/LapackEigenValues.o build/MultimodeFloquet.o src/main_qubit.f90  -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS) $(MKLFLAGS)

build/MultimodeFloquetRelease:  build/modes.o build/modes_release.o build/Modules.o build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o build/SetParameters.o  build/VarCRSPacking.o build/quick-sort-index-table.o build/sparse_utils.o  build/LapackEigenValues.o build/MultimodeFloquet.o src/main_Release.f90 
	$(GF) -o $@ build/modes.o build/modes_release.o build/Modules.o  build/Modules_release.o build/util.o build/delta_kr.o build/I_and_J_representations.o build/F_representation.o build/Floquet.o build/SetParameters.o  build/VarCRSPacking.o build/quick-sort-index-table.o build/sparse_utils.o  build/LapackEigenValues.o build/MultimodeFloquet.o src/main_Release.f90  -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS) $(MKLFLAGS)

build/modes.o: src/modes.f90
	$(GF) -c -o $@ src/modes.f90 

build/modes_release.o: src/modes_release.f90
	$(GF) -c -o $@ src/modes_release.f90 $(GFFLAGS)

build/Modules.o: src/Modules.f90
	$(GF) -c -o $@ src/Modules.f90 $(GFFLAGS)

build/Modules_release.o: build/Modules.o src/Modules_release.f90
	$(GF) -c -o $@ src/Modules_release.f90 $(GFFLAGS)

build/delta_kr.o: src/delta_kr.f90
	$(GF) -c -o $@ src/delta_kr.f90 $(GFFLAGS)

build/Floquet.o: build/Modules.o build/modes.o src/Floquet_init.f90
	$(GF) -c -o $@ src/Floquet_init.f90 $(GFFLAGS)

build/I_and_J_representations.o: src/I_and_J_representations.f90
	$(GF) -c  -o $@ src/I_and_J_representations.f90 $(GFFLAGS) 

build/F_representation.o: src/F_representation.f90
	$(GF) -c  -o $@ src/F_representation.f90 $(GFFLAGS) 

build/IJtoF_TransformationMatrix.o:  src/ClebshGordan_gsl_fun.o src/IJtoF_TransformationMatrix.f90
	$(GF)  -c -o $@ $(GFFLAGS) src/IJtoF_TransformationMatrix.f90 

src/funciones.mod: build/Modules.o  src/delta_kr.f90 src/Modules.f90 	
	$(GF) -c  -o build/delta_kr.o src/delta_kr.f90 $(GFFLAGS) 

build/SetParameters.o: src/modes.f90 src/SetParameters.f90
	$(GF)  -c -o $@  src/SetParameters.f90 

build/PerturbativeShift.o: build/Modules.o src/funciones.mod src/PerturbativeShift.f90
	$(GF) -c  -o $@ src/PerturbativeShift.f90 $(GFFLAGS)  

build/RF-Landscape.o: build/Modules.o src/RF-Landscape.f90
	$(GF) -c -o $@ src/RF-Landscape.f90 $(GFFLAGS)

build/ClebshGordan_gsl_fun.o:	src/ClebshGordan_gsl_fun.c
	$(CC) $(CFLAGS) -c -o $@ src/ClebshGordan_gsl_fun.c

build/LapackEigenValues.o:src/LapackEigenValues.f90
	$(GF) -c -o $@ src/LapackEigenValues.f90 $(GFFLAGS)

build/util.o: src/util.f90
	$(GF) -c -o $@ src/util.f90 $(GFFLAGS)

build/MultimodeFloquet.o: src/Modules.f90 build/modes.o src/MultimodeFloquet.f90
	$(GF) -c -o $@ src/MultimodeFloquet.f90 

build/quick-sort-index-table.o: src/quick-sort-index-table.f90
	$(GF) -o $@ -c src/quick-sort-index-table.f90

build/VarCRSPacking.o: src/VarCRSPacking.f90
	$(GF) -o $@ -c src/VarCRSPacking.f90

build/sparse_utils.o:src/sparse_utils.f90
	$(GF) -o $@ -c src/sparse_utils.f90

build/MultimodeHamiltonian_SP.o:src/MultimodeHamiltonian_SP.f90
	$(GF) -o $@ -c src/MultimodeHamiltonian_SP.f90

clean:
	rm build/*.o ./*mod build/Multi*
