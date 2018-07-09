!export LD_LIBRARY_PATH="/opt/intel/compilers_and_libraries_2017/linux/mkl/lib/intel64";
PROGRAM MULTIMODEFLOQUET

  USE ATOMIC_PROPERTIES
  USE TYPES
  USE SPARSE_INTERFACE
  USE SUBINTERFACE
  USE SUBINTERFACE_LAPACK
  USE FLOQUETINIT_ 
  USE ARRAYS 

  IMPLICIT NONE
  TYPE(MODE),       DIMENSION(:),   ALLOCATABLE :: FIELDS
  TYPE(ATOM)                                       ID
  INTEGER,          DIMENSION(:),   ALLOCATABLE :: MODES_NUM
  INTEGER                                       :: TOTAL_FREQUENCIES
  INTEGER                                       :: INFO,m,INDEX0
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: ENERGY,E_FLOQUET
  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: H__,U_F,U_AUX,U_B2D,U_F1,U_F2,U_F1_red,U_F2_red
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: P_AVG
  DOUBLE PRECISION                              :: T1,T2

  
  INTEGER :: FIELD_INDEX

  !PARAMETERS NEEDED TO DEFINE THE SPARSE MATRIX
  INTEGER,    DIMENSION(:), ALLOCATABLE :: ROW_INDEX,COLUMN
  COMPLEX*16, DIMENSION(:), ALLOCATABLE :: VALUES
  !PARAMETERS TO DIAGONALIZE THE SPARSE MATRIX WITH MKL
  DOUBLE PRECISION :: E_L,E_R
  INTEGER          :: D_MULTIFLOQUET,r

  INFO = 0
  CALL FLOQUETINIT('qubit','U',2,ID,INFO)
  ALLOCATE(ENERGY(SIZE(J_Z,1)))
  ALLOCATE(H__(SIZE(J_z,1),SIZE(J_z,1)))
  ALLOCATE(P_AVG(SIZE(J_z,1),SIZE(J_z,1)))
  ALLOCATE(U_AUX(SIZE(J_z,1),SIZE(J_z,1)))
  H__ = J_z
  CALL LAPACK_FULLEIGENVALUES(H__,TOTAL_STATES_LSI,Energy,INFO)
  DEALLOCATE(H__)
  DEALLOCATE(ENERGY)
  
  
  ALLOCATE(MODES_NUM(3))

  MODES_NUM(1) = 1 !(STATIC FIELD)
  MODES_NUM(2) = 2 !(DRIVING BY TWO HARMONICS)
  MODES_NUM(3) = 2 !(DRIVING BY TWO HARMONICS)
!  MODES_NUM(4) = 1 !(DRIVING BY TWO HARMONICS)
  
  TOTAL_FREQUENCIES = SUM(MODES_NUM,1)
  ALLOCATE(FIELDS(TOTAL_FREQUENCIES))
  DO m=1,TOTAL_FREQUENCIES
     ALLOCATE(FIELDS(m)%V(ID%D_BARE,ID%D_BARE))
  END DO
  
  FIELDS(1)%X    = 0.0
  FIELDS(1)%Y    = 0.0
  FIELDS(1)%Z    = 2.0
  FIELDS(1)%phi_x = 0.0
  FIELDS(1)%phi_y = 0.0
  FIELDS(1)%phi_z = 0.0
  FIELDS(1)%omega = 0.0
  FIELDS(1)%N_Floquet = 0

  FIELDS(2)%X     = 4.0
  FIELDS(2)%Y     = 0.0
  FIELDS(2)%Z     = 0.0
  FIELDS(2)%phi_x = 0.0
  FIELDS(2)%phi_y = 0.0
  FIELDS(2)%phi_z = 0.0
  FIELDS(2)%omega = 2.0
  FIELDS(2)%N_Floquet = 3

!!$  FIELDS(3)%X     = 8.0
!!$  FIELDS(3)%y     = 0.0
!!$  FIELDS(3)%z     = 0.0
!!$  FIELDS(3)%phi_x = 0.0
!!$  FIELDS(3)%phi_y = 0.0
!!$  FIELDS(3)%phi_z = 0.0
!!$  FIELDS(3)%omega = 10
!!$  FIELDS(3)%N_Floquet = 1
!!$

  FIELDS(3)%X     = 8.0
  FIELDS(3)%y     = 0.0
  FIELDS(3)%z     = 0.0
  FIELDS(3)%phi_x = 0.0
  FIELDS(3)%phi_y = 0.0
  FIELDS(3)%phi_z = 0.0
  FIELDS(3)%omega = 2.0
  FIELDS(3)%N_Floquet = 1


  FIELDS(4)%X     = 8.0
  FIELDS(4)%y     = 0.0
  FIELDS(4)%z     = 0.0
  FIELDS(4)%phi_x = 0.0
  FIELDS(4)%phi_y = 0.0
  FIELDS(4)%phi_z = 0.0
  FIELDS(4)%omega = 10.0
  FIELDS(4)%N_Floquet = 3


!!$  FIELDS(4)%X    = 10.0
!!$  FIELDS(4)%Y    = 0.0
!!$  FIELDS(4)%Z    = 0.0
!!$  FIELDS(4)%phi_x = 0.0
!!$  FIELDS(4)%phi_y = 0.0
!!$  FIELDS(4)%phi_z = 0.0
!!$  FIELDS(4)%omega = 10.0
!!$  FIELDS(4)%N_Floquet = 1

  FIELDS(5)%X    = 12.0
  FIELDS(5)%Y    =  0.0
  FIELDS(5)%Z    =  0.0
  FIELDS(5)%phi_x = 0.0
  FIELDS(5)%phi_y = 0.0
  FIELDS(5)%phi_z = 0.0
  FIELDS(5)%omega = 20.0
  FIELDS(5)%N_Floquet = 1
!!$
!!$  FIELDS(6)%X    = 1.0
!!$  FIELDS(6)%Y    = 1.0
!!$  FIELDS(6)%Z    = 0.0
!!$  FIELDS(6)%phi_x = 0.0
!!$  FIELDS(6)%phi_y = 0.0
!!$  FIELDS(6)%phi_z = 0.0
!!$  FIELDS(6)%omega = 30.0
!!$  FIELDS(6)%N_Floquet = 1
!!$
!!$  FIELDS(7)%X    = 1.0
!!$  FIELDS(7)%Y    = 1.0
!!$  FIELDS(7)%Z    = 1.0
!!$  FIELDS(7)%phi_x = 0.0
!!$  FIELDS(7)%phi_y = 0.0
!!$  FIELDS(7)%phi_z = 0.0
!!$  FIELDS(7)%omega = 100.0
!!$  FIELDS(7)%N_Floquet = 1

  D_MULTIFLOQUET = ID%D_BARE
  DO r=2,size(MODES_NUM)
     FIELD_INDEX = 2+SUM(MODES_NUM(2:r-1))
     D_MULTIFLOQUET = D_MULTIFLOQUET*(2*FIELDS(FIELD_INDEX)%N_Floquet+1)
  END DO
  write(*,*) D_MULTIFLOQUET

  CALL SETHAMILTONIANCOMPONENTS(ID,size(modes_num,1),total_frequencies,MODES_NUM,FIELDS,INFO)
!  call write_matrix(real(fields(1)%V))
!  call write_matrix(real(fields(2)%V))
  !--- FIND THE MULTIMODE FLOQUET SPECTRUM 
  CALL MULTIMODEFLOQUETMATRIX_SP(ID,SIZE(MODES_NUM,1),total_frequencies,MODES_NUM,FIELDS,VALUES,ROW_INDEX,COLUMN,INFO)
!  write(*,*) real(Values)
!  write(*,*) size(values,1), size(row_index,1), size(column,1)
  E_L = -40.0
  E_R =  40.0
  ALLOCATE(E_FLOQUET(D_MULTIFLOQUET))
  ALLOCATE(U_F(D_MULTIFLOQUET,D_MULTIFLOQUET))
  E_FLOQUET = 0.0
  U_F = 0.0
  !  WRITE(*,*)D_MULTIFLOQUET,SIZE(VALUES,1),SIZE(ROW_INDEX,1),SIZE(COLUMN,1),SIZE(U_F,1),SIZE(U_F,2)
  !  WRITE(*,*) REAL(VALUES)!,ROW_INDEX,COLUMN,E_L,E_R!,E_FLOQUET,U_F,INFO
  !DO r=1,SIZE(VALUES,1)
  !   WRITE(*,*) r,ROW_INDEX(r),COLUMN(r),REAL(VALUES(r))
  !END DO
  
  CALL MKLSPARSE_FULLEIGENVALUES(D_MULTIFLOQUET,SIZE(VALUES,1),VALUES,ROW_INDEX,COLUMN,E_L,E_R,E_FLOQUET,U_F,INFO)
  write(*,*) E_FLOQUET

!  CALL MULTIMODEFLOQUETMATRIX(ID,size(modes_num,1),total_frequencies,MODES_NUM,FIELDS,INFO)
!  call write_matrix(real(H_floquet))
!  WRITE(*,*)
!  WRITE(*,*)
!  ALLOCATE(E_FLOQUET(SIZE(H_FLOQUET,1)))
!  E_FLOQUET = 0.0  
!  CALL LAPACK_FULLEIGENVALUES(H_FLOQUET,SIZE(H_FLOQUET,1),E_FLOQUET,INFO)
!  U_F = H_FLOQUET ! FOURIER DECOMPOSITION OF THE STATES DRESSED BY MODE NUMBER 
!  DEALLOCATE(H_FLOQUET)
  !write(*,*) e_floquet
  ! Evaluate avg transition probabilities in the bare basis
  P_AVG = 0.0
  CALL MULTIMODETRANSITIONAVG(SIZE(U_F,1),size(MODES_NUM,1),FIELDS,MODES_NUM,U_F,E_FLOQUET,ID%D_BARE,P_AVG,INFO) 
  !CALL WRITE_MATRIX(P_AVG)
  
  !---  EVALUATE INSTANTANEOUS MULTIMODE FLOQUET TRANSFORMATION
  ALLOCATE(U_B2D(ID%D_BARE,SIZE(U_F,1)))
  T1 = 0.0
  write(*,*) size(modes_num,1)
  CALL MULTIMODEFLOQUETTRANSFORMATION(SIZE(U_F,1), SIZE(MODES_NUM,1),MODES_NUM,U_F,E_FLOQUET,ID%D_BARE,FIELDS,T1,U_B2D,INFO) 

  !--- EVALUATE TIME-EVOLUTION OPERATOR IN THE BARE BASIS
  T1 = 0.0
  DO r=1,1!256
     T2 = r*32.0*4.0*atan(1.0)/256
     CALL MULTIMODETIMEEVOLUTINOPERATOR(SIZE(U_F,1),SIZE(MODES_NUM,1),MODES_NUM,U_F,E_FLOQUET,ID%D_BARE,FIELDS,T1,T2,U_AUX,INFO) 
     CALL WRITE_MATRIX(ABS(U_AUX)**2)
     !        WRITE(*,*) t2,FIELDS(2)%OMEGA,ABS(U_AUX)**2
  END DO
  WRITE(*,*)
  
END PROGRAM MULTIMODEFLOQUET

