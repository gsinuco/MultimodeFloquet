
MODULE FLOQUETINIT2
  implicit none
  interface
     !  CALL FLOQUETINIT('qubit',id,'U',0.2D1,INFO)  

     subroutine testopt(one,two)
!       USE TYPES
       implicit none
!       TYPE(ATOM),       intent(inout)       :: ID
!       CHARACTER(LEN=*), INTENT(IN)          :: spin
!       character(len=*), intent(in),optional :: qubit
       double precision, intent(in),optional :: one
       integer,          intent(inout)       :: two
!       integer, intent(out) :: three
!       integer, intent(in), optional :: four
!       integer, intent(out), optional :: five
     end subroutine testopt
     SUBROUTINE TESTFLOQUETINIT(atomicspecie,one,two,three,four,five)
       IMPLICIT NONE
       
       CHARACTER(LEN=*), INTENT(IN),optional    :: ATOMICSPECIE
       integer, intent(in), optional :: one,two,three,four,five
     end subroutine TESTFLOQUETINIT
  end interface
END MODULE FLOQUETINIT2

!MODULE FLOQUETINIT3
!  IMPLICIT NONE
!CONTAINS
  subroutine testopt(one,two)
!    USE TYPES
    implicit none
!    TYPE(ATOM),       intent(inout)       :: ID
!    CHARACTER(LEN=*), INTENT(IN)          :: spin
!    character(len=*), intent(in),optional :: qubit
    double precision, intent(in),optional :: one
    integer,          intent(inout)       :: two

 !   if(present(qubit)) write(*,*) qubit
  end subroutine testopt
!END MODULE FLOQUETINIT3


PROGRAM MULTIMODEFLOQUET

  USE ATOMIC_PROPERTIES
  USE TYPES
  USE SUBINTERFACE
  USE SUBINTERFACE_LAPACK
  USE ARRAYS 
  USE FLOQUETINITINTERFACE
!  USE FLOQUETINIT2
 ! USE FLOQUETINIT3

  IMPLICIT NONE
  TYPE(MODE),       DIMENSION(:),   ALLOCATABLE :: FIELDS
  TYPE(ATOM)                                       ID
  INTEGER,          DIMENSION(:),   ALLOCATABLE :: MODES_NUM
  INTEGER                                          TOTAL_FREQUENCIES,D_BARE
  INTEGER                                          INFO,m,INDEX0,r
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: ENERGY,E_FLOQUET
  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: H__,U_F,U_AUX,U_B2D
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: P_AVG
  DOUBLE PRECISION                              :: T1,T2


  INTEGER :: one,two,three,four,five

!  OPEN(UNIT=3,FILE="qubit_oscillation_.dat",ACTION="WRITE")
!  OPEN(UNIT=4,FILE="qubit_avgtransition_.dat",ACTION="WRITE")



  INFO = 0
  CALL FLOQUETINIT(INFO,ID,'qubit')  
  CALL DEALLOCATEALL(ID)
  CALL FLOQUETINIT(INFO,ID,'Rb87','B')  
  CALL DEALLOCATEALL(ID)
  CALL FLOQUETINIT(INFO,ID,'spin')  
!  CALL FLOQUETINIT('qubit',ID,INFO)  
! CALL TESTFLOQUETINIT('spin',one,two,three,four,five)  
!    CALL TESTFLOQUETINIT(one,two,three,four,five)  
 
!  CALL FLOQUETINIT('qubit','U',ID,INFO)  
!  CALL testopt('spin',ID,one,two,'qubit')
!  CALL testopt('',ID,one,two)
!  CALL testopt(0.2D1,INFO)
!  CALL testopt(INFO)
!  CALL testopt(one,two,three)

!!$  D_BARE = ID%D_BARE
!!$  ALLOCATE(P_AVG(D_BARE,D_BARE))
!!$  ALLOCATE(U_AUX(D_BARE,D_BARE))
!!$!  ALLOCATE(U_B2D(D_BARE,D_BARE))
!!$
!!$  ALLOCATE(MODES_NUM(2))
!!$  
!!$  MODES_NUM(1) = 1 !(STATIC FIELD)
!!$  MODES_NUM(2) = 1 !(DRIVING BY ONE HARMONIC)
!!$  
!!$  TOTAL_FREQUENCIES = SUM(MODES_NUM,1)
!!$  ALLOCATE(FIELDS(TOTAL_FREQUENCIES))
!!$  DO m=1,TOTAL_FREQUENCIES
!!$     ALLOCATE(FIELDS(m)%V(ID%D_BARE,ID%D_BARE))
!!$  END DO
!!$  
!!$  FIELDS(1)%X    = 0.0
!!$  FIELDS(1)%Y    = 0.0
!!$  FIELDS(1)%Z    = 1.0
!!$  FIELDS(1)%phi_x = 0.0
!!$  FIELDS(1)%phi_y = 0.0
!!$  FIELDS(1)%phi_z = 0.0
!!$  FIELDS(1)%omega = 0.0
!!$  FIELDS(1)%N_Floquet = 0
!!$
!!$  FIELDS(2)%X     = 2.0
!!$  FIELDS(2)%Y     = 0.0
!!$  FIELDS(2)%Z     = 0.0
!!$  FIELDS(2)%phi_x = 0.0
!!$  FIELDS(2)%phi_y = 0.0
!!$  FIELDS(2)%phi_z = 0.0
!!$  FIELDS(2)%omega = 1.0
!!$  FIELDS(2)%N_Floquet = 20
!!$
!!$  DO m=1,1!28
!!$
!!$     ! --- SET DRIVING PARAMETERS 
!!$     FIELDS(2)%omega = 0.2 + (m-1)*2.0/128
!!$     CALL SETHAMILTONIANCOMPONENTS(ID,size(modes_num,1),total_frequencies,MODES_NUM,FIELDS,INFO)
!!$     
!!$     !--- FIND THE MULTIMODE FLOQUET SPECTRUM 
!!$     CALL MULTIMODEFLOQUETMATRIX(ID,size(modes_num,1),total_frequencies,MODES_NUM,FIELDS,INFO)
!!$     ALLOCATE(E_FLOQUET(SIZE(H_FLOQUET,1)))
!!$     ALLOCATE(U_F(SIZE(H_FLOQUET,1),SIZE(H_FLOQUET,1)))
!!$     E_FLOQUET = 0.0  
!!$     CALL LAPACK_FULLEIGENVALUES(H_FLOQUET,SIZE(H_FLOQUET,1),E_FLOQUET,INFO)
!!$     U_F = H_FLOQUET ! FOURIER DECOMPOSITION OF THE STATES DRESSED BY MODE NUMBER 
!!$     DEALLOCATE(H_FLOQUET)
!!$!     write(*,*) e_floquet
!!$     !--- EVALUATE THE AVERAGE TRANSITION PROBATILIBIES IN THE BARE BASIS
!!$     P_AVG = 0.0
!!$     CALL MULTIMODETRANSITIONAVG(SIZE(U_F,1),size(MODES_NUM,1),FIELDS,MODES_NUM,U_F,E_FLOQUET,ID%D_BARE,P_AVG,INFO)   
!!$     WRITE(4,*) FIELDS(2)%omega,P_AVG
!!$         
!!$     !--- EVALUATE TIME-EVOLUTION OPERATOR IN THE BARE BASIS
!!$     T1 = 0.0
!!$     DO r=1,1!28
!!$        T2 = r*32.0*4.0*atan(1.0)/128
!!$        CALL MULTIMODETIMEEVOLUTINOPERATOR(SIZE(U_F,1),SIZE(MODES_NUM,1),MODES_NUM,U_F,E_FLOQUET,ID%D_BARE,FIELDS,T1,T2,U_AUX,INFO) 
!!$        P_AVG = ABS(U_AUX)**2
!!$        WRITE(3,*) t2,FIELDS(2)%OMEGA,ABS(U_AUX)**2
!!$     END DO
!!$     WRITE(3,*)
!!$     DEALLOCATE(E_FLOQUET)
!!$     DEALLOCATE(U_F)
!!$  END DO
!!$  
  
END PROGRAM MULTIMODEFLOQUET
