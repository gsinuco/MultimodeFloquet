MODULE TYPES_C

  TYPE :: MODE_C
     DOUBLE PRECISION :: OMEGA
     COMPLEX*16       :: X,Y,Z
     DOUBLE PRECISION :: phi_x,phi_y,phi_z
     INTEGER          :: N_Floquet
     !COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: V
     !COMPLEX*16, DIMENSION(:),   ALLOCATABLE :: VALUES
     !INTEGER,    DIMENSION(:),   ALLOCATABLE :: ROW,COLUMN
  END TYPE MODE_C
  
  TYPE :: ATOM_C
     INTEGER          :: id_system
     INTEGER          :: D_BARE
     !DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: E_BARE
  END TYPE ATOM_C
  
END MODULE TYPES_C

MODULE MODES_4F
  ! THIS IS A GLOBAL DEFINITIO OF ATOM AND COUPLING PARAMETERS
  ! THEIR VALUES ARE INITIALISES BELOW, IN COUPLINGINIT_C USING AS INPUT
  ! THE VARIABLES PASSED BY C
  USE TYPES
  USE ISO_C_BINDING
  TYPE(MODE),DIMENSION(:),ALLOCATABLE :: COUPLING
  TYPE(ATOM)                          :: ATOM_
  INTEGER(C_INT),BIND(C,name="h_floquet_size") :: H_FLOQUET_SIZE
  LOGICAL COUPLINGALLOCATED 
  
  
  ! THESE ARE GLOBAL DEFINITIONS OF ARRAYS NEEDED TO REPRESENT THE MULTIMODE
  ! HAMILTONIAN REQUIRED BY THE MKL LIBRARY FOR SPARSE MATRICES.  
  COMPLEX*16, DIMENSION(:), ALLOCATABLE   :: VALUES__
  INTEGER,    DIMENSION(:), ALLOCATABLE   :: COLUMN__
  INTEGER,    DIMENSION(:), ALLOCATABLE   :: ROW_INDEX__

END MODULE MODES_4F


SUBROUTINE COUPLINGINIT_C(DB,NF,ATOM__C,COUPLING_C,INFO) 
  ! CUOPLING_C AND ATOM_C ARE C STRUCTURES

  USE TYPES_C
  USE MODES_4F
  IMPLICIT NONE
  
  INTEGER,                    INTENT(IN)    :: NF,DB
  TYPE(ATOM_C),               INTENT(IN)    :: ATOM__C
  TYPE(MODE_C),DIMENSION(NF), INTENT(IN)    :: COUPLING_C
  INTEGER,                    INTENT(INOUT) :: INFO
  
  INTEGER r

  IF(COUPLINGALLOCATED .EQV. .FALSE. ) THEN
     ALLOCATE(COUPLING(NF))
     COUPLINGALLOCATED = .TRUE.
     DO r=1,NF
       ! ALLOCATE(COUPLING(r)%VALUES(DB*DB))
       ! ALLOCATE(COUPLING(r)%COLUMN(DB*DB))
       ! ALLOCATE(COUPLING(r)%ROW(DB*DB))
        ALLOCATE(COUPLING(r)%V(DB,DB))
     END DO
     ALLOCATE(ATOM_%E_BARE(DB))
  END IF
  
  DO r=1,NF
     COUPLING(r)%OMEGA  = COUPLING_C(r)%OMEGA
     !COUPLING(r)%VALUES = 0.0
     !COUPLING(r)%ROW    = 0
     !COUPLING(r)%COLUMN = 0
     COUPLING(r)%X      = COUPLING_C(r)%X
     COUPLING(r)%Y      = COUPLING_C(r)%Y
     COUPLING(r)%Z      = COUPLING_C(r)%Z
     COUPLING(r)%phi_x  = COUPLING_C(r)%phi_x
     COUPLING(r)%phi_y  = COUPLING_C(r)%phi_y
     COUPLING(r)%phi_z  = COUPLING_C(r)%phi_z
     COUPLING(r)%N_Floquet = COUPLING_C(r)%N_Floquet
  END DO
  
  ATOM_%id_system = ATOM__C%id_system
  ATOM_%D_BARE    = ATOM__C%D_BARE
  
  INFO =0 

END SUBROUTINE COUPLINGINIT_C

