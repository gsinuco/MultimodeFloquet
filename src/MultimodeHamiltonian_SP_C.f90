SUBROUTINE MULTIMODEFLOQUETMATRIX_SP_C(ATOM__C,NM,NF,MODES_NUM,FIELDS_C,INFO)!VALUES_,ROW_INDEX_,COLUMN_,SP,INFO)
! THIS SUBROUTINE BUILDS THE MULTIMODE FLOQUET MATRIX

!ATOM_      (IN)    : type of quantum system
!NM         (IN)    : number of modes
!NF         (IN)    : number of driving fields
!MODES_NUM  (IN)    : vector indicating the number of harmonics of each driving field (mode)
!FIELDS_C     (IN)    : Fields

! THE FOLLOWING VARIABLES ARE DECLARED AS GLOBAL ALLOCATABLE ARRAYS. THIS SUBROUTINE SET THEIR VALUES AND SIZE.

!VALUES_    (OUT)   : Hamiltonian values
!ROW_INDEX_ (OUT)   : vector indicating the row position of values
!COLUMN_    (OUT)   : vector indicating the column position of the values
!INFO       (INOUT) : error flag. INFO=0 means there is no error

  USE TYPES_C             !(modes.f90)
  USE MERGINGARRAYS     !(utils.f90)
  USE SPARSE_INTERFACE  !(sparse_utils.f90)
  USE MODES_4F          ! DEFINED IN modes_C.f90, declares atom_,coupling, values__, row_index__, column__
  
  IMPLICIT NONE
  INTEGER                  ,            INTENT(IN)    :: NM,NF
  TYPE(MODE_C), DIMENSION(NF),          INTENT(INout) :: FIELDS_C
  TYPE(ATOM_C),                         INTENT(IN)    :: ATOM__C
  INTEGER,    DIMENSION(NM),            INTENT(IN)    :: MODES_NUM
  INTEGER,                              INTENT(INOUT) :: INFO
  
!  COMPLEX*16, DIMENSION(:), ALLOCATABLE  :: VALUES__
!  INTEGER,    DIMENSION(:), ALLOCATABLE  :: COLUMN__
!  INTEGER,    DIMENSION(:), ALLOCATABLE  :: ROW_INDEX__
  !INTEGER,                              INTENT(OUT)   :: SP
 
!    write(*,*) 'FORTRAN FLOQUETMATRIX_SP SAYS',NM,NF,MODES_NUM!, COULPLIG(3)%OMEGA
  IF (INFO.EQ.0) THEN      
    CALL MULTIMODEFLOQUETMATRIX_SP(ATOM_,NM,NF,MODES_NUM,COUPLING,VALUES__,ROW_INDEX__,COLUMN__,INFO)
!    WRITE(*,*) "FORTRAN MULTIMODEHAMILTONAISPC SAYS: SIZE(VALUES__,1) =)",SIZE(VALUES__,1),SIZE(ROW_INDEX__,1)  
    H_FLOQUET_SIZE = SIZE(ROW_INDEX__,1)-1
  END IF

END SUBROUTINE MULTIMODEFLOQUETMATRIX_SP_C ! _SP  sparse packing

