SUBROUTINE MULTIMODEFLOQUETMATRIX_C(ATOM__C,NM,NF,MODES_NUM,FIELD_C,INFO)
  !ID,size(modes_num,1),total_frequencies,MODES_NUM,FIELDS,INFO
  !  USE FLOQUET
  !ATOM_C type atom, -> dimension of the bare Hilbert space
  !NM -> number of modes
  !NF -> Number of Fields
  !MODES_NUM -> number of harmonics of each mode
  !FIELD_C -> Field couplings
  !INFO

  USE ARRAYS
  USE TYPES_C
  USE MODES_4F

  IMPLICIT NONE
  INTEGER,                     INTENT(IN)    :: NM,NF
  INTEGER,                     INTENT(INOUT) :: INFO
  INTEGER,      DIMENSION(NM), INTENT(IN)    :: MODES_NUM
  TYPE(MODE_C), DIMENSION(NF), INTENT(IN)    :: FIELD_C
  TYPE(ATOM_C),                INTENT(IN)    :: ATOM__C                     


  IF (INFO.EQ.0) THEN
     CALL MULTIMODEFLOQUETMATRIX(ATOM_,NM,NF,MODES_NUM,COUPLING,INFO)
     H_FLOQUET_SIZE = SIZE(H_FLOQUET,1)
  ELSE
     WRITE(*,*) "THERE IS AN ERROR UPSTREAM :INFO.NE.0", INFO
  END IF
END SUBROUTINE MULTIMODEFLOQUETMATRIX_C

