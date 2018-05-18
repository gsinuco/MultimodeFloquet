MODULE FLOQUET

  INTEGER,                 PARAMETER :: BASIS = 0 ! 1: use basis dressed by the first non-static mode 
                                                  ! 0: use the basis of Zeeman states
  INTEGER,                 PARAMETER :: MODES_NUM = 3
  INTEGER, DIMENSION(MODES_NUM)      :: N_Floquet = (/0,6,2/)!(/0,2,2/)!(/0,2,2,2/)  ! The first mode is static
  INTEGER                            :: BAND = 0    ! BAND STORAGE OF THE HAMILTONIAN 0: NO
                                                    !                                 1: YES  

  INTEGER,                 PARAMETER     :: MODES_NUM_DRESSING = 2 ! to set the dressed states, we need two modes: the static field and the dressing field
  INTEGER, DIMENSION(MODES_NUM_DRESSING) :: N_FLOQUET_DRESSING  = (/0,6/) 
END MODULE FLOQUET
