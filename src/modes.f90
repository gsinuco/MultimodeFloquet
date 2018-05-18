MODULE FEAST
  integer     fpm(128)
  real*8      Emin,Emax
  real*8      epsout
  integer     loop
  integer     M0 ! initial guess 
  integer     M1 ! total number of eigenvalues found
  integer     info_FEAST
  real*8,     DIMENSION(:),   ALLOCATABLE :: E, RES ! vector of eigenvalues
  complex*16, DIMENSION(:,:), ALLOCATABLE :: X      ! matrix with eigenvectore
END MODULE FEAST

MODULE TYPES

  TYPE :: MODE
     DOUBLE PRECISION :: OMEGA
     COMPLEX*16       :: X,Y,Z
     DOUBLE PRECISION :: phi_x,phi_y,phi_z
     INTEGER          :: N_Floquet
     COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: V
     COMPLEX*16, DIMENSION(:),   ALLOCATABLE :: VALUES
     INTEGER,    DIMENSION(:),   ALLOCATABLE :: ROW,COLUMN
  END TYPE MODE
  
  TYPE :: ATOM
     INTEGER          :: id_system
     INTEGER          :: D_BARE
     DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: E_BARE
  END TYPE ATOM

  TYPE :: HARMONIC_FACTORS
!     INTEGER,   DIMENSION(:),   ALLOCATABLE :: H 
     COMPLEX*16,DIMENSION(:,:), ALLOCATABLE :: U,U_r,U_AVG
     INTEGER,   DIMENSION(:),   ALLOCATABLE :: n
  END type HARMONIC_FACTORS

  TYPE :: MWCOUPLING
     COMPLEX*16,DIMENSION(:,:),ALLOCATABLE :: TOP
     COMPLEX*16,DIMENSION(:,:),ALLOCATABLE :: TOP_DAGGER
     COMPLEX*16,DIMENSION(:,:),ALLOCATABLE :: DC
     COMPLEX*16,DIMENSION(:,:),ALLOCATABLE :: DC_DAGGER
     COMPLEX*16,DIMENSION(:,:),ALLOCATABLE :: MW
     COMPLEX*16,DIMENSION(:,:),ALLOCATABLE :: MW_DAGGER
     COMPLEX*16,DIMENSION(:,:),ALLOCATABLE :: RF
     COMPLEX*16,DIMENSION(:,:),ALLOCATABLE :: RF_DAGGER
  END type MWCOUPLING
END MODULE TYPES

!!$MODULE FLOQUET
!!$
!!$  INTEGER,                 PARAMETER :: BASIS = 0 ! 1: use basis dressed by the first non-static mode 
!!$                                                  ! 0: use the basis of Zeeman states
!!$  INTEGER,                 PARAMETER :: MODES_NUM = 3
!!$  INTEGER, DIMENSION(MODES_NUM)      :: N_Floquet = (/0,3,2/)!(/0,4,2/)!(/0,5,3/)!(/0,2,2/)!(/0,2,2,2/)  ! The first mode is static
!!$  INTEGER                            :: BAND = 0    ! BAND STORAGE2OF THE HAMILTONIAN 0: NO
!!$                                                    !                                 1: YES  
!!$
!!$  INTEGER,                 PARAMETER     :: MODES_NUM_DRESSING = 2 ! to set the dressed states, we need two modes: the static field and the dressing field
!!$  INTEGER, DIMENSION(MODES_NUM_DRESSING) :: N_FLOQUET_DRESSING = (/0,3/)
!!$END MODULE FLOQUET
