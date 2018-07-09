SUBROUTINE WRITE_MATRIX(A)
! it writes a matrix of doubles nxm on the screen
  DOUBLE PRECISION, DIMENSION(:,:) :: A
  CHARACTER(LEN=65) STRING
  CHARACTER(LEN=55) aux_char
  integer :: aux

  aux = int(UBOUND(A,2))
!  write(*,*) aux
  write(aux_char,"(I4)") aux
  aux_char = trim(aux_char)
  write(string,"(A1,I4,A6)") "(",aux,"E15.6)"

  DO I = LBOUND(A,1), UBOUND(A,1)
     WRITE(*,string) (A(I,J), J = LBOUND(A,2), UBOUND(A,2))
  END DO
  WRITE(*,*)
  WRITE(*,*)
END SUBROUTINE WRITE_MATRIX

SUBROUTINE WRITE_MATRIX_INT(A)
!it writes a matrix of integer nxm on the screen
  INTEGER, DIMENSION(:,:) :: A
  WRITE(*,*)
  DO I = LBOUND(A,1), UBOUND(A,1)
     WRITE(*,*) (A(I,J), J = LBOUND(A,2), UBOUND(A,2))
  END DO
END SUBROUTINE WRITE_MATRIX_INT


SUBROUTINE COORDINATEPACKING(D,A,V,R,C,index,INFO)
  IMPLICIT NONE
  INTEGER,INTENT(IN):: D
  COMPLEX*16,DIMENSION(D,D),INTENT(IN)  :: A
  COMPLEX*16,DIMENSION(D*D),INTENT(OUT) :: V
  INTEGER, DIMENSION(D*D),  INTENT(OUT) :: R,C
  INTEGER, INTENT(OUT)   :: index
  INTEGER, INTENT(INOUT) :: INFO
  
  INTEGER I,J
  V=0
  R=0
  C=0
  
  index = 1
  DO I=1,D
     DO J=1,D
        IF(ABS(A(I,J)).GT.0) THEN
           V(index) = A(I,J)
           R(index) = I
           C(index) = J
           index = index+1
        END IF
     END DO
  END DO
  index = index-1
END SUBROUTINE COORDINATEPACKING

MODULE MERGINGARRAYS
  INTERFACE
     SUBROUTINE APPENDARRAYS(V,B,INFO)
       COMPLEX*16, DIMENSION(:),ALLOCATABLE, INTENT(INOUT) :: V
       COMPLEX*16, DIMENSION(:),INTENT(IN)    :: B
       INTEGER,                 INTENT(INOUT) :: INFO
     END SUBROUTINE APPENDARRAYS
     SUBROUTINE APPENDARRAYSI(V,B,INFO)
       INTEGER, DIMENSION(:),ALLOCATABLE, INTENT(INOUT) :: V
       INTEGER, DIMENSION(:),INTENT(IN)    :: B
       INTEGER,                 INTENT(INOUT) :: INFO
     END SUBROUTINE APPENDARRAYSI
  END INTERFACE
END MODULE MERGINGARRAYS

SUBROUTINE APPENDARRAYS(V,B,INFO)
  COMPLEX*16, DIMENSION(:),ALLOCATABLE, INTENT(INOUT) :: V
  COMPLEX*16, DIMENSION(:),INTENT(IN)    :: B
  INTEGER,                 INTENT(INOUT) :: INFO
  
  COMPLEX*16,DIMENSION(:),ALLOCATABLE :: tmp_arr
!  write(*,*) V
!  write(*,*) B
  ALLOCATE(tmp_arr(SIZE(V,1)+SIZE(B,1)))
  tmp_arr(1:SIZE(V,1))=V
  tmp_arr(SIZE(V,1)+1:SIZE(tmp_arr))=B
  DEALLOCATE(V)
  ALLOCATE(V(SIZE(tmp_arr)))
  V=tmp_arr
END SUBROUTINE APPENDARRAYS

SUBROUTINE APPENDARRAYSI(V,B,INFO)
  INTEGER, DIMENSION(:),ALLOCATABLE, INTENT(INOUT) :: V
  INTEGER, DIMENSION(:),INTENT(IN)    :: B
  INTEGER,                 INTENT(INOUT) :: INFO
  
  COMPLEX*16,DIMENSION(:),ALLOCATABLE :: tmp_arr
  
  ALLOCATE(tmp_arr(SIZE(V,1)+SIZE(B,1)))
  tmp_arr(1:SIZE(V,1))=V
  tmp_arr(SIZE(V,1)+1:SIZE(tmp_arr))=B
  DEALLOCATE(V)
  ALLOCATE(V(SIZE(tmp_arr)))
  V=tmp_arr
END SUBROUTINE APPENDARRAYSI


! An example of reading a simple control with labels.
!
! Jason Blevins <jrblevin@sdf.lonestar.org>
! Durham, May 6, 2008

! German sinuco <gsinuco@gmailc.om>
! Brighton, July 17, 2017

!program control_file

!!$subroutine control_file(filename,N,INFO)
!!$
!!$! reads the file containing the parameters 
!!$! of the experiment.
!!$
!!$! asign those parameters to field variables for the main calculation
!!$
!!$! Check there are enough parameters.
!!$
!!$  USE DCRFMWFIELDS
!!$
!!$  implicit none
!!$
!!$
!!$  INTEGER, INTENT(IN)     :: N
!!$  character(N),INTENT(INOUT) :: filename
!!$  INTEGER, INTENT(INOUT)  :: INFO
!!$
!!$  ! Input related variables
!!$  character(len=200) :: buffer, label
!!$  integer :: pos
!!$  integer, parameter :: fh = 15
!!$  integer :: ios = 0,ios2=0
!!$  integer :: line = 0
!!$  integer :: parameter_counter = 0
!!$
!!$  ! Control file variables
!!$  real :: pi
!!$  integer, dimension(5) :: vector
!!$
!!$  !DOUBLE PRECISION, DIMENSION(17):: FIELD_CONFIGURATION
!!$  ! DC, RF and MW parameters
!!$  !! DC FIELD
!!$  !  DOUBLE PRECISION B_DC_X, B_DC_Y, B_DC_Z   
!!$
!!$  !! RF FIELD
!!$  !  DOUBLE PRECISION  OMEGA_RF, PHI_RF_X,PHI_RF_Y,PHI_RF_Z
!!$  !  DOUBLE PRECISION  B_RF_X,B_RF_Y,B_RF_Z
!!$  
!!$  !! MW FIELD
!!$  !  DOUBLE PRECISION  OMEGA_MW, PHI_MW_X,PHI_MW_Y,PHI_MW_Z
!!$  !  DOUBLE PRECISION  B_MW_X,B_MW_Y,B_MW_Z
!!$  
!!$  filename = trim(filename)
!!$!  write(*,*) filename
!!$  open(fh, file=filename)
!!$  IF (ios.ne.0) THEN
!!$     WRITE(*,*) "ERROR OPEING THE FILE OF PARAMETERS"
!!$     RETURN
!!$  END IF
!!$  
!!$  ! ios is negative if an end of record condition is encountered or if
!!$  ! an endfile condition was detected.  It is positive if an error was
!!$  ! detected.  ios is zero otherwise.
!!$
!!$  do while (ios == 0)
!!$     read(fh,'(A)',iostat=ios) buffer!, iostat=ios) buffer
!!$!     write(*,*) ios
!!$     if (ios == 0) then
!!$        line = line + 1
!!$
!!$        ! Find the first instance of whitespace.  Split label and data.
!!$        pos = scan(buffer, '=')
!!$      !  write(*,*) 'POS ',pos
!!$        label = buffer(1:pos-1)
!!$        buffer = buffer(pos+1:)
!!$        label =  trim(label)
!!$        !write(*,*) "label ", label
!!$        select case (label)
!!$        case ('B_DC_X')
!!$           read(buffer, *, iostat=ios) B_DC_X
!!$           print *, '#Bdcx: ', B_DC_X
!!$           parameter_counter = parameter_counter + 1
!!$        case ('B_DC_Y')
!!$           read(buffer, *, iostat=ios) B_DC_Y
!!$           print *, '#Bdcy: ', B_DC_Y
!!$           parameter_counter = parameter_counter + 1
!!$        case ('B_DC_Z')
!!$           read(buffer, *, iostat=ios) B_DC_Z
!!$           print *, '#Bdcz: ', B_DC_Z
!!$           parameter_counter = parameter_counter + 1
!!$
!!$        case ('OMEGA_RF')
!!$           read(buffer, *, iostat=ios) omega_RF
!!$           print *, '#omegarf: ', omega_RF
!!$           parameter_counter = parameter_counter + 1
!!$        case ('B_RF_X')
!!$           read(buffer, *, iostat=ios) B_RF_X
!!$           print *, '#BRFX: ', B_RF_X
!!$           parameter_counter = parameter_counter + 1
!!$        case ('B_RF_Y')
!!$           read(buffer, *, iostat=ios) B_RF_Y
!!$           print *, '#BRFY: ', B_RF_Y
!!$           parameter_counter = parameter_counter + 1
!!$        case ('B_RF_Z')
!!$           read(buffer, *, iostat=ios) B_RF_Z
!!$           print *, '#BRFZ ', B_RF_Z
!!$           parameter_counter = parameter_counter + 1
!!$
!!$        case ('PHI_RF_X')
!!$           read(buffer, *, iostat=ios) PHI_RF_X
!!$           print *, '#PHIRFX: ', PHI_RF_X
!!$           parameter_counter = parameter_counter + 1
!!$        case ('PHI_RF_Y')
!!$           read(buffer, *, iostat=ios) PHI_RF_Y
!!$           print *, '#PHIRFY ', PHI_RF_Y
!!$           parameter_counter = parameter_counter + 1
!!$        case ('PHI_RF_Z')
!!$           read(buffer, *, iostat=ios) PHI_RF_Z
!!$           print *, '#PHIRFZ ', PHI_RF_Z
!!$           parameter_counter = parameter_counter + 1
!!$
!!$        case ('OMEGA_MW')
!!$           read(buffer, *, iostat=ios) omega_MW
!!$           print *, '#OMEGA MW: ', omega_MW
!!$           parameter_counter = parameter_counter + 1
!!$        case ('B_MW_X')
!!$           read(buffer, *, iostat=ios) B_MW_X
!!$           print *, '#BMWX: ', B_MW_X
!!$           parameter_counter = parameter_counter + 1
!!$        case ('B_MW_Y')
!!$           read(buffer, *, iostat=ios) B_MW_Y
!!$           print *, '#BMWY ', B_MW_Y
!!$           parameter_counter = parameter_counter + 1
!!$        case ('B_MW_Z')
!!$           read(buffer, *, iostat=ios) B_MW_Z
!!$           print *, '#BMWZ ', B_MW_Z
!!$           parameter_counter = parameter_counter + 1
!!$        case ('PHI_MW_X')
!!$           read(buffer, *, iostat=ios) PHI_MW_X
!!$           print *, '#PHIMWX: ', PHI_MW_X
!!$           parameter_counter = parameter_counter + 1
!!$        case ('PHI_MW_Y')
!!$           read(buffer, *, iostat=ios) PHI_MW_Y
!!$           print *, '#PHIMWY ', PHI_MW_Y
!!$           parameter_counter = parameter_counter + 1
!!$        case ('PHI_MW_Z')
!!$           read(buffer, *, iostat=ios) PHI_MW_Z
!!$           print *, '#PHIMWZ ', PHI_MW_Z
!!$           parameter_counter = parameter_counter + 1
!!$        case default
!!$!           print *, 'Skipping invalid label at line', line
!!$        end select
!!$     else
!!$        !        write(*,*) ios
!!$     end if
!!$  end do
!!$
!!$  INFO = 0.0
!!$  IF (parameter_counter.ne.17) then
!!$     write(*,*) 'error: not all required parameters have been defined'
!!$     INFO = -1
!!$     return
!!$  end if
!!$  
!!$END subroutine control_file
!!$
!!$
