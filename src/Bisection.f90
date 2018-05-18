SUBROUTINE BISECTION(x_L,x_R,x,Y,INFO)

!  USE funciones
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(INOUT) :: x_L,x_R
  DOUBLE PRECISION, INTENT(OUT)   :: x,Y
  INTEGER,          INTENT(INOUT) :: INFO

  !DOUBLE PRECISION, EXTERNAL :: FUN

  DOUBLE PRECISION :: Y_L,Y_R,x_L_,x_R_
  DOUBLE PRECISION :: x_aux,y_aux,check,pi
  INTEGER DONE,I_,q


  DONE = 0
!  pi = 4.0*atan(1.0)
!  write(*,*) '#',x_L,Y_L,x_R,Y_R,"me", info
!  WRITE(*,*) "X_L",info
  CALL FUN_X(X_L,Y_L,INFO)
  CALL FUN_X(X_R,Y_R,INFO)
!  Y_L = FUN(X_L,INFO)
!  WRITE(*,*) "X_R",info
!  Y_R = FUN(x_R,INFO)
  x_L_ = x_L
  X_R_ = x_R
!  write(*,*) '#',x_L,Y_L,x_R,Y_R,"me", info

  check = Y_L*Y_R
  I_ = 0.0
  
  IF(check.GT.0)  THEN

     !     WRITE(*,*) "#THE INTERVAL DOES NOT CONTAIN A ZERO"
     !     write(*,*) '#',x_L,Y_L,x_R,Y_R,"me", info
     INFO = -1
     RETURN
     
  ELSE
     
     DO WHILE (DONE.NE.1)
        
        x_aux     = 0.5*(x_R_ + x_L_)
        CALL FUN_X(X_aux,Y_aux,INFO)
        !Y_aux     = FUN(x_aux,INFO)
        check     = Y_L*Y_aux
        IF(check.LT.0) THEN
           x_R_  = x_aux
        ELSE
           x_L_  = x_aux
        END IF
        I_ = I_ + 1
        !write(*,*) I_, x_aux,Y_aux,ABS(x_R_ - x_L_),ABS(x_R_ - x_L_)/x_aux
        !IF(I_.EQ.128 .OR. ABS(x_R_ - x_L_).LT.1.0E-12) THEN
           IF(I_.EQ.32 .OR. ABS(Y_aux).LT.1.0E-4) THEN
           x_L_ = x_aux
           x_R_ = x_aux
           DONE =1
           info = 0
        END IF
        
     END DO
     INFO = 0
  END IF
  
  x =  0.5*(x_r_ + x_l_)
  X_R = X_R_
  X_L = X_L_
  !Y = FUN(X,INFO)
  CALL FUN_X(X,Y,INFO)
  
END SUBROUTINE BISECTION

SUBROUTINE BISECTION_ROUTINE4FUN(y_dim,           N_DIM,k,FmF,x_L,x_R,x,Y,INFO)
  !y_dim : dimension of Y
  !N_DIM: dimension of the vectors X_L,x_R
  !k    : perform a bisection along direcition k
  !x    : root along direction k
  !Y    : value of the function along, it must be close to zero

  USE funciones
  IMPLICIT NONE

  INTEGER,                            INTENT(IN)    :: y_dim,N_DIM,K,FmF
  DOUBLE PRECISION, DIMENSION(N_DIM), INTENT(INOUT) :: X_L,X_R
  DOUBLE PRECISION,                   INTENT(OUT)   :: Y,X
  INTEGER,                            INTENT(INOUT) :: INFO

  DOUBLE PRECISION, DIMENSION(N_DIM) :: X_L_,X_R_,X_AUX
  DOUBLE PRECISION :: Y_L,Y_R
  DOUBLE PRECISION :: y_aux,check,pi,nan
  INTEGER DONE,I_,q

  DOUBLE PRECISION :: rho,rho_L,rho_R,sin_theta,cos_theta,sin_phi,cos_phi,z_0,theta,phi
  DOUBLE PRECISION, DIMENSION(3) :: POSITION
  !DOUBLE PRECISION, DIMENSION(y_dim) :: ENERGY,SHIFT_FULL


  DONE = 0
  pi = 4.0*atan(1.0)
  nan  = 0.
  nan  = nan / nan  

  !write(*,*) x_L,x_R,info
  CALL FUN_SUBROUTINE(y_dim,n_dim,k,FmF,X_L,Y_L,INFO)
  CALL FUN_SUBROUTINE(y_dim,n_dim,k,FmF,X_R,Y_R,INFO)
  x_L_ = x_L
  X_R_ = x_R
  !write(*,*) x_L,x_R,Y_L,Y_R
  check = Y_L*Y_R
  I_ = 0.0
  IF(check.GT.0)  THEN
   !  WRITE(*,*) "#THE INTERVAL DOES NOT CONTAIN A ZERO"
   !  write(*,*) '#',x_L(k),Y_L,x_R(k),Y_R,"me", info,y_dim,n_dim,k,FmF
     INFO = -1
     x=-1.0!nan
     RETURN
  ELSE     

     DO WHILE (DONE.NE.1)
         
        x_aux     = 0.5*(x_R_ + x_L_)
        CALL FUN_SUBROUTINE(y_dim,n_dim,k,FmF,X_AUX,Y_AUX,INFO)
        check     = Y_L*Y_aux
        IF(check.LT.0) THEN 
           !write(*,*) '#',x_L,Y_L,x_R,Y_R,"me", info
          x_R_  = x_aux          
        ELSE
          x_L_  = x_aux
        END IF
        I_ = I_ + 1 
!        write(*,*) I_, x_aux(1),Y_aux!,Y_aux,ABS(x_R_(k) - x_L_(k)),ABS(x_R_(k) - x_L_(k))/x_aux(k),Y_aux
        IF(I_.EQ.128 .OR. ABS(x_R_(K) - x_L_(k)).LT.1.0E-14) THEN
        !IF(I_.EQ.64 .OR. ABS(Y_aux).LT.1.0E-6) THEN
           x_L_ = x_aux
           x_R_ = x_aux
           DONE =1               
        END IF
        
     END DO
 !    INFO = -2
  END IF

  x =  0.5*(x_r_(K) + x_l_(K))
  X_R = X_R_
  X_L = X_L_

  X_R(K) = X
!  write(*,*) info,x
  CALL FUN_SUBROUTINE(y_dim,n_dim,k,FmF,X_R,Y,INFO)
!  write(*,*) info


END SUBROUTINE BISECTION_ROUTINE4FUN

SUBROUTINE FUN_SUBROUTINE(y_dim,n_dim,k,FmF,X_AUX,Y_AUX,INFO)
 
  USE physical_constants
  USE ATOMIC_PROPERTIES
  USE FLOQUET
  USE subinterface_lapack 
  USE SUBINTERFACE
  
  USE TYPES
  USE ARRAYS 
  
  USE QUADRUPOLE_TOP_TRAP
  USE ANTENNAS

  IMPLICIT NONE
    
  INTEGER,                            INTENT(IN)    :: y_dim,N_DIM,K,FmF
  DOUBLE PRECISION, DIMENSION(N_DIM), INTENT(IN)    :: X_AUX
  DOUBLE PRECISION,                   INTENT(OUT)   :: Y_AUX
  INTEGER,                            INTENT(INOUT) :: INFO

  INTEGER,          DIMENSION(5)              :: STEPPING
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: SHIFT_FULL,ENERGY,ENERGY_R,ENERGY_L
  TYPE(MODE) FIELD(MODES_NUM)
  TYPE(ATOM) Rb87
  
  DOUBLE PRECISION :: nan,DELTARHO,PHI,RHO,THETA
  INTEGER :: m

! FIND THE POINT WHERE THE CURVATURE OF THE POTENTIAL LANDSCAPE IS MINIMAL

    IF(Fmf.GT.y_dim) THEN
    WRITE(*,*) "SEARCHING A MINIMUN IN A NON-EXISTENT DRESSED ENERGY"
    RETURN
    END IF

  nan  = 0.
  nan  = nan / nan  
  Rb87%D_bare = TOTAL_STATES_LSI
  DELTARHO = 1.0E-6/32!128
  STEPPING = 1
  !write(*,*) "funtobisect:",y_dim,n_dim,k,FmF,X_AUX,Y_AUX,INFO
  
  ALLOCATE(SHIFT_FULL(y_dim))
  ALLOCATE(ENERGY(y_dim))
  ALLOCATE(ENERGY_R(y_dim))
  ALLOCATE(ENERGY_L(y_dim))
  
  shift_full = 0.0
  energy     = 0.0
  energy_r   = 0.0
  energy_l   = 0.0
  
  !  WRITE(*,*) MODES_NUM,TOTAL_STATES_LSI
  
  DO m=1,MODES_NUM
     ALLOCATE(FIELD(m)%V(Rb87%D_bare,Rb87%D_bare))
     FIELD(m)%V = DCMPLX(0.0,0.0)
  END DO
  
  
  STEPPING(3) = 1
  STEPPING(2) = 1
  STEPPING(1) = 1

  rho   = X_AUX(1)
  theta = X_AUX(2)
  phi   = X_AUX(3)
!  do m=1,128
!  deltarho = 20E-6 + m*110E-6/128.0
!  ENERGY = 0.0
!  SHIFT_FULL  = 0.0
!  POSITION(1) = 0.0!deltarho*sin(theta)*cos(phi)
!  POSITION(2) = 0.0!deltarho*sin(theta)*sin(phi)
!  POSITION(3) = (deltarho)*cos(theta) + z_0
!  CALL RFandMW_POTENTIALLANDSCAPE(POSITION,ENERGY,SHIFT_FULL,INFO)
!  write(*,*) POSITION,ENERGY(3)
!  end do
! !write(*,*)
!  DELTARHO = 1.0E-6/32.0!/32!64
  ENERGY = 0.0
  SHIFT_FULL  = 0.0
  POSITION(1) = rho*sin(theta)*cos(phi)
  POSITION(2) = rho*sin(theta)*sin(phi)
  POSITION(3) = rho*cos(theta) + z_0
!  write(*,*) POSITION,rho,z_0,theta,phi
!  write(*,*) info
  CALL RFandMW_POTENTIALLANDSCAPE(y_dim,POSITION,ENERGY,SHIFT_FULL,INFO)
!  INFO = 3
!  CALL SETPARAMETERS_XYZ(position,FIELD,INFO)
!  CALL COUPLINGMATRICES(FIELD,INFO)
!  CALL MWShiftOfRFDressed(TOTAL_STATES_LSI,FIELD,ENERGY,SHIFT_FULL,INFO)
  ENERGY_R = ENERGY+SHIFT_FULL
!  write(*,*) POSITION,ENERGY_R(8)

  ENERGY = 0.0
  SHIFT_FULL  = 0.0
  POSITION(1) = (rho+DELTARHO)*sin(theta)*cos(phi)
  POSITION(2) = (rho+DELTARHO)*sin(theta)*sin(phi)
  POSITION(3) = (rho+DELTARHO)*cos(theta)+z_0
  !write(*,*) POSITION
  CALL RFandMW_POTENTIALLANDSCAPE(y_dim,POSITION,ENERGY,SHIFT_FULL,INFO)
!  INFO = 3
!  CALL SETPARAMETERS_XYZ(position,FIELD,INFO)
!  CALL COUPLINGMATRICES(FIELD,INFO)
!  CALL MWShiftOfRFDressed(TOTAL_STATES_LSI,FIELD,ENERGY,SHIFT_FULL,INFO)
  ENERGY_L = ENERGY+SHIFT_FULL
!  write(*,*) POSITION,ENERGY_L(8)

  Y_AUX =  -(ENERGY_R(FmF) - ENERGY_L(FmF))!/abs(ENERGY_R(3)+ENERGY_L(3))

!  write(*,*) rho,Y_aux,FmF

END SUBROUTINE FUN_SUBROUTINE

SUBROUTINE GRADIENTENERGYLANDSCAPE(y_dim,n_dim,k,FmF,X_AUX,Y_AUX,INFO)
 
  USE physical_constants
  USE ATOMIC_PROPERTIES
  USE FLOQUET
  USE subinterface_lapack 
  USE SUBINTERFACE
  
  USE TYPES
  USE ARRAYS 
  
  USE QUADRUPOLE_TOP_TRAP
  USE ANTENNAS

  IMPLICIT NONE
  
  INTEGER,                            INTENT(IN)    :: y_dim,N_DIM,K,FmF
  DOUBLE PRECISION, DIMENSION(N_DIM), INTENT(IN)    :: X_AUX
  DOUBLE PRECISION,                   INTENT(OUT)   :: Y_AUX
  INTEGER,                            INTENT(INOUT) :: INFO
  
  INTEGER,          DIMENSION(5)              :: STEPPING
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: SHIFT_FULL,ENERGY,ENERGY_R,ENERGY_L
  TYPE(MODE) FIELD(MODES_NUM)
  TYPE(ATOM) Rb87
  
  DOUBLE PRECISION :: nan,DELTARHO,PHI,RHO,THETA
  INTEGER :: m
  
  ! FIND THE POINT WHERE THE CURVATURE OF THE POTENTIAL LANDSCAPE IS MINIMAL

  IF(Fmf.GT.y_dim) THEN
     WRITE(*,*) "SEARCHING A MINIMUN IN A NON-EXISTENT DRESSED ENERGY"
     RETURN
  END IF
  
  nan  = 0.
  nan  = nan / nan  
  Rb87%D_bare = TOTAL_STATES_LSI
  DELTARHO = 1.0E-6/128
  STEPPING = 1
  
  ALLOCATE(SHIFT_FULL(y_dim))
  ALLOCATE(ENERGY(y_dim))
  ALLOCATE(ENERGY_R(y_dim))
  ALLOCATE(ENERGY_L(y_dim))
  
  shift_full = 0.0
  energy     = 0.0
  energy_r   = 0.0
  energy_l   = 0.0
  
  DO m=1,MODES_NUM
     ALLOCATE(FIELD(m)%V(Rb87%D_bare,Rb87%D_bare))
     FIELD(m)%V = DCMPLX(0.0,0.0)
  END DO
  
  
  STEPPING(3) = 1
  STEPPING(2) = 1
  STEPPING(1) = 1
  
  rho   = X_AUX(1)
  theta = X_AUX(2)
  phi   = X_AUX(3)
  
  ENERGY = 0.0
  SHIFT_FULL  = 0.0
  POSITION(1) = rho*sin(theta)*cos(phi)
  POSITION(2) = rho*sin(theta)*sin(phi)
  POSITION(3) = rho*cos(theta) + z_0
  CALL RFandMW_POTENTIALLANDSCAPE(y_dim,POSITION,ENERGY,SHIFT_FULL,INFO)
  ENERGY_R = ENERGY+SHIFT_FULL
  
  ENERGY = 0.0
  SHIFT_FULL  = 0.0
  POSITION(1) = (rho+DELTARHO)*sin(theta)*cos(phi)
  POSITION(2) = (rho+DELTARHO)*sin(theta)*sin(phi)
  POSITION(3) = (rho+DELTARHO)*cos(theta)+z_0
  CALL RFandMW_POTENTIALLANDSCAPE(y_dim,POSITION,ENERGY,SHIFT_FULL,INFO)
  ENERGY_L = ENERGY+SHIFT_FULL
  
  Y_AUX =  -(ENERGY_R(FmF) - ENERGY_L(FmF))!/abs(ENERGY_R(3)+ENERGY_L(3))
  
END SUBROUTINE GRADIENTENERGYLANDSCAPE


!SUBROUTINE BISECTION_ROUTINE4FUN_(fun_subroutine_,N_DIM,k,FmF,x_L,x_R,x,Y,INFO)
SUBROUTINE BISECTION_ROUTINE4FUN_(fun_subroutine_,y_dim,N_DIM,k,FmF,x_L,x_R,x,Y,INFO)

  !N_DIM: dimension of the vectors X_L,x_R
  !k    : perform a bisection along direcition k
  !x    : root along direction k
  !Y    : value of the function along, it must be close to zero

!  interface
!    subroutine FUN_SUBROUTINE2(n_dim,k,X_L,Y_L,INFO)
!     INTEGER,          INTENT(IN)    :: N_DIM,K
!     DOUBLE PRECISION, INTENT(IN),DIMENSION(N_DIM):: X_AUX
!     DOUBLE PRECISION, INTENT(OUT)    :: Y_AUX
!     INTEGER,          INTENT(INOUT) :: INFO
!    end subroutine
!  end interface
!  USE funciones
!  USE FUNCTIONFORBISECTION
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: Y_DIM,N_DIM,K,FmF
  DOUBLE PRECISION, INTENT(INOUT), DIMENSION(N_DIM) :: X_L,X_R
  DOUBLE PRECISION, INTENT(OUT)   :: Y,X
  INTEGER,          INTENT(INOUT) :: INFO

  DOUBLE PRECISION, DIMENSION(N_DIM) :: X_L_,X_R_,X_AUX
  DOUBLE PRECISION :: Y_L,Y_R
  DOUBLE PRECISION :: y_aux,check,pi,nan
  INTEGER DONE,I_,q

  DOUBLE PRECISION :: rho,rho_L,rho_R,sin_theta,cos_theta,sin_phi,cos_phi,z_0,theta,phi
  DOUBLE PRECISION, DIMENSION(3) :: POSITION
  DOUBLE PRECISION, DIMENSION(8) :: ENERGY,SHIFT_FULL

  DONE = 0
  pi = 4.0*atan(1.0)
  nan  = 0.
  nan  = nan / nan
!  write(*,*) y_dim,N_DIM,k,FmF,x_L,x_R,x,Y,INFO,pi,nan

  CALL FUN_SUBROUTINE_(y_dim,n_dim,k,FmF,X_L,Y_L,INFO)
  CALL FUN_SUBROUTINE_(y_dim,n_dim,k,FmF,X_R,Y_R,INFO)
  x_L_ = x_L
  X_R_ = x_R

  check = Y_L*Y_R
  I_ = 0.0
  IF(check.GT.0)  THEN
!     WRITE(*,*) "#THE INTERVAL DOES NOT CONTAIN A ZERO"
!     write(*,*) '#',x_L(k),Y_L,x_R(k),Y_R,"me", info
         INFO = -1
     x=-1.0!nan
     RETURN
  ELSE

     DO WHILE (DONE.NE.1)

        x_aux     = 0.5*(x_R_ + x_L_)
        CALL FUN_SUBROUTINE_(y_dim,n_dim,k,FmF,X_AUX,Y_AUX,INFO)
        !CALL FUN_SUBROUTINE_(n_dim,k,FmF,X_AUX,Y_AUX,INFO)
        check     = Y_L*Y_aux
        IF(check.LT.0) THEN
         !  write(*,*) '#',x_L,Y_L,x_R,Y_R,"me", info
          x_R_  = x_aux
        ELSE
          x_L_  = x_aux
        END IF
        I_ = I_ + 1
        !write(*,*) I_, x_aux,Y_aux!,Y_aux,ABS(x_R_(k) - x_L_(k)),ABS(x_R_(k) - x_L_(k))/x_aux(k),Y_aux
        IF(I_.EQ.128 .OR. ABS(x_R_(K) - x_L_(k)).LT.1.0E-14) THEN
        !IF(I_.EQ.64 .OR. ABS(Y_aux).LT.1.0E-6) THEN
           x_L_ = x_aux
           x_R_ = x_aux
           DONE =1
        END IF

     END DO
 !    INFO = -2
  END IF

  x =  0.5*(x_r_(K) + x_l_(K))
  X_R = X_R_
  X_L = X_L_

  X_R(K) = X
  !CALL FUN_SUBROUTINE_(n_dim,k,FmF,X_R,Y,INFO)
  CALL FUN_SUBROUTINE_(y_dim,n_dim,k,FmF,X_R,Y,INFO)

END SUBROUTINE BISECTION_ROUTINE4FUN_

SUBROUTINE FUN_SUBROUTINE__(n_dim,k,FmF,X_AUX,Y_AUX,INFO)

  USE physical_constants
  USE ATOMIC_PROPERTIES
  USE FLOQUET
  USE subinterface_lapack
  USE SUBINTERFACE

  USE TYPES
  USE ARRAYS

  USE QUADRUPOLE_TOP_TRAP
  USE ANTENNAS
  USE MAGNETIC_BOTTLE

  IMPLICIT NONE

  INTEGER,          INTENT(IN)    :: N_DIM,FmF,K
  DOUBLE PRECISION, INTENT(IN),DIMENSION(N_DIM):: X_AUX
  DOUBLE PRECISION, INTENT(OUT)    :: Y_AUX
  INTEGER,          INTENT(INOUT) :: INFO

  INTEGER,          DIMENSION(5)                :: STEPPING
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: SHIFT_FULL,ENERGY
  TYPE(MODE) FIELD(MODES_NUM)
  TYPE(ATOM) Rb87

  DOUBLE PRECISION :: nan,DELTARHO,RHO
  INTEGER :: m,t,n,N_

  DOUBLE PRECISION, DIMENSION(3):: X_L,X_R
  DOUBLE PRECISION :: theta,phi,cos_theta,sin_theta,cos_phi,sin_phi,rho_L,rho_R,v,ENERGY_R,ENERGY_L

! FIND THE POINT WHERE THE CURVATURE OF THE POTENTIAL LANDSCAPE IS MINIMAL

  nan  = 0.
  nan  = nan / nan
  Rb87%D_bare = TOTAL_STATES_LSI
  DELTARHO = 1.0E-6/2
  STEPPING = 1

  ALLOCATE(SHIFT_FULL(TOTAL_STATES_LSI))
  ALLOCATE(ENERGY(TOTAL_STATES_LSI))

  shift_full = 0.0
  energy     = 0.0

  !  WRITE(*,*) MODES_NUM,TOTAL_STATES_LSI

  DO m=1,MODES_NUM
     ALLOCATE(FIELD(m)%V(Rb87%D_bare,Rb87%D_bare))
     FIELD(m)%V = DCMPLX(0.0,0.0)
  END DO


  STEPPING(3) = 1
  STEPPING(2) = 1
  STEPPING(1) = 1

  OMEGA_MW = X_AUX(k)
  t=1
  n=1
  N_ = 64
  DO M=1,2!2 ! find the two potential wells

     phi   = 0.0
     theta = (M-1)*4.0*atan(1.0)

     sin_phi   = sin(phi)
     cos_phi   = cos(phi)
     sin_theta = sin(theta)
     cos_theta = cos(theta)

     rho_L = 10.0e-6
     rho_R = 200.0e-6

     x_L(1) = rho_L
     x_L(2) = theta
     x_L(3) = phi

     x_R(1) = rho_R
     x_R(2) = theta
     x_R(3) = phi

     CALL BISECTION_ROUTINE4FUN(24,3,1,FmF,x_L,x_R,rho,v,INFO)
     !write(*,*) 1E-6*(omega_mw-2*A/hbar)/(2*pi),rho,phi,theta,v

     POSITION(1) = rho*sin_theta*cos_phi
     POSITION(2) = rho*sin_theta*sin_phi
     POSITION(3) = rho*cos_theta + z_0
    ! write(*,*) (1E-6/(2*pi))*(OMEGA_MW-2*A/hbar),M,position
     CALL RFandMW_POTENTIALLANDSCAPE(24,POSITION,ENERGY,SHIFT_FULL,INFO)
     IF(M.EQ.1) ENERGY_L = ENERGY(FmF) + shift_full(FmF)
     IF(M.EQ.2) ENERGY_R = ENERGY(FmF) + shift_full(FmF)
     !WRITE(*,*) M,1e-6*(OMEGA_MW-2*A/HBAR)/(2*PI),ENERGY(8)+SHIFT_FULL(8)
   END DO

   Y_AUX = ENERGY_L - ENERGY_R
   !write(*,*) 1e-6*(omega_mw-2*A/hbar)/(2*pi),Y_AUX

END SUBROUTINE FUN_SUBROUTINE__

