!!$SUBROUTINE SETPARAMETERS_XYZ(POSITION,FIELD,INFO)       
!!$  
!!$  USE physical_constants
!!$  USE ATOMIC_PROPERTIES
!!$!  USE QUADRUPOLE_TOP_TRAP
!!$  USE FLOQUET
!!$  USE TYPES
!!$!  USE ANTENNAS
!!$!  USE MAGNETIC_BOTTLE
!!$  IMPLICIT NONE
!!$  DOUBLE PRECISION, DIMENSION(3), INTENT(IN)    :: POSITION
!!$  INTEGER,                        INTENT(INOUT) :: INFO  
!!$  TYPE(MODE),DIMENSION(MODES_NUM),INTENT(OUT)   :: FIELD
!!$
!!$  DOUBLE PRECISION, DIMENSION(3) :: POSITION_CART,POSITION_POLAR
!!$  DOUBLE PRECISION, DIMENSION(3) :: B_MW_CART,B_MW_POLAR
!!$  DOUBLE PRECISION, DIMENSION(3) :: B_RF_CART,B_RF_POLAR
!!$  DOUBLE PRECISION, DIMENSION(3) :: B_DC_CART,B_DC_POLAR
!!$  DOUBLE PRECISION :: x,y,z
!!$  DOUBLE PRECISION :: rho,theta,tau,phi,OmegaRF_Rabi!,OMEGA_TOP
!!$  DOUBLE PRECISION :: B_lin,B_plus,B_minus,Br,Btheta
!!$  INTEGER          :: m,n,o,p,N_,T
!!$
!!$
!!$
!!$  BAND = 0
!!$
!!$  !---------- Giromagnetic factors in the coupled representation F ---------------
!!$
!!$  J = L+S
!!$  F    = 2
!!$  gF_2 = g_J*(F*(F+1) - I*(I+1) +J*(J+1))/(2*F*(F+1)) + g_I*(F*(F+1) + I*(I+1) - J*(J+1))/(2*F*(F+1))!
!!$  F    = 1
!!$  gF_1 = g_J*(F*(F+1) - I*(I+1) +J*(J+1))/(2*F*(F+1)) + g_I*(F*(F+1) + I*(I+1) - J*(J+1))/(2*F*(F+1))!
!!$  G_F  = (g_J-g_I)/16.0
!!$
!!$
!!$
!!$!--------- MULTIMODE COUPLINGS
!!$!---------- THREE DIMENSIONAL FIELD CONFIGURATION
!!$
!!$
!!$
!!$!  IF(INFO.EQ.3) THEN
!!$     POSITION_CART = POSITION
!!$
!!$     x = POSITION(1)
!!$     y = POSITION(2)
!!$     z = POSITION(3)
!!$!     write(*,*) position
!!$     CALL CARTESSIANTOPOLAR(POSITION_CART,POSITION_POLAR)
!!$     if(abs(POSITION_POLAR(2)).lt.1E-6) POSITION_POLAR(2) = 0.0
!!$     RHO   = POSITION_POLAR(1)
!!$     theta = 2.0*atan(1.0) - POSITION_POLAR(2)
!!$     phi   = POSITION_POLAR(3)
!!$          
!!$
!!$     CALL RINGFIELD(POSITION_POLAR(1),POSITION_POLAR(2),Br,Btheta,R_MW,M_MW)
!!$!     write(*,*) "XYZ",x,y,z,rho,theta,phi,Br,Btheta
!!$     !write(*,*) POSITION, Br,Btheta
!!$     
!!$     B_MW_CART(1) = Br*cos(theta)*cos(phi) - Btheta*sin(theta)*cos(phi)
!!$     B_MW_CART(2) = Br*cos(theta)*sin(phi) - Btheta*sin(theta)*sin(phi)
!!$     B_MW_CART(3) = Br*sin(theta) + Btheta*cos(theta)
!!$     
!!$     !B_MW_CART = 10.0*B_MW_CART
!!$!     write(*,*) B_MW_CART
!!$     CALL RINGFIELD(POSITION_POLAR(1),POSITION_POLAR(2),Br,Btheta,R_RF,M_RF)
!!$     !write(*,*) rho,theta,phi,Br,Btheta
!!$     
!!$     B_RF_CART(1) = Br*cos(theta)*cos(phi) - Btheta*sin(theta)*cos(phi)
!!$     B_RF_CART(2) = Br*cos(theta)*sin(phi) - Btheta*sin(theta)*sin(phi)
!!$     B_RF_CART(3) = Br*sin(theta) + Btheta*cos(theta)
!!$ !    write(*,*) B_RF_CART
!!$          
!!$     !B_DC_CART(1) = alpha_TOP*x 
!!$     !B_DC_CART(2) = alpha_TOP*y 
!!$     !B_DC_CART(3) = -2*alpha_TOP*(z-z_0)
!!$
!!$     B_DC_CART(1) =  alpha_TOP*x + 3.0*A3*y*x
!!$     B_DC_CART(2) =        -A1 - 3.0*A3*y*y + 1.5*A3*(x*x+(z-z_0)*(z-z_0))
!!$     B_DC_CART(3) = -alpha_TOP*(z-z_0) + 3.0*A3*y*(z-z_0)
!!$     !write(*,*) rho,theta,phi,REAL(B_RF_CART)
!!$     !write(*,*) position,mu_B*0.5*sqrt(abs(B_DC_CART(1))**2+abs(B_DC_CART(2))**2+abs(B_DC_CART(3))**2)/(hbar),omega_RF
!!$  !   write(*,*) B_DC_CART
!!$
!!$     FIELD(1)%OMEGA =  0.0
!!$     FIELD(1)%Bx    =  B_DC_CART(1)  
!!$     FIELD(1)%phi_x =  0.0
!!$     FIELD(1)%By    =  B_DC_CART(2)
!!$     FIELD(1)%phi_y =  0.0
!!$     FIELD(1)%Bz    =  B_DC_CART(3)
!!$     FIELD(1)%phi_z =  0.0
!!$     
!!$!     omega_RF       =   mu_B*0.5*sqrt((alpha_TOP*50E-6)**2+(A1+1.5*A3*(50.0E-6)**2)**2)/hbar
!!$
!!$!     write(*,*) position,z_0,alpha_top,A1,A3,&
!!$!        & 1E-6*omega_RF/(2*pi),1E-6*abs(mu_B*0.5*sqrt(abs(FIELD(1)%Bx)**2+abs(FIELD(1)%By)**2+abs(FIELD(1)%Bz)**2)/hbar)/(2*pi)
!!$     theta          =   0.0
!!$     FIELD(2)%OMEGA =   omega_RF!0.5*alpha_TOP*50.0E-6*mu_B/hbar
!!$     FIELD(2)%Bx    =   B_RF_CART(1)
!!$     FIELD(2)%phi_x =   0.0
!!$     FIELD(2)%By    =   B_RF_CART(2)
!!$     FIELD(2)%phi_y =   0.0
!!$     FIELD(2)%Bz    =   B_RF_CART(3)
!!$     FIELD(2)%phi_z =   0.0
!!$
!!$!     write(*,*) position,z_0,alpha_top,A1,A3,&
!!$!        & 1E-6*omega_RF/(2*pi),1E-6*real(mu_B*0.5*sqrt(FIELD(2)%Bx**2+FIELD(2)%By**2+FIELD(2)%Bz**2)/hbar)/(2*pi)
!!$
!!$     phi            =   0.0
!!$     tau            =   0.0
!!$     FIELD(3)%OMEGA =   omega_MW!2.0*A/hbar! + 2*pi*(4.0E6 + 0.0*(t-1.0)*0.7E6/N_)
!!$     FIELD(3)%Bx    =   B_MW_CART(1)
!!$     FIELD(3)%phi_x =   0.0
!!$     FIELD(3)%By    =   B_MW_CART(2)
!!$     FIELD(3)%phi_y =   0.0
!!$     FIELD(3)%Bz    =   B_MW_CART(3)
!!$     FIELD(3)%phi_z =   0.0
!!$!     write(*,*) position,z_0,alpha_top,A1,A3,&
!!$!        & 1E-6*omega_RF/(2*pi),1E-6*real(mu_B*0.5*sqrt(FIELD(3)%Bx**2+FIELD(3)%By**2+FIELD(3)%Bz**2)/hbar)/(2*pi)
!!$
!!$!     write(*,*) 'set parameters',omega_RF,omega_MW
!!$!    write(*,*) field(1)%omega,field(2)%omega, field(3)%omega
!!$
!!$ ! END IF
!!$
!!$!     write(*,*) MODES_NUM,total_states_LSI
!!$  
!!$  DO m=1,MODES_NUM
!!$
!!$     ALLOCATE(FIELD(m)%V(TOTAL_STATES_LSI,TOTAL_STATES_LSI))
!!$     FIELD(m)%V = 0.0
!!$     FIELD(m)%Bx = FIELD(m)%Bx*exp(DCMPLX(0.0,1.0)*FIELD(m)%phi_x)
!!$     FIELD(m)%By = FIELD(m)%By*exp(DCMPLX(0.0,1.0)*FIELD(m)%phi_y)
!!$     FIELD(m)%Bz = FIELD(m)%Bz*exp(DCMPLX(0.0,1.0)*FIELD(m)%phi_z)
!!$!     WRITE(*,*) SIZE(FIELD(M)%V,1)
!!$  END DO
!!$
!!$
!!$ 
!!$END SUBROUTINE SETPARAMETERS_XYZ
!!$


! Initialize the field parameters as complex numbers

SUBROUTINE SETPARAMETERS(STEPPING,FIELD,INFO)       

! STEPPING : INPUT, INTEGER VECTOR USED TO STEP OVER PARAMETERS
! FIELD    : OUTPUT, FIELD-TYPE, INITIALIZEDF FIELD TYPE
! INFO     : INOUT, INTEGER, ERROR FLAGS
  
  USE physical_constants
  USE ATOMIC_PROPERTIES
!  USE QUADRUPOLE_TOP_TRAP
  USE FLOQUET
  USE TYPES


  
  IMPLICIT NONE
  INTEGER,DIMENSION(:),INTENT(IN)             :: STEPPING
  INTEGER,             INTENT(INOUT)          :: INFO  
  TYPE(MODE),DIMENSION(MODES_NUM),INTENT(OUT) :: FIELD

  DOUBLE PRECISION :: theta,tau,phi,OmegaRF_Rabi!,OMEGA_TOP
  DOUBLE PRECISION :: B_lin,B_plus,B_minus
  INTEGER          :: m,n,N_,T

  BAND = 0

  !---------- Giromagnetic factors in the coupled representation F ---------------

  J = L+S
  F    = 2
  gF_2 = g_J*(F*(F+1) - I*(I+1) +J*(J+1))/(2*F*(F+1)) + g_I*(F*(F+1) + I*(I+1) - J*(J+1))/(2*F*(F+1))!
  F    = 1
  gF_1 = g_J*(F*(F+1) - I*(I+1) +J*(J+1))/(2*F*(F+1)) + g_I*(F*(F+1) + I*(I+1) - J*(J+1))/(2*F*(F+1))!
  G_F  = (g_J-g_I)/16.0


  n  = STEPPING(1)
  t  = STEPPING(2)
  N_ = STEPPING(SIZE(STEPPING,1))


!--------- ONE QUBITS

  !OMEGA_TOP      = 2*pi*5.02E3

  FIELD(1)%OMEGA =  0.0
  FIELD(1)%X    =  0.0
  FIELD(1)%phi_x =  0.0
  FIELD(1)%Y    =  0.0
  FIELD(1)%phi_y =  0.0
  FIELD(1)%Z    =  1.0
  FIELD(1)%phi_z =  0.0


  theta          =   0.0
  FIELD(2)%OMEGA =   1.0
  FIELD(2)%X    =   0.125*2.0
  FIELD(2)%phi_x =   0.0*pi/2.0
  FIELD(2)%Y    =   0.0
  FIELD(2)%phi_y =   0.0
  FIELD(2)%Z    =   0.0
  FIELD(2)%phi_z =   0.0


  theta          =   0.0
  FIELD(3)%OMEGA =   0.0
  FIELD(3)%X    =   0.0!(0.125/2.0)/10.0
  FIELD(3)%phi_x =   0.0
  FIELD(3)%Y    =   0.0!(0.125/2.0)/10.0
  FIELD(3)%phi_y =   0.0
  FIELD(3)%Z    =   0.0 ! TO
  FIELD(3)%phi_z =   0.0
  IF(stepping(2).EQ.0) THEN

     FIELD(3)%OMEGA =   2.0*(2*N_FLOQUET(2)+3)*FIELD(2)%OMEGA
     FIELD(3)%X    =   0.0!(0.125/2.0)/10.0
     FIELD(3)%phi_x =   0.0
     FIELD(3)%Y    =   0.0!(0.125/2.0)/10.0
     FIELD(3)%phi_y =   0.0
     FIELD(3)%Z    =   0.0 ! TO
     FIELD(3)%phi_z =   0.0
  ELSE
     FIELD(3)%OMEGA =   FIELD(2)%OMEGA - 2.0*0.125 + 2.0*(stepping(3)-1)*2.0*0.125/N_!0.01 + (t-1)*0.125/N_!0.125/2.0 - 0.02 + 0.1*(t-1)/N_! + 0.125*(t-1.0)/N_ !0.01 + (t-1.0)*0.125/N_!0.05 + 0.03*(t-1.0)/N_!
     FIELD(3)%X    =   0.0!(0.125/2.0)/10.0
     FIELD(3)%phi_x =   0.0
     FIELD(3)%Y    =   0.125/4.0 !(0.125/2.0)/10.0
     FIELD(3)%phi_y =   0.0
     FIELD(3)%Z    =   0.0!0.125/4.0
     FIELD(3)%phi_z =   0.0*pi/3.0
  END IF
  


!!$!--------- two QUBITS
!!$
!!$  !OMEGA_TOP      = 2*pi*5.02E3
!!$
!!$  FIELD(1)%OMEGA =  0.0
!!$  FIELD(1)%X    =  0.0
!!$  FIELD(1)%phi_x =  0.0
!!$  FIELD(1)%Y    =  0.0
!!$  FIELD(1)%phi_y =  0.0
!!$  FIELD(1)%Z    =  1.0
!!$  FIELD(1)%phi_z =  0.0
!!$
!!$
!!$  theta          =   0.0
!!$  FIELD(2)%OMEGA =   1.0
!!$  FIELD(2)%X    =   0.1225
!!$  FIELD(2)%phi_x =   0.0
!!$  FIELD(2)%Y    =   0.0
!!$  FIELD(2)%phi_y =   0.0
!!$  FIELD(2)%Z    =   0.0
!!$  FIELD(2)%phi_z =   0.0
!!$  
!!$  phi            =   0.0
!!$  tau            =   0.0
!!$  FIELD(3)%OMEGA =   1.0*(t-1.0)*sqrt((1.0 - FIELD(2)%OMEGA)**2+ FIELD(2)%X**2)/N_
!!$  FIELD(3)%X    =   0.0
!!$  FIELD(3)%phi_x =   0.0
!!$  FIELD(3)%Y    =   0.0
!!$  FIELD(3)%phi_y =   0.0
!!$  FIELD(3)%Z    =   0.125/2.0
!!$  FIELD(3)%phi_z =   0.0

  DO m=1,MODES_NUM
     FIELD(m)%X = FIELD(m)%X*exp(DCMPLX(0.0,1.0)*FIELD(m)%phi_x)
     FIELD(m)%Y = FIELD(m)%Y*exp(DCMPLX(0.0,1.0)*FIELD(m)%phi_y)
     FIELD(m)%Z = FIELD(m)%Z*exp(DCMPLX(0.0,1.0)*FIELD(m)%phi_z)
  END DO

  
END SUBROUTINE SETPARAMETERS
