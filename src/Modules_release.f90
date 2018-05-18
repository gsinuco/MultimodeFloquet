MODULE ATOMIC_PROPERTIES
  USE physical_constants
  IMPLICIT NONE
  DOUBLE PRECISION :: L=0.0,  S = 0.5
  DOUBLE PRECISION :: mass_at = 87*amu
  DOUBLE PRECISION :: I,g_I,g_J
  DOUBLE PRECISION :: J,F,gf,mf
  DOUBLE PRECISION :: gF_2,gF_1,G_F
  DOUBLE PRECISION :: A,a_s,alpha_E
  INTEGER          :: Fup,Fdown,Ftotal
  INTEGER          :: Total_states_LSI
  
!!$  !electro
!!$  DOUBLE PRECISION, PARAMETER :: Ftotal      =  0.5
!!$  DOUBLE PRECISION, PARAMETER :: I           =  0.0
!!$  DOUBLE PRECISION, PARAMETER :: g_J         =  2.0
!!$  DOUBLE PRECISION, PARAMETER :: g_I         = -0.000995
!!$  DOUBLE PRECISION, PARAMETER :: A           =  2*pi*hbar*3.417341E9
!!$  DOUBLE PRECISION, PARAMETER :: a_s         =  5.29E-9 !5.77E-9
!!$  DOUBLE PRECISION, PARAMETER :: alpha_E     =  2*pi*hbar*0.0794*1E-4
!!$  DOUBLE PRECISION, parameter :: Temperature =  2.0E-6
!!$  INTEGER,          PARAMETER :: Fup         =  0
!!$  INTEGER,          PARAMETER :: Fdown       =  0
!!$  DOUBLE PRECISION            :: N_particles,gF_2,gF_1
!!$  INTEGER                     :: Total_states_LSI

  !87Rb
  DOUBLE PRECISION :: I_87Rb   =  1.5  
  DOUBLE PRECISION :: J_87Rb   =  0.5  
  DOUBLE PRECISION :: gJ_87Rb  =  2.0
  DOUBLE PRECISION :: gI_87Rb  = -0.000995
  DOUBLE PRECISION :: A_87Rb   =  2*pi*hbar*3.417341E9
  DOUBLE PRECISION :: a_s_87Rb = 5.77E-9
  DOUBLE PRECISION :: alpha_E_87Rb = 2*pi*hbar*0.0794*1E-4
  INTEGER          :: Fup_87Rb     =  2
  INTEGER          :: Fdown_87Rb   =  1

  !6Li
  DOUBLE PRECISION :: I_6Li   =  1.0  
  DOUBLE PRECISION :: J_6Li   =  0.5  
  DOUBLE PRECISION :: gJ_6Li  =  2.0
  DOUBLE PRECISION :: gI_6Li  = -0.000995
  DOUBLE PRECISION :: A_6Li   =  2*pi*hbar*152.137E6
  DOUBLE PRECISION :: a_s_6Li = 5.77E-9
  DOUBLE PRECISION :: alpha_E_6Li = 2*pi*hbar*0.0794*1E-4
  INTEGER          :: Fup_6Li     =  1
  INTEGER          :: Fdown_6Li   =  1

  !qubit
  DOUBLE PRECISION :: I_qubit   =  0.0
  DOUBLE PRECISION :: J_qubit   =  0.0  
  DOUBLE PRECISION :: gJ_qubit  =  1.0
  DOUBLE PRECISION :: gI_qubit  =  0.0
  DOUBLE PRECISION :: A_qubit   =  1.0
  DOUBLE PRECISION :: a_s_qubit =  0.0
  DOUBLE PRECISION :: alpha_E_qubit = 0.0
  INTEGER          :: Fup_qubit     =  1
  INTEGER          :: Fdown_qubit   =  1

  !spin
  DOUBLE PRECISION :: I_spin   =  0.0
  DOUBLE PRECISION :: J_spin   =  0.0  
  DOUBLE PRECISION :: gJ_spin  =  1.0
  DOUBLE PRECISION :: gI_spin  =  0.0
  DOUBLE PRECISION :: A_spin   =  1.0
  DOUBLE PRECISION :: a_s_spin =  0.0
  DOUBLE PRECISION :: alpha_E_spin = 0.0
  INTEGER          :: Fup_spin     =  1
  INTEGER          :: Fdown_spin   =  1

  !6Li  !DOUBLE PRECISION, PARAMETER :: I=1.0  
  !DOUBLE PRECISION, PARAMETER :: g_J = 2.0023010
  !DOUBLE PRECISION, PARAMETER :: g_I = -0.0004476540
  !DOUBLE PRECISION, PARAMETER :: A = 2*pi*hbar*152.137E6
  
END MODULE ATOMIC_PROPERTIES

