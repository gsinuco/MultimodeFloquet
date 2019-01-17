PROGRAM QUDITSTABILITY

  USE ATOMIC_PROPERTIES
  USE TYPES
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(8) :: E_DRESSED,magneticmoment
  DOUBLE PRECISION, DIMENSION(8,8) :: QUBITFREQUENCIES
  DOUBLE PRECISION  :: B_DC,B_plus,omega_plus,B_minus,omega_minus,B_ratio,RMS
  INTEGER           :: INFO
  CHARACTER*10200 printingstring 
  character*20 charaux
  TYPE(ATOM) :: ID
  INTEGER r,N_,i_,j_,l_
  DOUBLE PRECISION :: g_1,g_2

  OPEN(UNIT=3,file="Rb87_RFx2_QuditSensitivityB.dat", action="write")

  n_=64
  J    = 0.5
  I    = 1.5
  g_J  = 2.0
  g_I  = -0.000995

  F    = 2
  g_2 = g_J*(F*(F+1) - I*(I+1) +J*(J+1))/(2*F*(F+1)) + g_I*(F*(F+1) + I*(I+1) - J*(J+1))/(2*F*(F+1))!

  F    = 1
  g_1 = g_J*(F*(F+1) - I*(I+1) +J*(J+1))/(2*F*(F+1)) + g_I*(F*(F+1) + I*(I+1) - J*(J+1))/(2*F*(F+1))!

  INFO = 0
  CALL FLOQUETINIT('87Rb','B',0.2D1,ID,INFO)

  DO l_ =1,N_


     E_DRESSED = 0

     B_ratio = 1.0!2.0 - l_*2.0/N_

!     B_DC        = 0.5E-4!A
!     B_PLUS      = 0.05E-4!A
     B_DC        = 3.2E-4!B
     B_PLUS      = 0.2E-4!B
     B_MINUS     = B_ratio*B_PLUS
     OMEGA_PLUS  = mu_B*g_2*B_DC/hbar
     OMEGA_MINUS = (0.97 + (l_-1)*0.06/N_)*mu_B*abs(g_1)*B_DC/hbar

     !  printingstring = "# plot ""datos5.dat""  u 1 w l"
     !  printingstring = trim(printingstring)
     !  DO r=2,340
     !     printingstring = printingstring // """ u ", r , " r " 
     !     write(charaux,'(A,I3,A)') ", """" u ", r, " w l "
     !     write(*,*) charaux
     !     printingstring = trim(printingstring) // trim(charaux)
     !  END DO
     !  write(*,*) trim(printingstring)

     !write(*,*) '#',omega_plus,omega_minus,g_2,g_1
     DO r=1,N_
        !B_DC = 0.49E-4 + (r-1)*0.02E-4/N_!A 
        B_DC = 3.18E-4 + (r-1)*0.04E-4/N_!B
        !CALL RB87RFMWSPECTRUM(ID,B_DC,B_plus,omega_plus,B_minus,omega_minus,E_DRESSED,INFO)
        !WRITE(3,*) 1E4*B_DC,E_DRESSED
        CALL RB87RFMW_MAGNETICMOMENT(ID,B_DC,B_plus,omega_plus,B_minus,omega_minus,MAGNETICMOMENT,INFO)
        QUBITFREQUENCIES = 0.0
        RMS = 0.0
        DO i_ = 1,3
           DO j_ = 4,8
              QUBITFREQUENCIES(i_,j_) = abs(MAGNETICMOMENT(i_)-MAGNETICMOMENT(j_)) - 2.0
           END DO
           RMS = RMS + DOT_PRODUCT(QUBITFREQUENCIES(i_,:), QUBITFREQUENCIES(i_,:))
        END DO
        RMS = SQRT(RMS/15)
        WRITE(3,*) omega_plus/omega_minus,1E4*B_DC,RMS
     END DO
     WRITE(3,*)
  END DO

END PROGRAM QUDITSTABILITY


SUBROUTINE RB87RFMW_MAGNETICMOMENT(ID,B_DC,B_plus,omega_plus,B_minus,omega_minus,MAGNETICMOMENT,INFO)
  
  USE ATOMIC_PROPERTIES
  USE TYPES
  USE SUBINTERFACE
  USE SUBINTERFACE_LAPACK
  USE FLOQUETINIT_ 
  USE ARRAYS 
  
  
  IMPLICIT NONE
  TYPE(ATOM),                     INTENT(IN)    :: ID
  DOUBLE PRECISION, DIMENSION(8), INTENT(OUT)   :: MAGNETICMOMENT
  DOUBLE PRECISION,               INTENT(IN)    :: B_DC,B_plus,B_minus,omega_plus,omega_minus
  INTEGER,                        INTENT(INOUT) :: INFO
  
  DOUBLE PRECISION,parameter :: DB_ = 1.0E-12 ! 1 pT
  DOUBLE PRECISION, DIMENSION(8) :: E_QUDIT_DB,E_QUDIT
  
  
  E_qudit    = 0.0
  E_qudit_db = 0.0
  CALL RB87RFMWSPECTRUM(ID,B_DC,   B_plus,omega_plus,B_minus,omega_minus,E_QUDIT,   INFO)
  CALL RB87RFMWSPECTRUM(ID,B_DC+db_,B_plus,omega_plus,B_minus,omega_minus,E_QUDIT_DB,INFO)

  MagneticMoment = 1e-6*(A*(E_QUDIT_DB - E_QUDIT)/(2*pi*hbar))/(DB_) ! In MHz/T

!  WRITE(*,*) E_QUDIT_DB
!  WRITE(*,*) E_QUDIT
!  WRITE(*,*) 1E-6*MagneticMoment
  
END SUBROUTINE RB87RFMW_MAGNETICMOMENT 

SUBROUTINE RB87RFMWSPECTRUM(ID,B_DC,B_plus,omega_plus,B_minus,omega_minus,E_QUDIT,INFO)

  USE ATOMIC_PROPERTIES
  USE TYPES
  USE SUBINTERFACE
  USE SUBINTERFACE_LAPACK
  USE FLOQUETINIT_ 
  USE ARRAYS 


  IMPLICIT NONE
  TYPE(ATOM), INTENT(IN) ::               ID
  DOUBLE PRECISION, DIMENSION(8), INTENT(OUT) :: E_QUDIT
  DOUBLE PRECISION, INTENT(IN)    :: B_DC,B_plus,omega_plus,B_minus,omega_minus
  INTEGER,          INTENT(INOUT) :: INFO




  TYPE(MODE),       DIMENSION(:),   ALLOCATABLE :: FIELDS
!  TYPE(ATOM)                                       ID
  INTEGER,          DIMENSION(:),   ALLOCATABLE :: MODES_NUM
  INTEGER                                          TOTAL_FREQUENCIES,D_MULTIFLOQUET
  INTEGER                                          m,INDEX_UP,INDEX_DOWN,r,D_BARE
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: ENERGY,E_FLOQUET
  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: H__,U_F,U_AUX,U_B2D,U_F1,U_F2,U_F1_red,U_F2_red
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: P_AVG
  DOUBLE PRECISION                              :: T1,T2
 
  ! ===================================================
  !PARAMETERS REQUIRED TO DEFINE THE DRESSED BASIS
  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: U_FD           ! TRANSFORMATION FROM THE BARE TO THE DRESSED BASIS
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE   :: E_DRESSED      ! DRESSED SPECTRUM
  INTEGER                                       :: DRESSINGFIELDS, DRESSINGFLOQUETDIMENSION! NUMBER OF DRESSING FIELDS
  INTEGER,          DIMENSION(:), ALLOCATABLE   :: DRESSINGFIELDS_INDICES ! IDENTITY OF THE DRESSING FIELDS

  INTEGER                              :: TOTAL_FREQUENCIES_,NM_,FIELD_INDEX
  TYPE(MODE), DIMENSION(:),ALLOCATABLE :: FIELDS_
  INTEGER,    DIMENSION(:),ALLOCATABLE :: MODES_NUM_
  ! ===================================================


  INTEGER :: t,o,OFFSET_C,OFFSET_R,MANIFOLD_INDEX,index_
  DOUBLE PRECISION :: E_DOWN,E_UP
  DOUBLE PRECISION :: freq,pop

  INFO = 0
  !ALLOCATE(ENERGY(SIZE(J_Z,1)))
  !ALLOCATE(P_AVG(SIZE(J_z,1),SIZE(J_z,1)))
  !ALLOCATE(U_AUX(SIZE(J_z,1),SIZE(J_z,1)))
  
  D_BARE = ID%D_BARE
  ALLOCATE(MODES_NUM(3))

  MODES_NUM(1) = 1 !(STATIC FIELD)
  MODES_NUM(2) = 1 !(DRIVING BY SIGMA PLUS)
  MODES_NUM(3) = 1 !(DRIVING BY SIGMA MINUS)
  
  TOTAL_FREQUENCIES = SUM(MODES_NUM,1)
  ALLOCATE(FIELDS(TOTAL_FREQUENCIES))
  DO m=1,TOTAL_FREQUENCIES
     ALLOCATE(FIELDS(m)%V(ID%D_BARE,ID%D_BARE))
  END DO
  

  FIELDS(1)%X         = 0.0
  FIELDS(1)%Y         = 0.0
  FIELDS(1)%Z         = B_dc
  FIELDS(1)%phi_x     = 0.0
  FIELDS(1)%phi_y     = 0.0
  FIELDS(1)%phi_z     = 0.0
  FIELDS(1)%omega     = 0.0
  FIELDS(1)%N_Floquet = 0

  FIELDS(2)%X         = B_plus
  FIELDS(2)%Y         = B_plus
  FIELDS(2)%Z         = 0.0
  FIELDS(2)%phi_x     = 0.0
  FIELDS(2)%phi_y     = pi/2.0
  FIELDS(2)%phi_z     = 0.0
  FIELDS(2)%omega     = omega_plus
  FIELDS(2)%N_Floquet = 3
  
  FIELDS(3)%X         = B_minus
  FIELDS(3)%Y         = B_minus
  FIELDS(3)%Z         = 0.0
  FIELDS(3)%phi_x     = 0.0
  FIELDS(3)%phi_y     = -pi/2.0
  FIELDS(3)%phi_z     = 0.0
  FIELDS(3)%omega     = omega_minus
  FIELDS(3)%N_Floquet = 1
  


  DO m=1,TOTAL_FREQUENCIES
     FIELDS(m)%OMEGA = HBAR*FIELDS(m)%OMEGA/A        
  END DO

  D_MULTIFLOQUET = ID%D_BARE
  DO r=1,TOTAL_FREQUENCIES
     D_MULTIFLOQUET = D_MULTIFLOQUET*(2*FIELDS(r)%N_Floquet+1)
  END DO


  !=================================================================================
  !== MULTIMODE REPRESENTATION OF THE MW COUPLING IN THE BARE BASIS
  !=================================================================================
  !1. Build the components of the Hamiltonian
  CALL SETHAMILTONIANCOMPONENTS(ID,size(modes_num,1),total_frequencies,MODES_NUM,FIELDS,INFO)
  
  !2. Build the full multimode Floquet matrix
  CALL MULTIMODEFLOQUETMATRIX(ID,size(modes_num,1),total_frequencies,MODES_NUM,FIELDS,INFO)

  !3. DIAGONALISE THE MULTIMODE FLOQUET
  ALLOCATE(E_FLOQUET(SIZE(H_FLOQUET,1)))
  ALLOCATE(U_F(SIZE(H_FLOQUET,1),SIZE(H_FLOQUET,1)))
  E_FLOQUET = 0.0  
  CALL LAPACK_FULLEIGENVALUES(H_FLOQUET,SIZE(H_FLOQUET,1),E_FLOQUET,INFO)
  U_F = H_FLOQUET ! FOURIER DECOMPOSITION OF THE STATES DRESSED BY MODE NUMBER 
  DEALLOCATE(H_FLOQUET)
!!$!  write(*,*) e_floquet
  !4. SELECT THE EIGENSTATES   
  MANIFOLD_INDEX = ((2*FIELDS(3)%N_FLOQUET+1)*(2*FIELDS(2)%N_FLOQUET+1)*3 + 1)/2 
  DO r = MANIFOLD_INDEX-4,MANIFOLD_INDEX+4,4!(2*FIELDS(2)%N_FLOQUET+1)*3
!     write(*,*) r,E_FLOQUET(r)
     index_ = (r-manifold_index+4)/4 +1
     E_QUDIT(index_) = E_FLOQUET(r)
  END DO
  !WRITE(*,*)
  MANIFOLD_INDEX = (2*FIELDS(3)%N_FLOQUET+1)*(2*FIELDS(2)%N_FLOQUET+1)*3  &
       & + ((2*FIELDS(3)%N_FLOQUET+1)*(2*FIELDS(2)%N_FLOQUET+1)*5+1)/2 ! + 5*(FIELDS(2)%N_FLOQUET)
  DO r = MANIFOLD_INDEX-6,MANIFOLD_INDEX+6,3!(2*FIELDS(2)%N_FLOQUET+1)*3
     index_ =  (r-manifold_index+6)/3+1 + 3
     E_QUDIT(index_) = E_floquet(r)
!     write(*,*) r,E_FLOQUET(r)
  END DO

  

  DEALLOCATE(E_FLOQUET)
  DEALLOCATE(U_F)



!!$  !=================================================================================
!!$  !==== DEFINITION OF THE DRESSING FIELDS AND DRESSED BASIS AND VARIABLES NEEDED TO DEFINE THE MICROMOTION OPERATOR
!!$  !=================================================================================
!!$
!!$  DRESSINGFIELDS = 2   ! NUMBER OF DRESSING FIELDS
!!$  ALLOCATE(DRESSINGFIELDS_INDICES(DRESSINGFIELDS)) ! ARRAY THAT TELL US WHICH OF THE FIELD DEFINED ABOVE ARE THE DRESSING ONES
!!$  DRESSINGFIELDS_INDICES(1) = 1 
!!$  DRESSINGFIELDS_INDICES(2) = 2
!!$  ALLOCATE(H__(SIZE(E_DRESSED,1),SIZE(E_DRESSED,1)))
!!$
!!$  DO r=1,256
!!$     FIELDS(1)%Z = 3.38E-4 + (r-1)*0.04E-4/256
!!$     CALL MICROMOTIONFOURIERDRESSEDBASIS(ID,DRESSINGFIELDS_INDICES,MODES_NUM,FIELDS,U_FD,E_DRESSED,INFO)
!!$     !  write(*,*) e_dressed, size(E_dressed,1)
!!$     
!!$     
!!$     !=================================================================================
!!$     !== MULTIMODE REPRESENTATION OF THE MW COUPLING IN THE BARE BASIS
!!$     !=================================================================================
!!$     !1. Build the components of the Hamiltonian
!!$     CALL SETHAMILTONIANCOMPONENTS(ID,size(modes_num,1),total_frequencies,MODES_NUM,FIELDS,INFO)
!!$     
!!$     !2. Build the full multimode Floquet matrix
!!$     CALL MULTIMODEFLOQUETMATRIX(ID,size(modes_num,1),total_frequencies,MODES_NUM,FIELDS,INFO)
!!$     
!!$     !3. Select an off-diagonal section of the Floquet matrix.
!!$     OFFSET_C =  FIELDS(3)%N_FLOQUET    * (SIZE(E_DRESSED,1))
!!$     OFFSET_R = (FIELDS(3)%N_FLOQUET+1) * (SIZE(E_DRESSED,1))
!!$     H__      = H_FLOQUET(1+OFFSET_R:SIZE(E_DRESSED,1)+OFFSET_R,1+OFFSET_C:SIZE(E_DRESSED,1)+OFFSET_C)
!!$     
!!$     !4. Transform the couplings to the basis of dressed states
!!$     H__ = MATMUL(TRANSPOSE(CONJG(U_FD)),MATMUL(H__,U_FD))
!!$     
!!$     !5. Select pairs of dressed states
!!$     t = 2
!!$     o = 0
!!$     INDEX_DOWN = (2*Fdown+1)*FIELDS(2)%N_FLOQUET + t
!!$     E_DOWN = E_DRESSED(INDEX_DOWN)
!!$     DO t=1,5
!!$        INDEX_UP = (2*Fdown+1)*(2*FIELDS(2)%N_FLOQUET+1) + (2*Fup+1)*FIELDS(2)%N_FLOQUET + t + o*(2*Fup+1)
!!$        E_UP   = E_DRESSED(INDEX_UP)
!!$        write(3,*) abs(1E4*FIELDS(1)%Z),E_DOWN,E_UP,A*(E_UP-E_DOWN - 2.0)/(2*pi*hbar)!,&
!!$        !   & A*ABS(H__(INDEX_DOWN,INDEX_UP))/(2*pi*hbar),A*ABS(H__(INDEX_UP,INDEX_DOWN))/(2*pi*hbar)
!!$     END DO
!!$     DEALLOCATE(U_FD,E_DRESSED,H_FLOQUET)
!!$  !CALL WRITE_MATRIX(REAL(H__))
!!$  END DO
END SUBROUTINE RB87RFMWSPECTRUM

