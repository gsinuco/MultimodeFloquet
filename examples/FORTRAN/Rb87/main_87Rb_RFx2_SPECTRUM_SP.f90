!export LD_LIBRARY_PATH="/opt/intel/compilers_and_libraries_2017/linux/mkl/lib/intel64"; 


PROGRAM QUDITSTABILITY

  USE ATOMIC_PROPERTIES
  USE TYPES
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(8) :: E_DRESSED
  DOUBLE PRECISION  :: B_DC,B_plus,omega_plus,B_minus,omega_minus
  INTEGER           :: INFO

  double precision :: g_1,g_2

  J    = 0.5
  I    = 1.5
  g_J  = 2.0
  g_I  = -0.000995

  F    = 2
  g_2 = g_J*(F*(F+1) - I*(I+1) +J*(J+1))/(2*F*(F+1)) + g_I*(F*(F+1) + I*(I+1) - J*(J+1))/(2*F*(F+1))!

  F    = 1
  g_1 = g_J*(F*(F+1) - I*(I+1) +J*(J+1))/(2*F*(F+1)) + g_I*(F*(F+1) + I*(I+1) - J*(J+1))/(2*F*(F+1))!

    
  INFO = 0
  E_DRESSED = 0
  
  B_DC        = 1.0E-4
  B_PLUS      = 0.1E-4
  B_MINUS     = 0.1E-4
  OMEGA_PLUS  = mu_B*g_2*B_DC/hbar
  OMEGA_MINUS = mu_B*abs(g_1)*B_DC/hbar
  
  write(*,*) omega_plus,omega_minus,g_2,g_1
  CALL RB87RFMWSPECTRUM(B_DC,B_plus,omega_plus,B_minus,omega_minus,E_DRESSED,INFO)


END PROGRAM QUDITSTABILITY

SUBROUTINE RB87RFMWSPECTRUM(B_DC,B_plus,omega_plus,B_minus,omega_minus,E_QUDIT,INFO)

  USE ATOMIC_PROPERTIES
  USE TYPES
  USE SPARSE_INTERFACE
  USE SUBINTERFACE
  USE SUBINTERFACE_LAPACK
  USE FLOQUETINIT_ 
  USE ARRAYS 


  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(8), INTENT(OUT) :: E_QUDIT
  DOUBLE PRECISION, INTENT(IN)    :: B_DC,B_plus,omega_plus,B_minus,omega_minus
  INTEGER,          INTENT(INOUT) :: INFO
  !PARAMETERS NEEDED TO DEFINE THE SPARSE MATRIX
  INTEGER,    DIMENSION(:), ALLOCATABLE :: ROW_INDEX,COLUMN
  COMPLEX*16, DIMENSION(:), ALLOCATABLE :: VALUES




  TYPE(MODE),       DIMENSION(:),   ALLOCATABLE :: FIELDS
  TYPE(ATOM)                                       ID
  INTEGER,          DIMENSION(:),   ALLOCATABLE :: MODES_NUM
  INTEGER                                          TOTAL_FREQUENCIES,D_MULTIFLOQUET
  INTEGER                                          m,INDEX_UP,INDEX_DOWN,r,D_BARE
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: ENERGY,E_FLOQUET
  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: H__,U_F,U_AUX,U_B2D,U_F1,U_F2,U_F1_red,U_F2_red
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: P_AVG
  DOUBLE PRECISION                              :: T1,T2,E_L,E_R
 
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


  INTEGER :: t,o,OFFSET_C,OFFSET_R
  DOUBLE PRECISION :: E_DOWN,E_UP
  DOUBLE PRECISION :: freq,pop

  !OPEN(UNIT=3,file="Rb87_RFMW_SPECTRUM.dat", action="write")


  INFO = 0
!  CALL FLOQUETINIT('spin','U',0.1D1,ID,INFO)
  CALL FLOQUETINIT('87Rb','B',0.2D1,ID,INFO)
  ALLOCATE(ENERGY(SIZE(J_Z,1)))
  ALLOCATE(P_AVG(SIZE(J_z,1),SIZE(J_z,1)))
  ALLOCATE(U_AUX(SIZE(J_z,1),SIZE(J_z,1)))
  
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
  FIELDS(2)%Y         = 0.0!B_plus
  FIELDS(2)%Z         = 0.0
  FIELDS(2)%phi_x     = 0.0
  FIELDS(2)%phi_y     = 0.0!-pi/2.0
  FIELDS(2)%phi_z     = 0.0
  FIELDS(2)%omega     = omega_plus
  FIELDS(2)%N_Floquet = 3

  FIELDS(3)%X         = 2.0*B_minus
  FIELDS(3)%Y         = 0.0!B_minus
  FIELDS(3)%Z         = 0.0
  FIELDS(3)%phi_x     = 0.0
  FIELDS(3)%phi_y     = 0.0!pi/2.0
  FIELDS(3)%phi_z     = 0.0
  FIELDS(3)%omega     = omega_minus
  FIELDS(3)%N_Floquet = 2
  

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
  DO t=1,3
     DO r=1,size(FIELDS(1)%V,1)
        DO m=1,size(FIELDS(1)%V,1)
           IF (abs(FIELDS(t)%V(r,m)).LT.10E-8) FIELDS(t)%V(r,m) = 0.0
        END DO
     END DO
  END DO
  !2. Build the full multimode Floquet matrix
  CALL MULTIMODEFLOQUETMATRIX_SP(ID,SIZE(MODES_NUM,1),total_frequencies,MODES_NUM,FIELDS,VALUES,ROW_INDEX,COLUMN,INFO)
!  write(*,*) real(Values)
!  DO m=1,size(VALUES,1)
!     IF (REAL(abs(VALUES(M))).LT.10E-8) VALUES(M) = 0.0
!  END DO
!  write(*,*) real(Values)
  !3. DIAGONALISE THE MULTIMODE FLOQUET
  E_L =  -6.0
  E_R =   6.0
  ALLOCATE(E_FLOQUET(D_MULTIFLOQUET))
  ALLOCATE(U_F(D_MULTIFLOQUET,D_MULTIFLOQUET))
  E_FLOQUET = 0.0
  U_F = 0.0
  INFO = 1
  CALL MKLSPARSE_FULLEIGENVALUES(D_MULTIFLOQUET,SIZE(VALUES,1),VALUES,ROW_INDEX,COLUMN,E_L,E_R,E_FLOQUET,U_F,INFO)
!  U_F = H_FLOQUET ! FOURIER DECOMPOSITION OF THE STATES DRESSED BY MODE NUMBER 
!  DEALLOCATE(H_FLOQUET)
!  write(*,*) E_FLOQUET

  !4. SELECT THE EIGENSTATES 
  E_QUDIT = 0.0


!  DEALLOCATE(E_FLOQUET)
  !  DEALLOCATE(U_F)
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




