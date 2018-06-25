SUBROUTINE MULTIMODETIMEEVOLUTINOPERATOR(D,NM,MODES_NUM,U_F_MODES,E_MULTIFLOQUET,D_BARE,FIELD,T1,T2,U,INFO) 

  ! TIME EVOLUTION OPERATOR OF A MULTIMODE DRESSED SYSTEM. THE EVOLUTION OPERATOR IS WRITEN IN THE BASIS USED TO EXPRESS THE 
  ! MULTIMODE FLOQUET HAMILTONIAN
  ! U : MATRIX OF AMPLITUED OF PROBABILITIES FOR TRANSITIONS BETWEEN T1 TO T2
!!$  D              (IN)   : DIMENSION OF THE EXTENDED HILBERT SPACE (SIZE OF THE MULTIMODE FLOQUET MATRIX)
!!$  NM             (IN)   : NUMBER OF MODES            
!!$  MODES_NUM      (IN)   : VECTOR (NM) INDICATING THE NUMBER OF HARMONICS OF EACH MODE
!!$  U_F_MODES      (IN)   : TRANSFORMATION, DIMENSOON (D,D) 
!!$  E_MULTIFLOQUET (IN)   : MULTIMODE FLOQUET SPECTRUM
!!$  D_BARE         (IN)   : DIMENSION OF THE BARE HILBERT SPACE
!!$  FIELD          (IN)   : STRUCTURE DESCRIBING THE COUPLINGS
!!$  T1             (IN)   : INITIAL TIME
!!$  T2             (IN)   : FINAL TIME
!!$  U              (OUT)  : TRANFORMATION BETWEEN THE EXTENDED BARE BASIS AND THE FLOQUET STATES, DIMENSION (D_BARE,D)
!!$  INFO           (INOUT): (POSSIBLE) ERROR FLAG

  USE TYPES
  USE SUBINTERFACE_LAPACK


  IMPLICIT NONE
  INTEGER,                                    INTENT(IN)    :: D,D_BARE,NM ! DIMENSION OF THE MULTIMODE FLOQUET SPACE AND THE BARE BASIS
  INTEGER,                                    INTENT(INOUT) :: INFO
  INTEGER,          DIMENSION(NM),            INTENT(IN)    :: MODES_NUM
  TYPE(MODE),       DIMENSION(NM),            INTENT(IN)    :: FIELD  ! FIELDS PROPERTIES: FREQUENCY, AMPLITUDES AND PHASES
  DOUBLE PRECISION,                           INTENT(IN)    :: T1,T2  ! IN SECONDS
  DOUBLE PRECISION, DIMENSION(D),             INTENT(IN)    :: E_MULTIFLOQUET ! SET OF MULTIMODE FLOQUET ENERGIES, IN Hz, TO AVOID HBAR FACTORS
  COMPLEX*16,       DIMENSION(D,D),           INTENT(IN)    :: U_F_MODES   ! TRANFORMATION MATRIX BETWEEN DRESSED AND BARE BASIS
  COMPLEX*16,       DIMENSION(D_BARE,D_BARE), INTENT(OUT)   :: U           ! EVOLUTION OPERATOR U(T2,T1)

  COMPLEX*16, DIMENSION(D,D) :: U_DIAGONAL
  COMPLEX*16, DIMENSION(D_BARE,D) :: U_AUX

  INTEGER :: MULTIMODE_HARMONICS, n,i,j,m,index0,index1,FIELD_INDEX
  DOUBLE PRECISION, DIMENSION(NM) :: OMEGA_VEC
  INTEGER, DIMENSION(NM) :: N_FLOQUET
  TYPE(HARMONIC_FACTORS), DIMENSION(:),ALLOCATABLE:: U_MODES_n
  DOUBLE PRECISION :: PHASE

  MULTIMODE_HARMONICS   = D/D_BARE

  N_FLOQUET = 0
  OMEGA_VEC = 0

  DO n=2,NM
     FIELD_INDEX  = 2+SUM(MODES_NUM(2:n-1))
     N_FLOQUET(n) = FIELD(FIELD_INDEX)%N_Floquet
     OMEGA_VEC(n) = FIELD(FIELD_INDEX)%OMEGA 
     !     write(*,*) n,N_FLOQUET(n),field_index,modes_num(n)
     IF(modes_num(n).GT.N_FLOQUET(n)+1) THEN
        WRITE(*,*) "TO BUILD THE EXTENDED HAMILTONIAN THE NUMBER OF FLOQUET MODES MUST BE DEFINED"
        WRITE(*,*) "LARGER THAN THE NUMBER OF FIELD MODES"
        INFO = -10
     END IF
  END DO

  ALLOCATE(U_MODES_n(MULTIMODE_HARMONICS))

  DO n=1,MULTIMODE_HARMONICS
     ALLOCATE(U_MODES_n(n)%U(D_BARE,D))
     U_MODES_n(n)%U = 0.0
     ALLOCATE(U_MODES_n(n)%n(NM))
     U_MODES_n(n)%n = 0.0
  END DO

  DO i=0,MULTIMODE_HARMONICS-1
     U_MODES_n(i+1)%n = 0
     index0 = i
     DO j=2,NM
        U_MODES_n(i+1)%n(j)= -N_FLOQUET(j) + MOD(index0,(2*N_FLOQUET(j)+1))
        index0 = INT(index0/(2*N_FLOQUET(j)+1))        
     END DO
  END DO

  i  = 1
  DO index1=1,MULTIMODE_HARMONICS
     DO n=1,D
        j = n
        i = 1 + (index1-1)*D_BARE 
        DO m=1,D_BARE
           U_MODES_n(index1)%U(m,n) = U_F_MODES(i,j)
           i = i + 1           
        END DO
     END DO
  END DO

  index0 = (MULTIMODE_HARMONICS -1)/2 + 1
  U_AUX  =  U_MODES_n(index0)%U


  U_DIAGONAL = 0.0
  DO i=1,D
     U_DIAGONAL(i,i) = EXP(-DCMPLX(0.0,1.0)*E_MULTIFLOQUET(i)*(T2-T1))
  END DO

  DO index1=1,MULTIMODE_HARMONICS
     U_MODES_n(index1)%U = MATMUL(U_MODES_n(index1)%U,U_DIAGONAL)
  END DO

  U = 0.0
  DO index1=1,MULTIMODE_HARMONICS
     PHASE = DOT_PRODUCT(U_MODES_n(index1)%n,omega_vec)*T2
     U = U +  MATMUL(U_MODES_n(index1)%U,TRANSPOSE(CONJG(U_AUX)))*EXP(DCMPLX(0.0,1.0)*PHASE)
  END DO

END SUBROUTINE MULTIMODETIMEEVOLUTINOPERATOR

SUBROUTINE MULTIMODETIMEEVOLUTINOPERATOR_RESTRICTED(D,NF,U_F_MODES,E_MULTIFLOQUET,D_BARE,FIELD,T1,T2,U,INFO) 

  ! TIME EVOLUTION OPERATOR OF A MULTIMODE DRESSED SYSTEM. THE EVOLUTION OPERATOR IS WRITEN IN THE BASIS USED TO EXPRESS THE 
  ! MULTIMODE FLOQUET HAMILTONIAN
  ! U : MATRIX OF AMPLITUED OF PROBABILITIES FOR TRANSITIONS BETWEEN T1 TO T2
  !
  ! BUT THE SUM OVER FLOQUET EIGENVALUES IS RESTRICTED TO THOSE ACCURATELY DEFINED, FAR FROM THE EDGES OF THE MATRICES  
  ! THE SUBROUTINE HAS BEEN DESIGNED FOR 87Rb driven by RF+MW fields, AND IT IS NOT EXPECTED TO WORK FOR OTHER CONFIGURATIONS 
  ! D: THE DIMENSOIN OF THE MULTIMODE FLOQUET IS NOW DEFINED AS THE NUMBER OF FLOQUET MODES USED TO CALCULATE THE EVOLUTION OPERATOR
  ! the point of this is to make a comparison with MULTIMODETIMEEVOLUTIONOPERATOR AND CHECK WHETHER.
  ! MULTIMODETIMEEVOLUTIONOPERATOR SHOULD WORK FOR ANY ATOM/FIELD CONFIGURATION.
  ! MULTIMODETIMEEVOLUTIONOPEATOR_RESTRICTED IS DIFFICULT TO RESTRICT IN GENERAL SINCE WE NEED TO IDENTIFY THE "CENTRAL" MANIFOLD OF FLOQUET 
  ! STATES

!!$  D              (IN)   : DIMENSION OF THE EXTENDED HILBERT SPACE (SIZE OF THE MULTIMODE FLOQUET MATRIX)
!!$  NM             (IN)   : NUMBER OF MODES            
!!$  MODES_NUM      (IN)   : VECTOR (NM) INDICATING THE NUMBER OF HARMONICS OF EACH MODE
!!$  U_F_MODES      (IN)   : TRANSFORMATION, DIMENSOON (D,D) 
!!$  E_MULTIFLOQUET (IN)   : MULTIMODE FLOQUET SPECTRUM
!!$  D_BARE         (IN)   : DIMENSION OF THE BARE HILBERT SPACE
!!$  FIELD          (IN)   : STRUCTURE DESCRIBING THE COUPLINGS
!!$  T1             (IN)   : INITIAL TIME
!!$  T2             (IN)   : FINAL TIME
!!$  U              (OUT)  : TRANFORMATION BETWEEN THE EXTENDED BARE BASIS AND THE FLOQUET STATES, DIMENSION (D_BARE,D)
!!$  INFO           (INOUT): (POSSIBLE) ERROR FLAG


  USE TYPES
  USE SUBINTERFACE_LAPACK


  IMPLICIT NONE
  INTEGER,                                    INTENT(IN)    :: D,D_BARE,NF ! DIMENSION OF THE MULTIMODE FLOQUET SPACE AND THE BARE BASIS
  INTEGER,                                    INTENT(INOUT) :: INFO
  TYPE(MODE),       DIMENSION(NF),     INTENT(IN)    :: FIELD  ! FIELDS PROPERTIES: FREQUENCY, AMPLITUDES AND PHASES
  DOUBLE PRECISION,                           INTENT(IN)    :: T1,T2  ! IN SECONDS
  DOUBLE PRECISION, DIMENSION(D),             INTENT(IN)    :: E_MULTIFLOQUET ! SET OF MULTIMODE FLOQUET ENERGIES, IN Hz, TO AVOID HBAR FACTORS
  COMPLEX*16,       DIMENSION(D,D),           INTENT(IN)    :: U_F_MODES   ! TRANFORMATION MATRIX BETWEEN DRESSED AND BARE BASIS
  COMPLEX*16,       DIMENSION(D_BARE,D_BARE), INTENT(OUT)   :: U           ! EVOLUTION OPERATOR U(T2,T1)

  COMPLEX*16, DIMENSION(:,:),ALLOCATABLE :: U_DIAGONAL
  COMPLEX*16, DIMENSION(:,:),ALLOCATABLE :: U_AUX

  INTEGER :: MULTIMODE_HARMONICS, n,i,j,m,index0,index1,D_r,i_0,i_r,Fup,Fdown
  DOUBLE PRECISION, DIMENSION(NF) :: OMEGA_VEC
  TYPE(HARMONIC_FACTORS), DIMENSION(:),ALLOCATABLE:: U_MODES_n
  DOUBLE PRECISION :: PHASE
  DOUBLE PRECISION :: A,HBAR

  MULTIMODE_HARMONICS   = D/D_BARE
  IF(NF.GT.2) THEN
     D_r                   = (2*FIELD(3)%N_FLOQUET - 2)*(2*(FIELD(2)%N_FLOQUET-2)-1)*D_bare
  END IF
  Fup                   = 2
  Fdown                 = 1
  hbar                  = 1.054E-34 
  A                     = 2.0*4.0*ATAN(1.0)*hbar*3.417341E9

  DO i=1,NF
     OMEGA_VEC(i) = FIELD(i)%OMEGA
  END DO

  ALLOCATE(U_MODES_n(MULTIMODE_HARMONICS))

  DO n=1,MULTIMODE_HARMONICS
     ALLOCATE(U_MODES_n(n)%U(D_BARE,D))
     U_MODES_n(n)%U = 0.0
     ALLOCATE(U_MODES_n(n)%U_r(D_BARE,D_r))
     U_MODES_n(n)%U_r = 0.0
     ALLOCATE(U_MODES_n(n)%n(NF))
     U_MODES_n(n)%n = 0.0
  END DO

  DO i=0,MULTIMODE_HARMONICS-1
     U_MODES_n(i+1)%n = 0
     index0 = i
     DO j=2,NF
        U_MODES_n(i+1)%n(j)= -FIELD(j)%N_FLOQUET + MOD(index0,(2*FIELD(j)%N_FLOQUET+1))
        index0 = INT(index0/(2*FIELD(j)%N_FLOQUET+1))
     END DO
  END DO

  i  = 1
  DO index1=1,MULTIMODE_HARMONICS
     DO n=1,D
        j = n
        i = 1 + (index1-1)*D_BARE 
        DO m=1,D_BARE
           U_MODES_n(index1)%U(m,n) = U_F_MODES(i,j)
           i = i + 1           
        END DO
     END DO
  END DO

  index0 = (MULTIMODE_HARMONICS -1)/2 + 1


  !HOW MANY FLOQUET MODES AND WHICH ONES ARE WE GOING TO USE?
  !
  ! D_r  = (2*FIELD(3)%N_FLOQUET - 2)*(2*(FIELD(2)%N_FLOQUET-2)-1)*RB87%D_bare

  ! (2*FIELD(3)%N_FLOQUET - 2)       : number of mw manifolds, we neglect one up and one down 
  ! (2*(FIELD(2)%N_FLOQUET - 2) - 1) : number of complete (i.e. having Rb87%D_bare) RF manifolds (not in general, it depends on Bdc and omega_rf)
  ! Rb87%D_bare                : number of bare states

  !WITHIN EACH MW MANIFOLD, WE KEEP FLOQUETS MODES WITH INDICES 
  ! I \IN  I_0+1, I_0 +  (2*(FIELD(2)%N_FLOQUET - 2) - 1)*Rb87%D_bare    
  ! WITH 
  ! I_0 = 3*(2*Fdown+2*Fup+2)+1
  ! and the starting index for each manifold is
  ! (2*FIELD(2)%N_FLOQUET+1)*(2*Fdown+1)+(m)*(2*FIELD(2)%N_FLOQUET+1)*(2*Fdown+2*Fup+2)
  ! with m the manifold number, but we keep manifolds m=1,(2*FIELD(3)%N_FLOQUET - 2) 


  ALLOCATE(U_DIAGONAL(D_r,D_r))
  ALLOCATE(U_AUX(D_bare,D_r))


  n  = 1
  U_DIAGONAL = 0.0
  DO i=1,(2*FIELD(3)%N_FLOQUET - 2)
     I_0 = (2*FIELD(2)%N_FLOQUET+1)*(2*Fdown+1) + i*(2*FIELD(2)%N_FLOQUET+1)*(2*Fdown+2*Fup+2)
     DO j=1,(2*(FIELD(2)%N_FLOQUET-2)-1)*D_bare

        i_r = I_0 + 3*(2*Fdown+2*Fup+2)+1 + j 
        U_DIAGONAL(n,n) = EXP(-DCMPLX(0.0,1.0)*E_MULTIFLOQUET(i_r)*(T2-T1))
        n = n + 1
     END DO
  END DO

  DO index1=1,MULTIMODE_HARMONICS
     n = 1
     DO i=1,(2*FIELD(3)%N_FLOQUET - 2)
        I_0 = (2*FIELD(2)%N_FLOQUET+1)*(2*Fdown+1) + i*(2*FIELD(2)%N_FLOQUET+1)*(2*Fdown+2*Fup+2)
        DO j=1,(2*(FIELD(2)%N_FLOQUET-2)-1)*D_bare        
           i_r = I_0 + 3*(2*Fdown+2*Fup+2)+1 + j 
           U_MODES_n(index1)%U_r(:,n) = U_MODES_n(index1)%U(:,i_r)
           n = n+1
        END DO
     END DO
  END DO

  index0 = (MULTIMODE_HARMONICS -1)/2 + 1
  U_AUX  =  U_MODES_n(index0)%U_r


  DO index1=1,MULTIMODE_HARMONICS
     U_MODES_n(index1)%U_r = MATMUL(U_MODES_n(index1)%U_r,U_DIAGONAL)
  END DO

  U = 0.0
  DO index1=1,MULTIMODE_HARMONICS
     PHASE = DOT_PRODUCT(U_MODES_n(index1)%n,omega_vec)*T2
     U = U +  MATMUL(U_MODES_n(index1)%U_r,TRANSPOSE(CONJG(U_AUX)))*EXP(DCMPLX(0.0,1.0)*PHASE)
  END DO


END SUBROUTINE MULTIMODETIMEEVOLUTINOPERATOR_RESTRICTED




