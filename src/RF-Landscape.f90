SUBROUTINE ROTATION_TO_MAKE_Bdc_Bdc_Z(FIELD,INFO)

  ! rotation to make the static field pointing along the z direction
  ! field, in, field-type, field configuraiotn
  ! info, inout, error flag

  USE physical_constants
  USE FLOQUET
  USE TYPES

  !USE FIELD  
  IMPLICIT NONE
  TYPE(MODE), DIMENSION(MODES_NUM),INTENT(INOUT) :: FIELD
  INTEGER, INTENT(INOUT) :: INFO

  DOUBLE PRECISION :: b0
  COMPLEX*16       ::ax,ay,az,beff

  DOUBLE PRECISION :: cos_alpha,sin_alpha,cos_theta,sin_theta
  DOUBLE PRECISION, DIMENSION(3,3) :: R_ALPHA,R_theta

  COMPLEX*16, DIMENSION(3)   :: B_AUX

  DOUBLE PRECISION :: Bo_x,Bo_y,Bo_z!,B_rf_x,B_rf_y,B_rf_z,B_rf,B_rf_ort,B_rf_par
  INTEGER i

  Bo_x   = FIELD(1)%Bx
  Bo_y   = FIELD(1)%By
  Bo_z   = FIELD(1)%Bz

  !  B_rf_x = FIELD(2)%Bx
  !  B_rf_y = FIELD(2)%By
  !  B_rf_z = FIELD(2)%Bz


  IF(SQRT(Bo_x**2+Bo_y**2) .gt. 0.0) THEN
     cos_alpha = Bo_x/sqrt(Bo_x**2+Bo_y**2)
     sin_alpha = Bo_y/sqrt(Bo_x**2+Bo_y**2)
  ELSE
     cos_alpha = 1.0
     sin_alpha = 0.0
  END IF
  cos_theta = Bo_z/Sqrt(Bo_x**2+Bo_y**2+Bo_z**2) 
  sin_theta = sqrt(Bo_x**2 + Bo_y**2)/sqrt(Bo_x**2+Bo_y**2+Bo_z**2)

  R_alpha = 0.0 
  R_theta = 0.0

  R_alpha(1,1) =  cos_alpha
  R_alpha(2,2) =  cos_alpha
  R_alpha(1,2) =  sin_alpha
  R_alpha(2,1) = -sin_alpha
  R_alpha(3,3) =  1

  R_theta(1,1) =  cos_theta
  R_theta(1,3) = -sin_theta
  R_theta(2,2) =  1.0
  R_theta(3,1) =  sin_theta
  R_theta(3,3) =  cos_theta

!!$  !This is mod B static
!!$  b0  = sqrt( Bo_x**2 + Bo_y**2 + Bo_z**2 );
!!$  !This is B rf orthogonal to the static field.
!!$  ax = (Bo_y*B_rf_z - Bo_z*B_rf_y)/b0;
!!$  ay = (Bo_z*B_rf_x - Bo_x*B_rf_z)/b0;
!!$  az = (Bo_x*B_rf_y - Bo_y*B_rf_x)/b0;        
!!$  beff = sqrt( ax**2 + ay**2 + az**2);
!!$ 
!!$  B_rf_ort = abs(beff)
!!$  B_rf     = sqrt(abs(B_rf_x)**2 + abs(B_rf_y)**2 +abs(B_rf_z)**2)
!!$  B_rf_par = sqrt(abs(abs(B_rf)**2 - abs(beff)**2))
!!$  
!!$  Bo_x = 0.0
!!$  Bo_y = 0.0
!!$  Bo_z = b0
!!$
!!$  B_aux(1) = B_rf_x
!!$  B_aux(2) = B_rf_y
!!$  B_aux(3) = B_rf_z
!!$
!!$  B_AUX = MATMUL(R_theta,MATMUL(R_alpha,B_aux))
!!$
!!$  B_rf_x = sqrt(abs(B_aux(1))**2 + ABS(B_AUX(2))**2) !B_rf_ort
!!$  B_rf_y = 0.0
!!$  B_rf_z = ABS(B_AUX(3))!B_rf_par 

  DO i=1,MODES_NUM

     B_aux(1) = FIELD(i)%Bx
     B_aux(2) = FIELD(i)%By
     B_aux(3) = FIELD(i)%Bz
     !     WRITE(*,*) i,FIELD(i)%Bx, FIELD(i)%By, FIELD(i)%Bz

     B_AUX = MATMUL(R_theta,MATMUL(R_alpha,B_aux))

     FIELD(i)%Bx = B_aux(1) 
     FIELD(i)%By = B_aux(2) 
     FIELD(i)%Bz = B_aux(3) 
     !     WRITE(*,*) i,FIELD(i)%Bx, FIELD(i)%By, FIELD(i)%Bz


  END DO

END SUBROUTINE ROTATION_TO_MAKE_Bdc_Bdc_Z


SUBROUTINE DRESSED_ENERGIES_RWA(D,H_0,V,OMEGA,UR,ENERGY,INFO)
  ! GENERIC PROCEDURE TO OBTAIN DRESSED ENERGIES IN A ROTATING FRAME
  ! D   : MATRIX DIMENSION
  ! H_0 : DIAGONAL MATRIX WITH BARE ENERGIES
  ! V   : COUPLING MATRIX. THE TIME DEPENEDENT COUPLING IS OF THE FORM 
  !       EXP(I OMEGA T) V + EXP(-I OMEGA T) V^DAGGER
  ! UR : INPUT   DIAGONAL MATRIX REPRESENTATION OF THE GENERATOR OF THE ROTATION, J, WITH EXP(-I OMEGA J)
  !      OUTPUT  TRANSFORMATION BETWEEN BARE AND DRESSED BASIS 
  ! ENERGY: DRESSED ENERGIES  

  USE physical_constants
  USE subinterface_lapack 

  IMPLICIT NONE
  INTEGER,                          INTENT(IN)    :: D
  INTEGER,                          INTENT(INOUT) :: INFO
  DOUBLE PRECISION,                 INTENT(IN)    :: OMEGA
  COMPLEX*16,       DIMENSION(D,D), INTENT(IN)    :: H_0,V
  COMPLEX*16,       DIMENSION(D,D), INTENT(INOUT) :: UR
  DOUBLE PRECISION, DIMENSION(D),   INTENT(OUT)   :: ENERGY 

  INTEGER :: j,i,harmonic_plus,harmonic_minus,index_up,index_down

  COMPLEX*16, DIMENSION(D,D) :: VDC,VDC_DAGGER,HROTATING

  !BUILD THE COUPLING MATRIX WITH STATTIC COMPONENETS IN THE ROTATING FRAME  
  !   harmonic_plus  =  2 + (UR(j,j)-U(i,i)) HARMONICS GENERATED FROM THE FACTOR EXP I OMEGA T)
  !   harmonic_minus = -2 + (UR(j,j)-U(i,i)) HARMONICS GENERATED FROM THE FACTOR EXP(-I OMEGA T)

  HROTATING = 0.0
  VDC       = 0.0

  DO j=1,D
     DO i=1,D

        harmonic_plus  =  2 + (UR(j,j)-UR(i,i))
        harmonic_minus = -2 + (UR(j,j)-UR(i,i))

        IF(harmonic_plus.EQ.0) THEN
           IF((i.NE.j) .AND. (j.GT.63) .AND.(i.LE.63)) then
              index_up   = -10+INT((j-63-1)/5) !MOD(j-63-1,5)+1 
              index_down = -10+INT((i-1)/3)    !MOD(i-1,3)+1
           END IF
           IF((i.NE.j) .AND. (j.LE.63) .AND.(i.GT.63)) then
              index_up   = -10+INT((i-63-1)/5) !MOD(i-63-1,5)+1
              index_down = -10+INT((j-1)/3)    !MOD(j-1,3)+1
           END IF

           !IF(((index_up-index_down).eq.-3) .OR. ((index_up-index_down).eq.-1) .OR. &
           !     & ((index_up-index_down).eq.1) .OR. ((index_up-index_down).eq.3)) THEN
           VDC(j,i)        = V(j,i)
           !END IF
        END IF

        IF(harmonic_minus.EQ.0) THEN
           IF((i.NE.j) .AND. (j.GT.63) .AND.(i.LE.63)) then
              index_up   = -10+INT((j-63-1)/5) !MOD(j-63-1,5)+1 
              index_down = -10+INT((i-1)/3)    !MOD(i-1,3)+1
           END IF
           IF((i.NE.j) .AND. (j.LE.63) .AND.(i.GT.63)) then
              index_up   = -10+INT((i-63-1)/5) !MOD(i-63-1,5)+1
              index_down = -10+INT((j-1)/3)    !MOD(j-1,3)+1
           END IF

           !IF(((index_up-index_down).eq.-3) .OR. ((index_up-index_down).eq.-1) .OR. &
           !     & ((index_up-index_down).eq.1) .OR. ((index_up-index_down).eq.3)) THEN
           VDC_DAGGER(j,i) = CONJG(V(i,j))
           !END IF
        END IF

     END DO
     HROTATING(j,j) = H_0(j,j) - OMEGA*0.5*UR(j,j)
  END DO

  HROTATING = HROTATING + VDC + VDC_DAGGER

  INFO  = 0
  CALL LAPACK_FULLEIGENVALUES(HROTATING,D,ENERGY,INFO) ! after this, H_AUX is the matrix of eigenvectors
  UR = HROTATING  

END SUBROUTINE DRESSED_ENERGIES_RWA

SUBROUTINE RF_DRESSED_ENERGIES_RWA(D,FIELD,ENERGY,UF,INFO)
  !    PROCEDURE TO OBTAIN RF DRESSED ENERGIES IN A ROTATING FRAME
  ! D,IN   : MATRIX DIMENSION
  ! FIELD,IN: FIELD CONFIGURAITON
  ! ENERGY, OUT  : DRESSED ENERGY
  ! UF: OUT TRANSFORMATION BETWEEN BARE AND DRESSED BASIS
  !INFO: ERROR FLAG

  USE physical_constants
  USE ATOMIC_PROPERTIES
  USE FLOQUET
  USE ARRAYS             ! Arrays need to define the Floquet Matrix
  USE subinterface       ! To lapack and subroutines for representation of I and J operator
  USE subinterface_lapack 
  USE TYPES
  ! Modules defined in delta_kr.f90:
  !USE funciones

  IMPLICIT NONE
  INTEGER,                         INTENT(IN)    :: d
  DOUBLE PRECISION, DIMENSION(D),  INTENT(OUT)   :: ENERGY
  INTEGER,                         INTENT(INOUT) :: INFO
  COMPLEX*16,DIMENSION(D,D),       INTENT(OUT)   :: UF
  TYPE(MODE), DIMENSION(MODES_NUM),INTENT(IN)    :: FIELD

  INTEGER                            :: r,p,alpha_,n
  DOUBLE PRECISION                   :: nan,OMEGA_MW,W_RF
  COMPLEX*16, DIMENSION(D,D)         :: HAMILTONIAN_RF
  DOUBLE PRECISION, DIMENSION(D)     :: ENERGY_RFDRESSED
  DOUBLE PRECISION, DIMENSION(D+1,D) :: ENERGY_SHIFT
  CHARACTER (LEN=200)  :: NOTE

  if(modes_num.gt.2) then
     OMEGA_MW     = FIELD(3)%OMEGA
  end if
  w_RF         = FIELD(2)%OMEGA
  ENERGY_SHIFT = 0.0
  n =  info
  !----- Counting the number of states  

  F_t = 0
  DO r=1,3
     F_t(r,r) = 1
  END DO
  DO r=4,8
     F_t(r,r) = -1
  END DO

  nan = 0.
  nan = nan / nan

  !  INFO = 0

  !------- Build the Hyperfine, static field and driving terms in the JImImJ basis:      
  !------ Buid the basis of Zeeman shift and transform the terms of the Hamiltonian
  HAMILTONIAN =  FIELD(1)%V
  !------ Define the RF coupling
  !----- Transform the RF coupling to the basis of exact Zeeman splitted states
  H_RF        = FIELD(2)%V/2.0
  H_RF_DAGGER = TRANSPOSE(CONJG(H_RF))
  !H_AUX       = nint(REAL(MATMUL(TRANSPOSE(CONJG(U_ZEEMAN)),MATMUL(I_z+J_z,U_ZEEMAN))))            
  H_AUX = 0.0

  H_AUX(1,1) =  1.0
  H_AUX(2,2) =  0.0
  H_AUX(3,3) = -1.0
  H_AUX(4,4) = -2.0
  H_AUX(5,5) = -1.0
  H_AUX(6,6) =  0.0
  H_AUX(7,7) =  1.0
  H_AUX(8,8) =  2.0

  !  CALL WRITE_MATRIX(real(H_AUX))
  !  CALL WRITE_MATRIX_INT(F_t)
  !WRITE(*,*)
  !CALL WRITE_MATRIX(real(HAMILTONIAN))
  !----- Perform transformation to the RF rotating frame:
  DO r=1,Total_states_LSI
     HAMILTONIAN(r,r) = HAMILTONIAN(r,r) + H_AUX(r,r)*F_t(r,r)*hbar*w_rf/A
  END DO
  DO r=1,total_states_lsi
     DO p=1,total_states_lsi        
        H_w(r,p) = -int( F_t(r,r)*H_AUX(r,r) - F_t(p,p)*H_AUX(p,p)) ! this is the one I used for all calculations for NJP2015, but I think it is wrong. See notes
        ! When using the correct matrix (defined below), the amplitude of the Rabi oscillations 
        ! reduces, which indicates population of undesired states. However, this can be correct. 
        ! It depends on how to define m-m', or m+m'
        ! 19th January 2016: This is correct! I deleted the wrong formulas previously appearing below
        !-----------------------------------------------------  
        !  in the rotating frame, the couplings are:
        !  V(r,p) (new) = exp(i*w_rf*H_w(r,p))*V(r,p) (old)
        !-----------------------------------------------------
     END DO
  END DO
  !  call write_matrix_int(h_w)
  alpha_ = -1
  H_ALPHA        = DCMPLX(0.0,0.0)
  H_ALPHA_DAGGER = DCMPLX(0.0,0.0)
  DO r=1,Total_states_LSI
     DO p=1,Total_states_LSI
        IF((H_w(r,p) .EQ. alpha_)) THEN
           H_ALPHA(r,p)        = H_RF(r,p)                    
        END IF
        IF((H_w(r,p).EQ.-alpha_)) THEN
           H_ALPHA_DAGGER(r,p) = H_RF_DAGGER(r,p)
        END IF
     END DO
  END DO

  !--- Build the Hamiltonian corresponding to RWA of the RF dressing
  HAMILTONIAN = HAMILTONIAN + 1.0*H_ALPHA + 1.0*H_ALPHA_DAGGER 
  HAMILTONIAN_RF = HAMILTONIAN

  !----- CALCULATE RF DRESSED ENERGIES, USING RWA
  !     IF(alpha_.EQ.-4) THEN
  U_RF     = HAMILTONIAN_RF
  CALL LAPACK_FULLEIGENVALUES(U_RF,Total_states_LSI,ENERGY,INFO) ! after this, H_AUX is the matrix of eigenvectors
!  if(info.ne.1) write(*,*) total_states_lsI
  !HAMILTONIAN = MATMUL(TRANSPOSE(CONJG(U_RF)),MATMUL(HAMILTONIAN_RF,U_RF))
  H_AUX = U_RF
  UF    = U_RF
  !  write(*,*) '#RF dressed RWA energies: ', Energy

END SUBROUTINE RF_DRESSED_ENERGIES_RWA

SUBROUTINE RFandMW_DRESSED_ENERGIES_RWA(D,q,FIELD,ENERGY,SHIFT,INFO)

  !D     : Hilbert space dimension
  !q     : must be set to 1
  !FIELD : Field configuration
  !ENERGY: RF-Dressed energies
  !SHIFT : Energy shift due to MW field
  !INFO  : success flag

  USE physical_constants
  USE ATOMIC_PROPERTIES
  USE FLOQUET
  USE ARRAYS             ! Arrays need to define the Floquet Matrix
  USE subinterface       ! To lapack and subroutines for representation of I and J operator
  USE subinterface_lapack
  USE TYPES
  ! Modules defined in delta_kr.f90:
  !USE funciones

  IMPLICIT NONE
  INTEGER,                         INTENT(IN)    :: d,q
  DOUBLE PRECISION, DIMENSION(D),  INTENT(OUT)   :: ENERGY,SHIFT
  INTEGER,                         INTENT(INOUT) :: INFO
  TYPE(MODE), DIMENSION(MODES_NUM),INTENT(IN)    :: FIELD

  INTEGER                            :: r,p,alpha_,n
  DOUBLE PRECISION                   :: nan,OMEGA_MW,W_RF,frequency
  COMPLEX*16, DIMENSION(D,D)         :: HAMILTONIAN_RF
  DOUBLE PRECISION, DIMENSION(D)     :: ENERGY_RFDRESSED,shiftcandidate,SHIFT_
  DOUBLE PRECISION, DIMENSION(D+1,D) :: ENERGY_SHIFT
  INTEGER,          DIMENSION(D,D)   :: UJ

  CHARACTER (LEN=200)  :: NOTE

  OMEGA_MW     = FIELD(3)%OMEGA
  w_RF         = FIELD(2)%OMEGA
  ENERGY_SHIFT = 0.0
  n =  info

  !----- Counting the number of states

  F_t = 0
  DO r=1,3
     F_t(r,r) = 1
  END DO
  DO r=4,8
     F_t(r,r) = -1
  END DO

  nan = 0.
  nan = nan / nan
  !write(*,*) "jere",info
  CALL RF_DRESSED_ENERGIES_RWA(SIZE(ENERGY_RFDRESSED,1),FIELD,ENERGY_RFDRESSED,U_RF,INFO)
  if(info.ne.0)write(*,*) field(1)%Bx,field(2)%By,field(3)%Bz
  HAMILTONIAN_RF = 0.0
  DO r=1,size(energy_rfdressed,1)
     HAMILTONIAN_RF(r,r) = ENERGY_RFDRESSED(r)
  END DO
!  write(*,*) -ENERGY_RFDRESSED(1:3)+ENERGY_RFDRESSED(8) - 2.0 
  !------ SET THE MW FIELD
  H_RF        = FIELD(3)%V/2.0
  H_RF_DAGGER = TRANSPOSE(CONJG(H_RF))

  SHIFT = 0.0

  IF(q.EQ.1) then
     DO alpha_ = -3,-3
        Energy = 0.0
        energy_rfdressed = 0.0
        
        !---- DEFINE A MATRIX TO SHIFT UPWARDS THE LOWER HYPERFINE MANIFOLD BY THE MW FREQUENCY
        H_AUX = 0.0

        DO r=1,2*Fdown+1
           H_aux(r,r) = -1
        END DO
        DO r=2*Fdown+1+1,Total_states_LSI
           H_aux(r,r) =  0!1
        END DO
        
        HAMILTONIAN = HAMILTONIAN_RF - 2.0*0.5*(hbar*OMEGA_MW/A + alpha_*hbar*w_RF/A)*H_AUX
        DO r=1,size(energy_rfdressed,1)
           ENERGY_RFDRESSED(r) = HAMILTONIAN(r,r)
        END DO 
        
        ! FIND THE COUPLINGS IN THE ROTATING FRAME CORRESPONDING TO EXP(%i*(OMEGA+ALPHA*W_RF)t)
        H_ALPHA = DCMPLX(0.0,0.0)
        H_ALPHA_DAGGER = DCMPLX(0.0,0.0)
        DO r=1,Total_states_LSI
           DO p=1,Total_states_LSI
              IF((H_w(r,p) .EQ. alpha_)) THEN
                 H_ALPHA(r,p)        = H_RF(r,p)
              END IF
              IF((H_w(r,p) .EQ. -alpha_)) THEN
                 H_ALPHA_DAGGER(r,p) = H_RF_DAGGER(r,p)
              END IF
           END DO
        END DO
        
        ! TRANSFORM THE MICROWAVE COUPLINGS TO THE DRESSED BASIS
        H_ALPHA        = MATMUL(TRANSPOSE(CONJG(U_RF)),MATMUL(H_ALPHA,U_RF))
        H_ALPHA_DAGGER = MATMUL(TRANSPOSE(CONJG(U_RF)),MATMUL(H_ALPHA_DAGGER,U_RF))
        
        !---- KEEP THE DC COMPONENT OF THE MW COUPLING IN THE MW-ROTATING FRAME

        UJ = 0
        DO r=1,2*Fdown+1
           UJ(r,r) = -1
        END DO
        DO r=2*Fdown+1+1,Total_states_LSI
           UJ(r,r) =  1
        END DO
        frequency = hbar*(OMEGA_MW + alpha_*w_rf)/A

        H_MW = 0.0
        DO r=1,Total_states_LSI
           DO p=1,Total_states_LSI
              IF(UJ(r,r).EQ.-1 .AND. UJ(p,p).EQ. 1) H_MW(r,p) = H_ALPHA(r,p)
              IF(UJ(r,r).EQ. 1 .AND. UJ(p,p).EQ.-1) H_MW(r,p) = H_ALPHA_DAGGER(r,p)
           END DO
        END DO
        
        H_AUX          = HAMILTONIAN + H_MW
        
        !CALL WRITE_MATRIX(ABS(H_AUX))
        CALL LAPACK_FULLEIGENVALUES(HAMILTONIAN,Total_states_LSI,ENERGY_RFDRESSED,INFO) ! ENERGY_RFDRESSED are the RF dressed energies
        CALL LAPACK_FULLEIGENVALUES(H_AUX,Total_states_LSI,ENERGY,INFO) ! ENERGY are the MW+RF dressed energies
!        write(*,*) abs(ENERGY_RFDRESSED(8)- ENERGY_RFDRESSED(1:3)), H_MW(1,8)        
        

        !shift = energy-energy_rfdressed
        ENERGY_SHIFT(alpha_+5,:) = ENERGY!-ENERGY_RFDRESSED
        
     END DO
     SHIFT = 0.0
     DO alpha_=-3,-3
        SHIFT = SHIFT + ENERGY_SHIFT(alpha_+5,:) ! this is actually the energy shift, but energy shift can only be defined in perturvative regimes. when the
                                                 ! MW frequency in close to a resonance, the original character of the states is lost and a shift does not reflect
                                                 ! the effect of the second field. For this reason, and for the Balancing problem, we restrict the sum to -3,
                                                 ! since the MW frequency is near this value. Note also that, with this, ENERGY is the double dressed energy
                                                 ! from which a single state should be chosen.
     END DO
     !CALL RF_DRESSED_ENERGIES_RWA(SIZE(ENERGY_RFDRESSED,1),FIELD,ENERGY_RFDRESSED,U_RF,INFO)
     ENERGY = shift!ENERGY_RFDRESSED
  ELSE
     CALL LAPACK_FULLEIGENVALUES(HAMILTONIAN_RF,Total_states_LSI,ENERGY_RFDRESSED,INFO) ! after this, H_AUX is the matrix of eigenvectors
     ENERGY = ENERGY_RFDRESSED
  END IF
  
END SUBROUTINE RFANDMW_DRESSED_ENERGIES_RWA

SUBROUTINE RFandMW_DRESSED_ENERGIES_STATES_RWA(D,q,FIELD,ENERGY_MODES,U,INFO)

  !    PROCEDURE TO OBTAIN RF DRESSED ENERGIES IN A ROTATING FRAME
  ! D,IN   : MATRIX DIMENSION
  ! q, in, : integer, step parameters
  ! FIELD,IN: FIELD CONFIGURAITON
  ! ENERGY_MODES, OUT  : DRESSED ENERGY
  ! U: OUT TRANSFORMATION BETWEEN BARE AND DRESSED BASIS
  !INFO: ERROR FLAG

  USE physical_constants
  USE ATOMIC_PROPERTIES
  USE FLOQUET
  USE ARRAYS             ! Arrays need to define the Floquet Matrix
  !  USE subinterface       ! To lapack and subroutines for representation of I and J operator
  USE subinterface_lapack 
  USE TYPES
  ! Modules defined in delta_kr.f90:
  !USE funciones

  IMPLICIT NONE
  INTEGER,                              INTENT(IN)    :: d,q
  DOUBLE PRECISION, DIMENSION(D,9),     INTENT(OUT)   :: ENERGY_MODES
  INTEGER,                              INTENT(INOUT) :: INFO
  TYPE(MODE), DIMENSION(MODES_NUM),     INTENT(IN)    :: FIELD
  TYPE(HARMONIC_FACTORS), DIMENSION(9), INTENT(OUT)   :: U

  INTEGER                            :: r,p,alpha_,n
  DOUBLE PRECISION, DIMENSION(D)     :: ENERGY
  DOUBLE PRECISION                   :: nan,OMEGA_MW,W_RF
  COMPLEX*16,       DIMENSION(D,D)   :: HAMILTONIAN_RF,U_F
  DOUBLE PRECISION, DIMENSION(D,D)   :: U_UP
  INTEGER,          DIMENSION(D,D)   :: UJ
  DOUBLE PRECISION, DIMENSION(D)     :: ENERGY_RFDRESSED
  DOUBLE PRECISION, DIMENSION(D+1,D) :: ENERGY_SHIFT
  CHARACTER (LEN=200)  :: NOTE

  !  DO r=1,9
  !     ALLOCATE(U(r)%U(D,D))
  !  END DO

  !  DO r=1,9
  !     write(*,*) size(U(r)%U,1)
  !  END DO

  if(modes_num.gt.2) then
     OMEGA_MW     = FIELD(3)%OMEGA
  end if
  w_RF         = FIELD(2)%OMEGA

  nan = 0.
  nan = nan / nan

  r = info
  r = int((r-1)/5) - 3
  n = info - (r+3)*5 - 1

  !  write(*,*) info,r,n,q
  CALL RF_DRESSED_ENERGIES_RWA(SIZE(ENERGY_RFDRESSED,1),FIELD,ENERGY_RFDRESSED,U_RF,INFO)

  !CALL WRITE_MATRIX(ABS(U_RF))
  HAMILTONIAN_RF = 0.0
  DO r=1,size(energy_rfdressed,1)
     HAMILTONIAN_RF(r,r) = ENERGY_RFDRESSED(r)
  END DO

  !------ SET THE MW FIELD  
  H_RF = 0.0
  H_RF_DAGGER = 0.0
  if(modes_num.gt.2) then
     H_RF        = FIELD(3)%V/2.0
     H_RF_DAGGER = TRANSPOSE(CONJG(H_RF))
  END if
  !  WRITE(*,*) "H_RF"
  !  CALL WRITE_MATRIX(ABS(H_RF))        
  !  WRITE(*,*) "H_RF_DAGGER"
  !  CALL WRITE_MATRIX(ABS(H_RF_DAGGER))

  DO alpha_ = -4,4

     ! FIND THE COUPLINGS IN THE RF-ROTATING FRAME CORRESPONDING TO EXP(%i*(OMEGA+ALPHA*W_RF)t)
     ! exp(i (omega_MW+alpha*omega_RF)t) H_alpha + exp(-i(omega_MW+alpha(omega_RF)) H_alpha_dagger
     H_ALPHA = DCMPLX(0.0,0.0)
     H_ALPHA_DAGGER = DCMPLX(0.0,0.0)
     !     IF(ALPHA_.EQ.0) THEN
     !        CALL WRITE_MATRIX(1.0D0*H_W)
     !        WRITE(*,*) "H_RF"
     !        CALL WRITE_MATRIX(ABS(H_RF))        
     !        WRITE(*,*) "H_RF_DAGGER"
     !        CALL WRITE_MATRIX(ABS(H_RF_DAGGER))
     !     END IF
     DO r=1,Total_states_LSI
        DO p=1,Total_states_LSI
           IF((H_w(r,p) .EQ. alpha_)) THEN
              H_ALPHA(r,p)        = H_RF(r,p)                    
           END IF
           IF((H_w(r,p) .EQ. -alpha_)) THEN
              H_ALPHA_DAGGER(r,p) = H_RF_DAGGER(r,p)
           END IF
        END DO
     END DO

     !---------------------- TRANSFORM THE  MW COUPLINGS IN THE BASIS OF RF DRESSED STATES

     !     IF(ALPHA_.EQ.0) THEN
     !        WRITE(*,*) "H_ALPHA"
     !        CALL WRITE_MATRIX(ABS(H_ALPHA))
     !         WRITE(*,*) "H_ALPHA_DAGGER"
     !       CALL WRITE_MATRIX(ABS(H_ALPHA_DAGGER))
     !     END IF

     H_ALPHA           = MATMUL(TRANSPOSE(CONJG(U_RF)),MATMUL(H_ALPHA,U_RF))
     H_ALPHA_DAGGER    = MATMUL(TRANSPOSE(CONJG(U_RF)),MATMUL(H_ALPHA_DAGGER,U_RF))


     !---- KEEP THE DC COMPONENT OF THE MW COUPLING IN THE MW-ROTATING FRAME
     UJ = 0
     DO r=1,2*Fdown+1
        UJ(r,r) = -1
     END DO
     DO r=2*Fdown+1+1,Total_states_LSI
        UJ(r,r) =  1
     END DO

     H_MW = 0.0
     DO r=1,Total_states_LSI
        DO p=1,Total_states_LSI           
           IF(UJ(r,r).EQ.-1 .AND. UJ(p,p).EQ. 1) H_MW(r,p) = H_ALPHA(r,p)       
           IF(UJ(r,r).EQ. 1 .AND. UJ(p,p).EQ.-1) H_MW(r,p) = H_ALPHA_DAGGER(r,p)       
        END DO
     END DO

     !---- DEFINE A MATRIX TO SHIFT UPWARDS THE LOWER HYPERFINE MANIFOLD BY THE MW FREQUENCY
     !H_AUX = 0.0
     !DO r=1,2*Fdown+1
     !   H_AUX(r,r) = 0.5*(hbar*OMEGA_MW/A + alpha_*hbar*w_RF/A)
     !END DO
     !DO r=2*Fdown+1+1,Total_states_LSI
     !   H_AUX(r,r) = -0.5*(hbar*OMEGA_MW/A + alpha_*hbar*w_RF/A)
     !END DO
     HAMILTONIAN = HAMILTONIAN_RF - 0.5*(hbar*OMEGA_MW/A + alpha_*hbar*w_RF/A)*UJ 



     !---------------------- DRESSED STATES
     H_AUX          = HAMILTONIAN + H_MW!1.0*H_ALPHA + 1.0*H_ALPHA_DAGGER 

     !     IF(alpha_.eq.3) THEN
     !        WRITE(*,*) "MW COUPLING H1"
     !        CALL WRITE_MATRIX(REAL(H_MW))
     !        CALL WRITE_MATRIX(AIMAG(H_MW))
     !        CALL WRITE_MATRIX(A*ABS(H_MW)/(2*PI*HBAR))
     !     END IF

     !CALL LAPACK_FULLEIGENVALUES(H_AUX,Total_states_LSI,ENERGY,INFO) ! after this, H_AUX is the matrix of eigenvectors
     !ENERGY_SHIFT(alpha_+5,:) = ENERGY-ENERGY_RFDRESSED ! DEFINE AS THE DIFFERENCE TO CALCULATE ENERGY SHIFTS
     !ENERGY_SHIFT(alpha_+5,:) = ENERGY ! DEFINED AS THE RF+MW DRESSED ENERGY TO CALCULATE TRANSITION PROBABILITIES/TIME-EVOLUTION OPERATOR
     ENERGY_MODES(:,alpha_+5) = ENERGY
     U(alpha_+5)%U            = H_AUX ! TRANSFORMATION BETWEEN THE DOUBLE DRESSED BASIS TO THE RF-DRESSED BASIS
!     write(*,*) alpha_
!     CALL WRITE_MATRIX(ABS(U(alpha_+5)%U))
     !    if(alpha_.eq.0) write(*,*) ENERGY
     !U(alpha_+5)%U  = MATMUL(U_RF,H_AUX) ! TRANSFORMATION BETWEEN THE DOUBLE DRESSED BASIS TO THE BARE (THOUGH RF ROTATING) BASIS     
     !write(*,*) "#RF mode: ",alpha_,"MW dressing of RF dressed states: ", Energy

  END DO


END SUBROUTINE RFANDMW_DRESSED_ENERGIES_STATES_RWA



SUBROUTINE RFandMW_DRESSED_ENERGIES_STATES_RWA2(D,q,FIELD,ENERGY_MODES,U,INFO)

  !    PROCEDURE TO OBTAIN RF DRESSED ENERGIES IN A ROTATING FRAME
  ! D,IN   : MATRIX DIMENSION
  ! q, in, : integer, step parameters
  ! FIELD,IN: FIELD CONFIGURAITON
  ! ENERGY_MODES, OUT  : DRESSED ENERGY
  ! U: OUT TRANSFORMATION BETWEEN BARE AND DRESSED BASIS
  !INFO: ERROR FLAG

  USE physical_constants
  USE ATOMIC_PROPERTIES
  USE FLOQUET
  USE ARRAYS             ! Arrays need to define the Floquet Matrix
  !  USE subinterface       ! To lapack and subroutines for representation of I and J operator
  USE subinterface_lapack 
  USE TYPES
  ! Modules defined in delta_kr.f90:
  !USE funciones

  IMPLICIT NONE
  INTEGER,                              INTENT(IN)    :: d,q
  DOUBLE PRECISION, DIMENSION(D,9),     INTENT(OUT)   :: ENERGY_MODES
  INTEGER,                              INTENT(INOUT) :: INFO
  TYPE(MODE), DIMENSION(MODES_NUM),     INTENT(IN)    :: FIELD
  TYPE(HARMONIC_FACTORS), DIMENSION(9), INTENT(OUT)   :: U

  INTEGER                            :: r,p,alpha_,n
  DOUBLE PRECISION, DIMENSION(D)     :: ENERGY
  DOUBLE PRECISION                   :: nan,OMEGA_MW,W_RF
  COMPLEX*16,       DIMENSION(D,D)   :: HAMILTONIAN_RF,U_F
  DOUBLE PRECISION, DIMENSION(D,D)   :: U_UP
  INTEGER,          DIMENSION(D,D)   :: UJ
  DOUBLE PRECISION, DIMENSION(D)     :: ENERGY_RFDRESSED
  DOUBLE PRECISION, DIMENSION(D+1,D) :: ENERGY_SHIFT
  CHARACTER (LEN=200)  :: NOTE

  !  DO r=1,9
  !     ALLOCATE(U(r)%U(D,D))
  !  END DO

  !  DO r=1,9
  !     write(*,*) size(U(r)%U,1)
  !  END DO

  if(modes_num.gt.2) then
     OMEGA_MW     = FIELD(3)%OMEGA
  end if
  w_RF         = FIELD(2)%OMEGA

  nan = 0.
  nan = nan / nan

  r = info
  r = int((r-1)/5) - 3
  n = info - (r+3)*5 - 1

  !  write(*,*) info,r,n,q
  CALL RF_DRESSED_ENERGIES_RWA(SIZE(ENERGY_RFDRESSED,1),FIELD,ENERGY_RFDRESSED,U_RF,INFO)

  !CALL WRITE_MATRIX(ABS(U_RF))
  HAMILTONIAN_RF = 0.0
  DO r=1,size(energy_rfdressed,1)
     HAMILTONIAN_RF(r,r) = ENERGY_RFDRESSED(r)
  END DO

  !------ SET THE MW FIELD  
  H_RF = 0.0
  H_RF_DAGGER = 0.0
  if(modes_num.gt.2) then
     H_RF        = FIELD(3)%V/2.0
     H_RF_DAGGER = TRANSPOSE(CONJG(H_RF))
  END if
  !  WRITE(*,*) "H_RF"
  !  CALL WRITE_MATRIX(ABS(H_RF))        
  !  WRITE(*,*) "H_RF_DAGGER"
  !  CALL WRITE_MATRIX(ABS(H_RF_DAGGER))

  DO alpha_ = -4,4

     ! FIND THE COUPLINGS IN THE RF-ROTATING FRAME CORRESPONDING TO EXP(%i*(OMEGA+ALPHA*W_RF)t)
     ! exp(i (omega_MW+alpha*omega_RF)t) H_alpha + exp(-i(omega_MW+alpha(omega_RF)) H_alpha_dagger
     H_ALPHA = DCMPLX(0.0,0.0)
     H_ALPHA_DAGGER = DCMPLX(0.0,0.0)
     !     IF(ALPHA_.EQ.0) THEN
     !        CALL WRITE_MATRIX(1.0D0*H_W)
     !        WRITE(*,*) "H_RF"
     !        CALL WRITE_MATRIX(ABS(H_RF))        
     !        WRITE(*,*) "H_RF_DAGGER"
     !        CALL WRITE_MATRIX(ABS(H_RF_DAGGER))
     !     END IF
     DO r=1,Total_states_LSI
        DO p=1,Total_states_LSI
           IF((H_w(r,p) .EQ. alpha_)) THEN
              H_ALPHA(r,p)        = H_RF(r,p)                    
           END IF
           IF((H_w(r,p) .EQ. -alpha_)) THEN
              H_ALPHA_DAGGER(r,p) = H_RF_DAGGER(r,p)
           END IF
        END DO
     END DO

     !---------------------- TRANSFORM THE  MW COUPLINGS IN THE BASIS OF RF DRESSED STATES

     !     IF(ALPHA_.EQ.0) THEN
     !        WRITE(*,*) "H_ALPHA"
     !        CALL WRITE_MATRIX(ABS(H_ALPHA))
     !         WRITE(*,*) "H_ALPHA_DAGGER"
     !       CALL WRITE_MATRIX(ABS(H_ALPHA_DAGGER))
     !     END IF

     H_ALPHA           = MATMUL(TRANSPOSE(CONJG(U_RF)),MATMUL(H_ALPHA,U_RF))
     H_ALPHA_DAGGER    = MATMUL(TRANSPOSE(CONJG(U_RF)),MATMUL(H_ALPHA_DAGGER,U_RF))


     !---- KEEP THE DC COMPONENT OF THE MW COUPLING IN THE MW-ROTATING FRAME
     UJ = 0
     DO r=1,2*Fdown+1
        UJ(r,r) = -1
     END DO
     DO r=2*Fdown+1+1,Total_states_LSI
        UJ(r,r) =  1
     END DO

     H_MW = 0.0
     DO r=1,Total_states_LSI
        DO p=1,Total_states_LSI           
           IF(UJ(r,r).EQ.-1 .AND. UJ(p,p).EQ. 1) H_MW(r,p) = H_ALPHA(r,p)       
           IF(UJ(r,r).EQ. 1 .AND. UJ(p,p).EQ.-1) H_MW(r,p) = H_ALPHA_DAGGER(r,p)       
        END DO
     END DO

     !---- DEFINE A MATRIX TO SHIFT UPWARDS THE LOWER HYPERFINE MANIFOLD BY THE MW FREQUENCY
     !H_AUX = 0.0
     !DO r=1,2*Fdown+1
     !   H_AUX(r,r) = 0.5*(hbar*OMEGA_MW/A + alpha_*hbar*w_RF/A)
     !END DO
     !DO r=2*Fdown+1+1,Total_states_LSI
     !   H_AUX(r,r) = -0.5*(hbar*OMEGA_MW/A + alpha_*hbar*w_RF/A)
     !END DO
     HAMILTONIAN = HAMILTONIAN_RF - 0.5*(hbar*OMEGA_MW/A + alpha_*hbar*w_RF/A)*UJ 



     !---------------------- DRESSED STATES
     H_AUX          = HAMILTONIAN + H_MW!1.0*H_ALPHA + 1.0*H_ALPHA_DAGGER 

     !     IF(alpha_.eq.3) THEN
     !        WRITE(*,*) "MW COUPLING H1"
     !        CALL WRITE_MATRIX(REAL(H_MW))
     !        CALL WRITE_MATRIX(AIMAG(H_MW))
     !        CALL WRITE_MATRIX(A*ABS(H_MW)/(2*PI*HBAR))
     !     END IF

     CALL LAPACK_FULLEIGENVALUES(H_AUX,Total_states_LSI,ENERGY,INFO) ! after this, H_AUX is the matrix of eigenvectors
     !ENERGY_SHIFT(alpha_+5,:) = ENERGY-ENERGY_RFDRESSED ! DEFINE AS THE DIFFERENCE TO CALCULATE ENERGY SHIFTS
     !ENERGY_SHIFT(alpha_+5,:) = ENERGY ! DEFINED AS THE RF+MW DRESSED ENERGY TO CALCULATE TRANSITION PROBABILITIES/TIME-EVOLUTION OPERATOR
     ENERGY_MODES(:,alpha_+5) = ENERGY
     U(alpha_+5)%U            = H_AUX ! TRANSFORMATION BETWEEN THE DOUBLE DRESSED BASIS TO THE RF-DRESSED BASIS
!     write(*,*) alpha_
!     CALL WRITE_MATRIX(ABS(U(alpha_+5)%U))
     !    if(alpha_.eq.0) write(*,*) ENERGY
     !U(alpha_+5)%U  = MATMUL(U_RF,H_AUX) ! TRANSFORMATION BETWEEN THE DOUBLE DRESSED BASIS TO THE BARE (THOUGH RF ROTATING) BASIS     
     !write(*,*) "#RF mode: ",alpha_,"MW dressing of RF dressed states: ", Energy

  END DO


END SUBROUTINE RFANDMW_DRESSED_ENERGIES_STATES_RWA2

SUBROUTINE MWShiftOfRFDressed(N,FIELD,ENERGY,SHIFT_FULL,INFO)

  ! all defined in Modules.f90
  USE physical_constants
  USE ATOMIC_PROPERTIES
  USE FLOQUET
  USE ARRAYS             ! Arrays need to define the Floquet Matrix
  USE subinterface       ! To lapack and subroutines for representation of I and J operator
  USE subinterface_lapack 
  ! defined in delta_kr.f90
  USE TYPES

  IMPLICIT NONE

  INTEGER,                       INTENT(IN)      :: N    
  INTEGER,                       INTENT(INOUT)   :: INFO
  DOUBLE PRECISION, DIMENSION(N),INTENT(INOUT)   :: SHIFT_FULL,ENERGY
  TYPE(MODE), DIMENSION(MODES_NUM),INTENT(INOUT) :: FIELD

  INTEGER :: alpha_,p,q,r
  DOUBLE PRECISION :: OMEGA_MW,W_RF,FREQUENCY
  DOUBLE PRECISION, DIMENSION(N) :: SHIFT

  !  write(*,*) "omega mw",field(3)%omega,field(1)%Bz,Field(2)%Bx,Field(3)%Bx
  H_RF = 0.0
  !  call write_matrix(real(H_RF))


  !  WRITE(*,*)"# RF LANDSCAPE",INFO
  CALL RF_DRESSED_ENERGIES_RWA(Total_states_LSI,FIELD,ENERGY,U_RF,INFO)
  !  call write_matrix(real(H_AUX))
  ! TRANSFORM THE COUPLING MATRICES TO THE BASIS OF EXACT ZEEMAN SHIFT STATES.

  H_RF = 0.0
  H_RF_DAGGER = 0.0
  if(modes_num.gt.2) then
     H_RF        = 0.5*FIELD(3)%V
     H_RF_DAGGER = TRANSPOSE(CONJG(H_RF))
     OMEGA_MW = FIELD(3)%OMEGA
  END if
  !  CALL WRITE_MATRIX(ABS(H_RF))

  w_rf     = FIELD(2)%OMEGA
  !  WRITE(*,*) W_RF,Omega_mw
  !  call write_matrix(real(H_RF))
  SHIFT_FULL = 0
  IF(w_rf.GT.0D0) THEN
     SHIFT = 0.0
     DO alpha_ = -4,4
        ! FIND THE COUPLINGS IN THE ROTATING FRAME CORRESPONDING TO EXP(%i*(OMEGA+ALPHA*W_RF)t)
        H_ALPHA = DCMPLX(0.0,0.0)
        H_ALPHA_DAGGER = DCMPLX(0.0,0.0)
        DO r=1,Total_states_LSI
           DO p=1,Total_states_LSI
              IF((H_w(r,p) .EQ. alpha_)) THEN
                 H_ALPHA(r,p)        = H_RF(r,p)                    
              END IF
              IF((H_w(r,p) .EQ. -alpha_)) THEN
                 H_ALPHA_DAGGER(r,p) = H_RF_DAGGER(r,p)
              END IF
           END DO
        END DO
        ! TRANSFORM THE MICROWAVE COUPLINGS TO THE DRESSED BASIS
        H_ALPHA        = MATMUL(TRANSPOSE(CONJG(H_AUX)),MATMUL(H_ALPHA,H_AUX))
        H_ALPHA_DAGGER = MATMUL(TRANSPOSE(CONJG(H_AUX)),MATMUL(H_ALPHA_DAGGER,H_AUX))                           
        !EVALUATE PERTURBATIVE SHIFTS OF THE MICROWAVE DRESSING.
        frequency = hbar*(OMEGA_MW + alpha_*w_rf)/A
        !CALL PerturbativeShift2ndOrder(Total_states_LSI,Energy,H_ALPHA,H_ALPHA_DAGGER,frequency,SHIFT,INFO)   
        CALL PerturbativeShiftRWA(Total_states_LSI,Energy,H_ALPHA,H_ALPHA_DAGGER,frequency,SHIFT,INFO)   
        !        WRITE(*,*) ALPHA_,SHIFT(3),frequency,OMEGA_MW , alpha_*w_rf,A
        SHIFT_FULL = SHIFT_FULL+SHIFT
     END DO
  ELSE
     H_ALPHA        = H_RF
     H_ALPHA_DAGGER = H_RF_DAGGER
     frequency      = hbar*OMEGA_MW/A
     SHIFT          = ENERGY
     SHIFT_FULL     = SHIFT
     !     CALL PerturbativeShift2ndOrder(Total_states_LSI,ENERGY,H_ALPHA,H_ALPHA_DAGGER,frequency,SHIFT,INFO)
     !     CALL CLOCK_FREQUENCY_FORTAGH(Bo_z,OMEGA_MW,dressing_field(dressing),OMEGA_CLOCK)
  END IF

END SUBROUTINE MWShiftOfRFDressed



SUBROUTINE RFandMW_POTENTIALLANDSCAPE(N,POSITION_,ENERGY,SHIFT_FULL,INFO)

  !N          : Number of dressed energies
  !POSITION_  : calculated at this position
  !ENERGY     : dressed energies, as many as N
  !SHIFT_FULL : energy shift, it depends on the approach
  !INFO       : success/failure flag

  USE physical_constants
  USE ATOMIC_PROPERTIES
  USE FLOQUET
  USE subinterface_lapack
  USE SUBINTERFACE

  USE TYPES
  USE ARRAYS
  USE QUADRUPOLE_TOP_TRAP
  USE MAGNETIC_BOTTLE

  IMPLICIT NONE
  INTEGER,                       INTENT(IN)  :: N
  DOUBLE PRECISION, DIMENSION(3),INTENT(IN)  :: POSITION_
  DOUBLE PRECISION, DIMENSION(N),INTENT(OUT) ::ENERGY,SHIFT_FULL
  INTEGER, INTENT(INOUT) :: INFO

  DOUBLE PRECISION, DIMENSION(N):: ENERGY_RFDRESSED
  TYPE(MODE) FIELD(MODES_NUM)
  DOUBLE PRECISION :: BDC

  INTEGER M
  INFO = 0
  !write(*,*) "1", info
  CALL SETPARAMETERS_XYZ(POSITION_,FIELD,INFO)
  !write(*,*) "1",info
  CALL COUPLINGMATRICES(FIELD,INFO)
  !write(*,*) "2",info
  SHIFT_FULL = 0.0
  ENERGY     = 0.0
  
  !    do m=1,3
  !       write(*,*) m,field(m)%omega!,FIELD(m)%Bx,FIELD(m)%By,FIELD(m)%Bz, size(field(m)%v,1), position_   
  !   call write_matrix( abs(field(3)%V))
  !    END DO
  !CALL MWShiftOfRFDressed(TOTAL_STATES_LSI,FIELD,ENERGY,SHIFT_FULL,INFO)
  !write(*,*) "1",info
  CALL RFandMW_DRESSED_ENERGIES_RWA(Total_states_LSI, 1,FIELD,ENERGY,SHIFT_FULL,INFO) ! ENERGY SHIFT, MW field,
  if(info.ne.0)  WRITE(*,*) TOTAL_STATES_LSI,FIELD(2)%Bx,FIELD(2)%By,field(2)%Bz

  !CALL RFandMW_DRESSED_ENERGIES_FLOQUET(Total_states_LSI,3*Total_states_LSI,FIELD,ENERGY,INFO)! ENERGY SHIFT, MW field,
  


  ! BDC = A1+3*A3*POSITION(2)**2+(0.5*ALPHA_TOP**2/A1 - 1.5*A3)*(POSITION(1)**2+(z_0-POSITION(3))**2)
  !write(*,*) position(1),position(2),position(3)-z_0,Energy(8),shift_full(8),&
  !& sqrt(abs(field(1)%Bx)**2+abs(field(1)%By)**2+abs(field(1)%Bz)**2)*(mu_B*0.5)*1E-6/(2*pi*hbar),&
  !& sqrt(abs(field(3)%Bx)**2+abs(field(3)%By)**2+abs(field(3)%Bz)**2)*(mu_B*0.5)*1E-3/(hbar*2*pi),&
  !&1E-6*omega_RF/(2*pi),BDC*mu_B*0.5*1E-6/(2*pi*hbar)
  !write(*,*) energy(3),shift_full(3)

END SUBROUTINE RFandMW_POTENTIALLANDSCAPE


SUBROUTINE RFandMW_DRESSED_ENERGIES_FLOQUET(D,q,FIELD,ENERGY,INFO)! ENERGY SHIFT, MW field,


  !D      : Dimension of the bare Hilbert space
  !q      : number of floquet manifolds X D
  !FIELD  : Field configuration
  !ENERGY : Multimode Floquet enerygies, as many as q
  !INFO   : Success flag

  USE ATOMIC_PROPERTIES
  USE FLOQUET
  USE TYPES
  USE SUBINTERFACE
  USE SUBINTERFACE_LAPACK
  USE ARRAYS


  IMPLICIT NONE
  INTEGER,                           INTENT(IN)    :: D,q
  DOUBLE PRECISION, DIMENSION(q),    INTENT(OUT)   :: ENERGY
  INTEGER,                           INTENT(INOUT) :: INFO
  TYPE(MODE), DIMENSION(MODES_NUM),  INTENT(IN)    :: FIELD


  TYPE(ATOM) Rb87
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: E_MULTIFLOQUET
  INTEGER :: INDEX0,INDEX1

  Rb87%D_bare = D

  !    DO index0=1,3
  !       call write_matrix(abs(field(index0)%V))
  !       write(*,*) field(index0)%omega
  !    end do

  CALL MULTIMODEFLOQUETMATRIX(Rb87,FIELD,INFO)
  !call write_matrix(abs(H_floquet))

  ALLOCATE(E_MULTIFLOQUET(SIZE(H_FLOQUET,1)))
  IF(BAND.EQ.0) THEN
     CALL LAPACK_FULLEIGENVALUES(H_FLOQUET,SIZE(H_FLOQUET,1),E_MULTIFLOQUET,INFO)
  END IF

  index0 = (2*N_Floquet(2)+1)*(2*Fdown+1) + N_Floquet(3)*(2*N_Floquet(2)+1)*(2*Fup+2*Fdown+2) &
       & + (2*N_Floquet(2)+1)*(2*Fup+2*Fdown+2)/2 -2*(2*Fup+2*Fdown+2)
  index1 =  index0 + 3*(2*Fup+1+2*Fdown+1)
  !    write(*,*) index0,index1!,index0-index1,size(H_floquet,1),SIZE(ENERGY,1),Fup,Fdown,N_floquet(2),N_floquet(3)!E_MULTIFLOQUET
  !              !           WRITE(*,*) 1E-6*(FIELD(3)%OMEGA-2*A/HBAR)/(2*pi),REAL(E_MULTIFLOQUET(index0:index1))
  !
  !    index0  =1
  !    index1  =8
  ENERGY = E_MULTIFLOQUET(index0:index1)
  DEALLOCATE(H_FLOQUET)


END SUBROUTINE RFandMW_DRESSED_ENERGIES_FLOQUET


SUBROUTINE STATICMAGNETICFIELD(y_dim,n_dim,k,FmF,X_AUX,BDC,INFO)
!SUBROUTINE GRADIENTENERGYLANDSCAPE(y_dim,n_dim,k,FmF,X_AUX,Y_AUX,INFO)

  !N          : Number of dressed energies
  !POSITION_  : calculated at this position
  !ENERGY     : dressed energies, as many as N
  !SHIFT_FULL : energy shift, it depends on the approach
  !INFO       : success/failure flag

  USE physical_constants
  USE ATOMIC_PROPERTIES
  USE FLOQUET
  USE subinterface_lapack
  USE SUBINTERFACE

  USE TYPES
 ! USE ARRAYS
  !USE QUADRUPOLE_TOP_TRAP
  !USE MAGNETIC_BOTTLE

  INTEGER,                            INTENT(IN)    :: y_dim,N_DIM,K,FmF
  DOUBLE PRECISION, DIMENSION(N_DIM), INTENT(IN)    :: X_AUX
  DOUBLE PRECISION,                   INTENT(OUT)   :: BDC
  INTEGER,                            INTENT(INOUT) :: INFO

  INTEGER M
  DOUBLE PRECISION :: rho,theta,phi
  !double precision, dimension(3) :: position
  TYPE(MODE) FIELD(MODES_NUM)

  INFO = 0
  
  !rho   = X_AUX(1)
  !theta = X_AUX(2)
  !phi   = X_AUX(3)
  
  !POSITION(1) = rho*sin(theta)*cos(phi)
  !POSITION(2) = rho*sin(theta)*sin(phi)
  !POSITION(3) = rho*cos(theta) + z_0
  
  
  CALL SETPARAMETERS_XYZ(X_AUX,FIELD,INFO)
  BDC = SQRT(ABS(FIELD(1)%Bx)**2 + ABS(FIELD(1)%By)**2 + ABS(FIELD(1)%Bz)**2)

END SUBROUTINE STATICMAGNETICFIELD


SUBROUTINE GRADIENTSTATICMAGNETICFIELD(y_dim,n_dim,k,FmF,X_AUX,GBDC,INFO)
   !SUBROUTINE GRADIENTENERGYLANDSCAPE(y_dim,n_dim,k,FmF,X_AUX,Y_AUX,INFO)

  !N          : Number of dressed energies
  !POSITION_  : calculated at this position
  !ENERGY     : dressed energies, as many as N
  !SHIFT_FULL : energy shift, it depends on the approach
  !INFO       : success/failure flag

  USE physical_constants
  USE ATOMIC_PROPERTIES
  USE FLOQUET
  USE subinterface_lapack
  USE SUBINTERFACE

  USE TYPES
 ! USE ARRAYS
  !USE QUADRUPOLE_TOP_TRAP
  !USE MAGNETIC_BOTTLE

  INTEGER,                            INTENT(IN)    :: y_dim,N_DIM,K,FmF
  DOUBLE PRECISION, DIMENSION(N_DIM), INTENT(IN)    :: X_AUX
  DOUBLE PRECISION,                   INTENT(OUT)   :: GBDC
  INTEGER,                            INTENT(INOUT) :: INFO

  INTEGER M
  DOUBLE PRECISION :: rho,theta,phi
  DOUBLE PRECISION, DIMENSION(3) :: POSITION_L, POSITION_R
  DOUBLE PRECISION :: DELTAZ
  TYPE(MODE) FIELD(MODES_NUM)

  INFO = 0
  DELTAZ = 1.0E-6

  POSITION_L = X_AUX   
  CALL SETPARAMETERS_XYZ(POSITION_L,FIELD,INFO)
  BDC_L = SQRT(ABS(FIELD(1)%Bx)**2 + ABS(FIELD(1)%By)**2 + ABS(FIELD(1)%Bz)**2)


  POSITION_R = X_AUX
  POSITION_R(3) = X_AUX(3) + DELTAZ
  CALL SETPARAMETERS_XYZ(POSITION_R,FIELD,INFO)
  BDC_R = SQRT(ABS(FIELD(1)%Bx)**2 + ABS(FIELD(1)%By)**2 + ABS(FIELD(1)%Bz)**2)

  GBDC = (BDC_R-BDC_L)/deltaz
  !write(*,*) gbdc
  
END SUBROUTINE GRADIENTSTATICMAGNETICFIELD
