!!$PROGRAM PerturbativeShift2nTest
!!$
!!$  INTEGER                          :: N,resonances,INFO
!!$  DOUBLE PRECISION                 :: frequency
!!$  DOUBLE PRECISION, DIMENSION(2)   :: E,DeltaE
!!$  COMPLEX*16, DIMENSION(2,2)       :: DeltaV,sigma_x,DeltaV_dagger
!!$  
!!$  resonances   = 0
!!$  INFO         = 0
!!$  N            = 2
!!$
!!$  E(1)         =  1.0
!!$  E(2)         = -1.0
!!$  frequency    = 0.1
!!$  DeltaE       = 0.0
!!$
!!$  sigma_x(1,1) = 0.0
!!$  sigma_x(1,2) = 0.5
!!$  sigma_x(2,1) = 0.5
!!$  sigma_x(2,2) = 0.0
!!$  
!!$  Amplitude    = 0.01
!!$
!!$  DeltaV        = (Amplitude/2.0)*sigma_x
!!$  DeltaV_dagger = (Amplitude/2.0)*sigma_x
!!$  Energy_RWA   = sqrt((E(1) - E(2) - frequency)**2 + abs(Amplitude/2.0)**2)
!!$  
!!$  CALL PerturbativeShift2ndOrder(N,E,deltaV,DeltaV_dagger,frequency,deltaE,INFO)
!!$
!!$  write(*,*) E(1)-E(2) - frequency + deltaE(1), E(2) - E(1) + frequency + deltaE(2), Energy_RWA, - Energy_RWA 
!!$  
!!$END PROGRAM PerturbativeShift2nTest

SUBROUTINE PerturbativeShift2ndOrder(N,E,DeltaV,DeltaV_dagger,frequency,DeltaE,INFO)
  !Here we calculate the second order energy shifts (DeltaE) of the spectrum (E), associated to the coupling DeltaV and according to standard second order perturbation theory. 
  
! the Hamiltonian reads:

! H = \sum_i^N E_i |i><i| + exp(%i frequency t) DeltaV + exp(-%i frequency t) DeltaV_dagger (1)

! and the energy shift are:

!!!!!! ORIGINAL EXPRESSION: DeltaE_i = - \sum_m <m|DeltaV|i> <i|DeltaV_dagger|m>/(E_m - E_i - h \nu) + <i|DeltaV|m><m|DeltaV_dagger|i>/(Em - E_i + h \nu) (2)
!!!!!! BUT THEN, MODIFIED  ON THE 19 OCTOBER 2016: I think the correct expression for the energy shift is:

! DeltaE_i = - \sum_m <i|DeltaV_dagger|m><m|DeltaV|i>/(E_m - E_i + h \nu) + <i|DeltaV|m><m|DeltaV_dagger|i>/(Em - E_i - h \nu) (2)

!N   : (Integer) number of energy levels
!E(N): (double precision array of dimension N) Unperturbed energy levels
!DeltaV: (Complex*16 matrix of NxN). Matrix of couplings
!DeltaV_dagger: as in Eq. (1)
!frequency: (Double precision) Frequency of the coupling field, in the same units used for the unperturbed energy levels
!DeltaE: (Double precision, array of dimension N) Energy level shifts
!resonances: (?). To store resonances. for the moment it is an interger = 0. Maybe to be implemented as output in deltaV

!RESONANT CONTRIBUTIONS: RWA. (STILL TO DO. OCT 2015)
! Resonances occurs when the denominator of the contributions is close to zero. As a test to include the effect of 
! resonant driving we:
! 1. Evaluate if the contribution to the shift exceeds the 10% of the current value of the energy
! 2. Check that the coupling is much smaller than the energy gap
! 3. define the resonant shift as RWA:
!
!   Delta ~ sqrt(detuning^2 + coupling^2) - detuning




  IMPLICIT NONE
  INTEGER,                          INTENT(IN)    :: N
  INTEGER,                          INTENT(INOUT) :: INFO
  DOUBLE PRECISION, DIMENSION(N),   INTENT(IN)    :: E
  COMPLEX*16,       DIMENSION(N,N), INTENT(IN)    :: DeltaV,DeltaV_dagger
  DOUBLE PRECISION,                 INTENT(IN)    :: frequency
  DOUBLE PRECISION, DIMENSION(N),   INTENT(INOUT) :: DeltaE
 
  
  INTEGER I,J
  DOUBLE PRECISION pi,AUX
  DOUBLE PRECISION, DIMENSION(N) :: SHIFT_THRESHOLD,auxE

  pi = 4.0*atan(1.0)
  
  SHIFT_THRESHOLD = 0.0
  DO I=1,N
     DO J=1,N
        auxE(I) = E(I)
        IF(J.NE.I) THEN
           auxE(J) = ABS(E(I) - E(J))
           !write(*,*) abs(E(I)-E(J))
        END IF
     END DO
     !write(*,*)
     SHIFT_THRESHOLD(I) = 0.1*MINVAL(abs(auxE),1)
  END DO

!  write(*,*) shift_threshold
!  write(*,*) N,info,E,deltaV
  if(info.ne.0) info=-1
  DELTAE = 0.0
  DO I=1,N
     DO J=1,N
        !IF(((abs(DeltaV(I,J))**2/(E(J) - E(I) - frequency)).LT.E(I)) .AND. &
        !     ((abs(DeltaV(I,J))**2/(E(J) - E(I) + frequency)).LT.E(I))) THEN
        IF(J.NE.I) THEN
           !write(*,*) deltaV(J,I)
           !IF(ABS(DeltaV(J,I)).GT.0.0E-6 .AND. ABS(DeltaV_dagger(I,J)).GT.0.0E-6 .AND. (E(J) - E(I) + frequency)) THEN 
           !              Delta_aux = ABS(ABS((DeltaV(J,I)*DeltaV_dagger(I,J)))/(E(J) - E(I) + frequency))
           !              IF(Delta_aux .LT. 0.1*E(I)) THEN
           !DeltaE(I) = DeltaE(I) -  ABS((DeltaV(J,I)*DeltaV_dagger(I,J)))/(E(J) - E(I) + frequency)
           !IF(I.EQ.3) WRITE(*,*) I,J,E(j),E(I), E(J)-E(I)+FREQUENCY
!!$           ELSE
!!$              DeltaE(I) = DeltaE(I) + 0.5*SQRT((E(J) - E(I) + frequency)**2 + ABS(4.0*(DeltaV(J,I)*DeltaV_dagger(I,J)))) &
!!$                   & 0.5*(E(I)-E(J) - frequency)
!!$           END IF
           !END IF
           !IF(ABS(DeltaV(I,J)).GT.0.0E-6 .AND. ABS(DeltaV_dagger(J,I)).GT.0.0E-6) THEN
           !DeltaE(I) = DeltaE(I) - ABS((deltaV(I,J)*DeltaV_dagger(J,I)))/(E(J) - E(I) - frequency)
           !IF(I.EQ.3) WRITE(*,*) I,J,E(j),E(I), (E(J)-E(I)-FREQUENCY), &
           !     & ABS((deltaV(I,J)*DeltaV_dagger(J,I)))/(E(J)-E(I)-FREQUENCY)
           !END IF
           !           write(*,*) I,J,1E-6*(E(J) - E(I) - frequency)/(2*pi), 1E-6*(E(J) - E(I) + frequency)/(2*pi)
           !DeltaE(i) = E(i)
           AUX = ABS((DeltaV(J,I)*DeltaV_dagger(I,J)))/(E(J) - E(I) + frequency)
           IF(ABS(AUX).LT.shift_threshold(I)) DeltaE(I) = DeltaE(I) + &
                & ABS((DeltaV(J,I)*DeltaV_dagger(I,J)))/(E(I) - E(J) - frequency)

           AUX = ABS((DeltaV(I,J)*DeltaV_dagger(J,I)))/(E(J) - E(I) - frequency)
           IF(ABS(AUX).LT.shift_threshold(I)) DeltaE(I) = DeltaE(I) + &
                & ABS((deltaV(I,J)*DeltaV_dagger(J,I)))/(E(I) - E(J) + frequency)
        END IF
        !           DeltaE(I) = DeltaE(I) - ( (DeltaV(J,I)*DeltaV_dagger(I,J))/(E(J) - E(I) + frequency) +&
        !                & (deltaV(I,J)*DeltaV_dagger(J,I))/(E(J) - E(I) - frequency))
        !write(*,*) abs(DeltaV(I,J)),abs(DeltaV(J,I))
        !ELSE
        !   DeltaE(I) = 0.0
        !END IF
     END DO
  END DO
  !  WRITE(*,*) deltaE(:)
END SUBROUTINE PerturbativeShift2ndOrder


SUBROUTINE PerturbativeShiftRWA(N,E,DeltaV,DeltaV_dagger,frequency,DeltaE,INFO)
  !Here we calculate energy shifts (DeltaE) of the spectrum (E), associated to the coupling DeltaV and according to standard RWA approx.. 
  
! the Hamiltonian reads:

! H = \sum_i^N E_i |i><i| + exp(%i frequency t) DeltaV + exp(-%i frequency t) DeltaV_dagger (1)

! and the energy shift are defined by:

! if E[i] < E[j]
! DeltaE_i = -0.5*sqrt(detuning^2 + rabi^2) + 0.5*detuning
! if E[i] > E[j]
! DeltaE_i =  0.5*sqrt(detuning^2 + rabi^2) - 0.5*detuning

!N   : (Integer) number of energy levels
!E(N): (double precision array of dimension N) Unperturbed energy levels
!DeltaV: (Complex*16 matrix of NxN). Matrix of couplings
!DeltaV_dagger: as in Eq. (1)
!frequency: (Double precision) Frequency of the coupling field, in the same units used for the unperturbed energy levels
!DeltaE: (Double precision, array of dimension N) Energy level shifts
!resonances: (?). To store resonances. for the moment it is an interger = 0. Maybe to be implemented as output in deltaV

!RESONANT CONTRIBUTIONS: RWA. (STILL TO DO. OCT 2015)
! Resonances occurs when the denominator of the contributions is close to zero. As a test to include the effect of 
! resonant driving we:
! 1. Evaluate if the contribution to the shift exceeds the 10% of the current value of the energy
! 2. Check that the coupling is much smaller than the energy gap
! 3. define the resonant shift as RWA:
!
!   Delta ~ sqrt(detuning^2 + coupling^2) - detuning




  IMPLICIT NONE
  INTEGER,                          INTENT(IN)    :: N
  INTEGER,                          INTENT(INOUT) :: INFO
  DOUBLE PRECISION, DIMENSION(N),   INTENT(IN)    :: E
  COMPLEX*16,       DIMENSION(N,N), INTENT(IN)    :: DeltaV,DeltaV_dagger
  DOUBLE PRECISION,                 INTENT(IN)    :: frequency
  DOUBLE PRECISION, DIMENSION(N),   INTENT(INOUT) :: DeltaE
 
  
  INTEGER I,J,signfactor
  DOUBLE PRECISION pi,AUX
  DOUBLE PRECISION, DIMENSION(N) :: SHIFT_THRESHOLD,auxE

  pi = 4.0*atan(1.0)
  
  SHIFT_THRESHOLD = 0.0
  DO I=1,N
     DO J=1,N
        auxE(I) = E(I)
        IF(J.NE.I) THEN
           auxE(J) = ABS(E(I) - E(J))
           !write(*,*) abs(E(I)-E(J))
        END IF
     END DO
     !write(*,*)
     SHIFT_THRESHOLD(I) = 0.1*MINVAL(abs(auxE),1)
  END DO

  
  !write(*,*) shift_threshold
  !write(*,*) N,info,E,deltaV
  if(info.ne.0) info=-1
  DELTAE = 0.0
  DO I=1,N
     DO J=1,N
        IF(J.NE.I) THEN
           IF(E(I)<E(J)) THEN
               
              IF(0.5*(E(J) - E(I) - FREQUENCY) .LT. 0)  SIGNFACTOR = -1
              IF(0.5*(E(J) - E(I) - FREQUENCY) .ge. 0)  SIGNFACTOR =  1
!!$              write(*,*) I,J,0.5*(E(J) - E(I) - FREQUENCY),ABS(DELTAV(J,I)),&
!!$                   & -0.5*SQRT((E(J) - E(I) - FREQUENCY)**2 + &
!!$                   & ABS(2.0*DELTAV(J,I))**2) +  SIGNFACTOR*0.5*(E(J) - E(I) - FREQUENCY), &
!!$                   & ABS((DeltaV(J,I)*DeltaV_dagger(I,J)))/(E(I) - E(J) - frequency) + &
!!$                   & ABS((deltaV(I,J)*DeltaV_dagger(J,I)))/(E(I) - E(J) + frequency)
              
              DELTAE(I) = DELTAE(I) - 0.5*SQRT((E(J) - E(I) - FREQUENCY)**2 + ABS(2.0*DELTAV(J,I))**2) &
                   & + SIGNFACTOR*0.5*(E(J) - E(I) - FREQUENCY)
           ELSE
              IF(0.5*(E(I) - E(J) - FREQUENCY) .LT. 0)  SIGNFACTOR =  1
              IF(0.5*(E(I) - E(J) - FREQUENCY) .GE. 0)  SIGNFACTOR = -1
!!$              write(*,*) I,J,0.5*(E(I) - E(J) - FREQUENCY),ABS(DELTAV(I,J)), &
!!$                   & 0.5*SQRT((E(I) - E(J) - FREQUENCY)**2 + &
!!$                   & ABS(2.0*DELTAV(I,J))**2) + SIGNFACTOR*0.5*(E(I) - E(J) - FREQUENCY), &
!!$                   & ABS((DeltaV(J,I)*DeltaV_dagger(I,J)))/(E(I) - E(J) - frequency) + &
!!$                   & ABS((deltaV(I,J)*DeltaV_dagger(J,I)))/(E(I) - E(J) + frequency)
              DELTAE(I) = DELTAE(I) + 0.5*SQRT((E(I) - E(J) - FREQUENCY)**2 + ABS(2.0*DELTAV(I,J))**2)&
                   & + SIGNFACTOR*0.5*(E(I) - E(J) - FREQUENCY)
           END IF           
        END IF
     END DO
!     write(*,*)
  END DO
  !WRITE(*,*) deltaE(:)
END SUBROUTINE PerturbativeShiftRWA



