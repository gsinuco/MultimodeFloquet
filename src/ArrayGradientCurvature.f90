! two subroutines to obtain the gradient and curvature of a 1d function defined as an array

MODULE FUNATX_
  TYPE :: FUNATX
     INTEGER          :: index
     DOUBLE PRECISION :: fn,grad,curv
  END TYPE FUNATX
END MODULE FUNATX_

SUBROUTINE QUICK_SORT_I_T(v,index_t,N)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  
  !INTEGER, DIMENSION(N),INTENT(INOUT) :: v
  DOUBLE PRECISION, DIMENSION(N),INTENT(INOUT) :: v
  INTEGER, DIMENSION(N),INTENT(INOUT) :: index_t

  INTEGER, PARAMETER :: NN=2500, NSTACK=500
  real :: a, cpu
  integer :: trial
  INTEGER :: k,i,j,jstack,l,r,istack(NSTACK),indice
  
  jstack=0
  l=1
  r=n
  
  do i =1,N
     index_t(i) = i
  end do
  do
     if (r-l < NN) then
        do j=l+1,r
           indice = index_t(j)
           a=v(j)
           do i=j-1,l,-1
              if (v(index_t(i)) <= a) exit
              index_t(i+1) = index_t(i)
           end do
           index_t(i+1) = indice
        end do
        if (jstack == 0) exit
        r=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
     else
        k=(l+r)/2
        call swap(index_t(k),index_t(l+1))
        if (v(index_t(l)) > v(index_t(r))) call swap(index_t(l),index_t(r))
        if (v(index_t(l+1))>v(index_t(r))) call swap(index_t(l+1),index_t(r))        
        if (v(index_t(l))>v(index_t(l+1))) call swap(index_t(l),index_t(l+1))
        
        i=l+1
        j=r
        indice = index_t(l+1)
        a=v(indice)
        do
           do
              i=i+1
              if (v(index_t(i)) >= a) exit
           end do
           do
              j=j-1
              if (v(index_t(j)) <= a) exit
           end do
           if (j < i) exit
           call swap(index_t(i),index_t(j))
        end do
        index_t(l+1)=index_t(j)
        index_t(j)=indice
        jstack=jstack+2
        if (jstack > NSTACK) then
           print *, "NSTACK too small", jstack
        end if
        if (r-i+1 >= j-l) then
           istack(jstack)=r
           istack(jstack-1)=i
           r=j-1
        else
           istack(jstack)=j-1
           istack(jstack-1)=l
           l=i
        end if
     end if
  end do
  
END SUBROUTINE QUICK_SORT_I_T

subroutine swap(a,b)
  !real, intent(inout) :: a, b
  !real :: c
  integer, intent(inout) :: a,b
  integer :: c
  c = a; a = b; b = c
end subroutine swap


SUBROUTINE firstderivative(N,y,dy,info)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  DOUBLE PRECISION, DIMENSION(N), INTENT(IN) :: Y
  DOUBLE PRECISION, DIMENSION(N), INTENT(OUT) :: DY
  INTEGER, INTENT(INOUT) :: INFO

  INTEGER i
 
  
  DY(1) = Y(2) - Y(1)
  i = 1
  DO i=2,N-1
     DY(i) = (Y(i+1) - Y(i-1))/2.0
  END DO
  i = N
  DY(N) = Y(N) - Y(N-1)
  INFO = 1
  IF(N.GT.2) INFO = 0

END SUBROUTINE firstderivative


SUBROUTINE SECONdDerivative(N,y,ddy,info)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  DOUBLE PRECISION, DIMENSION(N), INTENT(IN) :: Y
  DOUBLE PRECISION, DIMENSION(N), INTENT(OUT) :: dDY
  INTEGER, INTENT(INOUT) :: INFO
  
  INTEGER i;
               
  DDY(1) = (Y(3)-2.0*y(2)+Y(1))
  DDY(2) = (Y(4)-2.0*y(3)+y(2))
  DO i=3,N-2
     DDY(i) = (-Y(i+2) + 16*Y(i+1) - 30*Y(i) + 16*Y(i-1) - Y(i-2))/(12.0)
     !if i==j Hessian[i][j] = (-functiontocall(x_h1) + 16*functiontocall(x_h2) -30*functiontocall(x) + 16*functiontocall(x_h3) - functiontocall(x_h4))/(12*step[i]*step[i]);
     
  END DO
  DDY(N-1) =(Y(i+1) - 2*y(i) + y(i-1))
  DDY(N)   = 0.5*Y(N-2) + 0.5*Y(N) - Y(N-1) 

  INFO = 1
  IF(N.GT.2) INFO = 0
!  DO i=1,N
!     write(*,*)  ddY(i)
!  END DO
  
END SUBROUTINE SECONdDerivative

SUBROUTINE EQUILIBRIUMPOINTSOFANARRAY(N,y,min_,info)

  !We find the absolute equilibrium points of an array representing a 1D potential energy landscape
  ! 
  USE FUNATX_
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  DOUBLE PRECISION, DIMENSION(N), INTENT(IN) :: Y
  INTEGER, INTENT(OUT) :: min_
  INTEGER, INTENT(INOUT) :: INFO
  
  DOUBLE PRECISION, DIMENSION(N) :: ddy
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: y_aux
  INTEGER, DIMENSION(N) ::index_m,index_minArray
  DOUBLE PRECISION :: min__,dy

  INTEGER i,j,l,ell, i_l,max_
  LOGICAL, DIMENSION(:), ALLOCATABLE :: true_minArray
  integer,dimension(1) :: inte

  ALLOCATE(y_aux(N))

  !NORMALIZE THE FUNCTION TO MAKE ITS RANGE BETWEEN [0,1]
  min_  = MINLOC(Y,DIM=1)
  max_  = MAXLOC(Y,DIM=1)
  Y_aux = Y - Y(min_)
  Y_aux = Y_AUX/Y_AUX(maxloc(Y_aux,dim=1))

  !CALCULATE THE ARRAY OF SECOND DERIVATIVES
  CALL SECONdDerivative(N,y_aux,ddy,info)  

  !MAKE A LIST OF POSITIONS AT WHICH THE SECOND DERIVATE IS LARGER THAN ZERO
  !START AT 2 AND END AT N-2 BECAUSE ONLY IN THAT INTERVAL, THE SECOND DEVIRATE IS WELL DEFINED.
  index_m  = 0
  j = 1
  DO i=2,N-2
     IF(ddy(i).GT.0) THEN
        index_m(j) = i
        j          = j+1 
     END IF
  END DO
  ! AFTER THIS LOOP, J-1 IS THE NUMBER OF POINTS IN THE DOMAIN WITH A POSITIVE SECOND DERIVATIVE

  ! FIND THE LIMITING POSITIONS OF THE INTERVALES (IN THE DOMAIN) WHERE THE SECOND DERIVATE IS POSITIVE
  ! AND IDENTIFY THE MINIMUN OF EACH INTERVAL
  i_l = 1
  ell = 1
  DO i=1,j-2
     IF(ABS(index_m(i+1) - index_m(i)).NE.1) THEN
        index_minArray(ell) = index_m(i_l)-1+MINLOC(y_aux(index_m(i_l):index_m(i)),1)
        i_l = i +1
        ell = ell +1
     END IF
  END DO
  IF(index_m(i_l).LT.index_m(j-1)) THEN
     index_minArray(ell) = index_m(i_l)-1+MINLOC(y_aux(index_m(i_l):index_m(j-1)),1)
     ell = ell +1
  END IF


  ! NOW CHECK THAT THE IDENTIFIED MINIMA ARE ACTUAL MINIMUN, COMPARING THE VALUES OF THE FUNCITION AT THE LEFT AND RIGHT. BOTH SHOULD BE LARGER THAN THE 
  ! VALUE AT THE POSITION
  ell = ell -1  !(number of identified minima)
  ALLOCATE(true_minArray(ell))
  true_minArray = .FALSE. ! ASSUME THEY ARE NOT
  DO i=1,ell
     IF(y_aux(index_minArray(i)).LT.y_aux(index_minArray(i)+1) .AND. y_aux(index_minArray(i)).LT.y_aux(index_minArray(i)-1)) THEN
        ! YES, ell is a minimum
 !       WRITE(*,*) "YES", i,ell
        true_minArray(i) = .TRUE.
     ELSE
        ! NO, ell is not a minimun
  !      WRITE(*,*) "NO", i,ell
        true_minArray(i) = .FALSE.
     END IF
  END DO
  
  ! from all the true minima, select the one who has the smallest value 
  inte = index_minArray(MINLOC(Y_AUX(index_minArray(1:ell)),MASK=true_minArray(1:ell) .EQV. .TRUE.) )
  min_ = inte(1)

  !NOW CHECK THAT THERE IS NO A MINIMA AT THE EDGES OF THE DOMAIN, I.E. Y(1),Y(2) OR Y(N-1), Y (N-2)
  dy = ABS(Y(2)-Y(1))
!  WRITE(*,*) Y(1),Y(min_)
!  WRITE(*,*) abs(ABS(Y(2)-Y(1)))/ABS(0.5*(Y(min_+1)-Y(min_-1)))
  IF( Y(1)  .LT.Y(min_) .AND. ABS(Y(2)-Y(1))/ABS(0.5*(Y(min_+1)-Y(min_-1))).LT.10) min_ = 1
  IF(Y(2)  .LT.Y(min_) .AND. ABS(0.5*(Y(3)-Y(1)))/ABS(0.5*(Y(min_+1)-Y(min_-1))).LT.10) min_ = 2
  IF(Y(N-1).LT.Y(min_) .AND. ABS(0.5*(Y(N)-Y(N-2)))/ABS(0.5*(Y(min_+1)-Y(min_-1))).LT.10) min_ = N-1
  IF(Y(N)  .LT.Y(min_) .AND. ABS(Y(N)-Y(N-1))/ABS(0.5*(Y(min_+1)-Y(min_-1))).LT.10) min_ = N


END SUBROUTINE EQUILIBRIUMPOINTSOFANARRAY

!!$PROGRAM EQUILIBRIUMPOINTSOFANARRAY_TEST
!!$  
!!$  !We find the absolute equilibrium points of an array representing a 1D potential energy landscape
!!$  ! 
!!$  !After calculating the first derivate, we sort dy. 
!!$  USE FUNATX_
!!$  IMPLICIT NONE
!!$
!!$  INTEGER N,i,info,min_,max_
!!$  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Y,DY,DDY
!!$
!!$  DOUBLE PRECISION :: pi,phase
!!$
!!$  pi = 4.0*atan(1.0)
!!$
!!$  info = 0
!!$  N = 128
!!$  ALLOCATE(Y(N))
!!$  ALLOCATE(DY(N))
!!$  ALLOCATE(DDY(N))
!!$
!!$  phase = 0.0*2.5*pi/2.0
!!$  DO i=1,N
!!$     Y(i) = (I-n/2)**2!0.36*cos(2.0*i*pi/N + phase) + 0.28*cos(6.0*i*pi/N+phase) + 0.35*sin(2.0*i*pi/N+phase) &
!!$          !&+ 0.5*cos(7.0*i*pi/N+phase) !+ 0.75*sin(3.65*i*pi/N)
!!$     write(*,*) Y(i)
!!$  END DO
!!$  write(*,*)
!!$  write(*,*)
!!$  
!!$  CALL EQUILIBRIUMPOINTSOFANARRAY(N,Y,min_,info)
!!$  write(*,*) "The minimun is located at:", min_,Y(min_)
!!$  CALL EQUILIBRIUMPOINTSOFANARRAY(N,-Y,max_,info)
!!$  write(*,*) "The maximun is located at:",max_,Y(max_)
!!$
!!$END PROGRAM EQUILIBRIUMPOINTSOFANARRAY_TEST
