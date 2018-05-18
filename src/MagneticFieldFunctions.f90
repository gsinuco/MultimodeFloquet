SUBROUTINE CARTESSIANTOPOLAR(A,R_)

  IMPLICIT NONE
  DOUBLE PRECISION,DIMENSION(3), INTENT(IN)  :: A
  DOUBLE PRECISION,DIMENSION(3), INTENT(OUT) :: R_

  DOUBLE PRECISION :: r_xy,x,y,z,r,theta,phi
  !--------- position in polar coordinates --------
  
  x  = A(1)
  y  = A(2)
  z  = A(3)

  r    = sqrt(x*x+y*y+z*z)
  r_xy = sqrt(x*x+y*y)
  if(r_xy.gt.0 .and. z.gt.0) theta = atan(r_xy/z)
  if(r_xy.gt.0 .and. z.lt.0) theta = atan(r_xy/z) + 4.0*ATAN(1.0)
  if(z.eq.0 .and. r_xy.gt.0) theta = 2.0*ATAN(1.0)              
  if(r_xy.eq.0 .and. z.lt.0) theta = 4.0*ATAN(1.0)
  if(r_xy.eq.0 .and. z.gt.0) Then
!     write(*,*) "yo"
     theta = 0.0
  end if
  if(r_xy.eq.0 .and. z.eq.0) theta = 0.0

  if(x.gt.0 .and. y.gt.0) phi = atan(y/x)
  if(x.gt.0 .and. y.lt.0) phi = atan(y/x) + 8.0*ATAN(1.0)
  if(x.lt.0 .and. y.lt.0) phi = atan(y/x) + 4.0*ATAN(1.0)
  if(x.lt.0 .and. y.gt.0) phi = atan(y/x) + 4.0*ATAN(1.0)
  if(y.eq.0 .and. x.gt.0) phi = 0.0*ATAN(1.0)              
  if(y.eq.0 .and. x.lt.0) phi = 4.0*ATAN(1.0)
  if(x.eq.0 .and. y.eq.0) phi = 0.0


  R_(1) = R
  R_(2) = THETA
  R_(3) = PHI

END SUBROUTINE CARTESSIANTOPOLAR

SUBROUTINE POLARComponentsToCartesianComponent(position,B_POL,B_CART,INFO)

  !INFO=1 :POSITION IN CARTESIAN COORDINATES
  !INFO=2 :POSITION IN POLAR COORDINATES

  ! USING MATLAB DEFINITIONS
  ! https://www.mathworks.com/help/phased/ref/sph2cartvec.html
  
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(3), INTENT(OUT)   :: B_CART
  DOUBLE PRECISION, DIMENSION(3), INTENT(IN)    :: B_POL,POSITION
  INTEGER,                        INTENT(INOUT) :: INFO

  DOUBLE PRECISION, DIMENSION(3) :: R_
  DOUBLE PRECISION               :: r,theta,phi
  DOUBLE PRECISION, DIMENSION(3,3) :: M

  IF(INFO.EQ.1) THEN
     
     CALL CARTESSIANTOPOLAR(POSITION,R_)
     r     = R_(1)
     theta = R_(2)
     phi   = R_(3)

     M(1,1) = -sin(phi)
     M(1,2) = -sin(theta)*cos(phi)
     M(1,3) =  cos(theta)*cos(phi)

     M(2,1) =  cos(phi)
     M(2,2) = -sin(theta)*sin(phi)
     M(2,3) =  cos(theta)*sin(phi)

     M(3,1) =  0
     M(3,2) =  cos(theta)
     M(3,3) =  sin(theta)

     B_CART = MATMUL(M,B_POL)

  END IF
 
  IF(INFO.EQ.2) THEN
     
     r     = POSITION(1)
     theta = POSITION(2)
     phi   = POSITION(3)

     M(1,1) = -sin(phi)
     M(1,2) = -sin(theta)*cos(phi)
     M(1,3) =  cos(theta)*cos(phi)

     M(2,1) =  cos(phi)
     M(2,2) = -sin(theta)*sin(phi)
     M(2,3) =  cos(theta)*sin(phi)

     M(3,1) =  0
     M(3,2) =  cos(theta)
     M(3,3) =  sin(theta)

     B_CART = MATMUL(M,B_POL)

  END IF
  
 
END SUBROUTINE POLARCOMPONENTSTOCARTESIANCOMPONENT



MODULE VectorPotentialFunction
  IMPLICIT NONE
CONTAINS
  
  FUNCTION A_phi(r,theta,RingRadius)
    
    USE physical_constants 
    
    IMPLICIT NONE
    DOUBLE PRECISION :: r,theta,RingRadius
    DOUBLE PRECISION :: A_phi
    DOUBLE PRECISION :: kE,Elliptic_K,Elliptic_E   
    

    kE = sqrt(abs(4*RingRadius*r*sin(theta)/(RingRadius**2 + r**2+2*r*RingRadius*sin(theta))))
    !write(*,*) kE,RingRadius,r,theta
    CALL CompleteK(kE,Elliptic_K)
    CALL CompleteE(kE,Elliptic_E)
!    write(*,*) R,THETA,kE, Elliptic_K,ELLIPTIC_E
         
    !Imax =1E-9
    IF(ABS(KE).GT.0) THEN
       A_phi = ((mu_cero/(4*pi))*(4*RingRadius)/sqrt(RingRadius**2+r**2+2*r*RingRadius*sin(theta)))*&
            & (((2-kE**2)*Elliptic_K - 2*Elliptic_E)/kE**2)    
    END IF
!    write(*,*) kE, Elliptic_K,Elliptic_E,A_phi    
    
  END FUNCTION A_phi
END MODULE VectorPotentialFunction


SUBROUTINE RINGFIELD(r,theta,Br,Btheta,RingRadius,M)
  ! This is the magnetic field produced by a current flowing in an ideal ring.
  ! The ring lies at the origin on the xy plane.
  USE VectorPotentialFunction
  USE physical_constants
  
  
  IMPLICIT NONE
  
  DOUBLE PRECISION, INTENT(OUT) :: Br,Btheta
  DOUBLE PRECISION, INTENT(IN)  :: r,theta,RingRadius,M
  DOUBLE PRECISION current,dr,dtheta,Tfinal,w_t,alpha_mf,r_0 

  dr        = 1.0E-10
  dtheta    = 2.0*pi/3600
  
  current   = M/(pi*RingRadius**2) 

  Btheta    = 0.0
  Br        = 0.0
  !write(*,*)r,theta

  IF(r.GT.0 .and. abs(theta).gt.0) THEN
!     Btheta = -current*(A_phi(r,theta,RingRadius) + r*(A_phi(r+dr,theta,RingRadius) -A_phi(r,theta,RingRadius))/dr)/r
     Btheta = -current*(1.0*A_phi(r,theta,RingRadius)/r + 1.0*(A_phi(r+dr,theta,RingRadius) - A_phi(r,theta,RingRadius))/dr)
     !write(*,*) "it's me, Btheta",r,theta,Btheta

  END IF
  IF(r.GT.0 .and. abs(theta).eq.0) THEN
     Br = (mu_cero/(4*pi))*(2*pi*RingRadius**2)*current/(r**2+RingRadius**2)**1.5
   !  write(*,*) "its me theta.eq.0", Br
  END IF
  
  IF(r.GT.0 .AND. theta.NE.2.0*ATAN(1.0) .AND.  theta.NE.6.0*ATAN(1.0) .and. theta.ne.0) THEN  
     Br =  (current/r)*(1.0*A_phi(r,theta,RingRadius)/(tan(theta)) + &
          &   1.0*(A_phi(r,theta+dtheta,RingRadius) - A_phi(r,theta,RingRadius))/(dtheta))

    ! write(*,*) "it's me2, Br",Br
  END IF
  
  IF(r.EQ.0)       Br = mu_cero*current/(2*RingRadius)
  
  IF(r.NE.0 .AND. theta.EQ.2.0*ATAN(1.0)) Br =  0.0

  IF(r.NE.0 .AND. theta.EQ.6.0*ATAN(1.0)) Br =  0.0
  
END SUBROUTINE RINGFIELD

