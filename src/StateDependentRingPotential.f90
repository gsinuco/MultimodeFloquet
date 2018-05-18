PROGRAM RINGPOTENTIAL

! PROGRAM TO GENERATE SURFACE DATA TO SHOW SCHEMATICALLY A STATE-DEPENDENT POTENTIAL LANDSCAPE  BY RF DRESSING THE INDUCTIVELY COUPLED RING TRAP

INTEGER m,n,N_
DOUBLE PRECISION r, theta,x,y,V_plus,V_minus,pi

pi = 4.0*atan(1.0)
N_ = 256
DO m=1,N_+1
   theta = (m-1.0)*2.0*pi/N_
   DO n=1,N_
      r = 5.0+ (n-1.0)*5.0/N_
      V_plus  = ((r-7.5)**2) + 3.0*cos(theta + pi/5.0 + pi)
      V_minus = ((r-7.5)**2) + 3.0*cos(theta - pi/5.0 + 2.0*pi)
      write(*,*) r*cos(theta),r*sin(theta),V_plus,V_minus
   END DO
   write(*,*)
END DO
reset
set cbrange [0:9]
#set size square
unset key
set colorbox
set cbtics -0,3
unset label 
#set label "{/Symbol m}K" at 11,11 font "Helvetica,40"
unset box
unset border
unset xtics
unset ytics
set label "(b)" at -11,11 font "Helvetica,40"
#set terminal postscript eps color enhanced "Helvetica,30" size 3,3
set terminal png transparent font "Helvetica,27" size 680,680
set output "V_Ring_up.png"
set pm3d at b map
unset surf
splot 'datos.dat' u 1:2:($3+3) 
unset output

unset label 
set cbrange [0:9]
#set size square
unset key
unset colorbox
unset label 
unset box
unset border
unset xtics
unset ytics
set label "(a)" at -11,11 font "Helvetica,40"
set output "V_Ring_down.png"
splot 'datos.dat' u 1:2:($4+3) w pm3d
unset output

set term wxt 1

END PROGRAM RINGPOTENTIAL
