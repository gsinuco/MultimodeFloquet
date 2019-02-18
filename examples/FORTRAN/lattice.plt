set term wxt 2
#set term qt 2
reset
set view 90,0
set xlabel "{/Symbol w}" offset 0,-1
set zlabel "Quasienergy"
unset key
set xrange [0.2:0.8]
set zrange [-1:1]
splot "oneDlattice.dat" i 1 w d

reset
set term wxt 1
#set term qt 1
set origin 0,0
set size 2.8,1
set terminal postscript eps color enhanced "Helvetica,30"
set output "ShakenLatticeA-C.eps"
set multiplot



set origin 0,0
set size 1,1
set yrange [-0.0:0.5]
set xrange [0.4:0.8]
set xlabel "Driving frequency {/Symbol w}   (V_0)" offset 0,0.25
set ylabel "Dressed energy   (V_0)" offset 1,0
set xtics 0.4,0.1
set ytics 0,0.1
unset key
set label "(a)" at 0.425,0.47
#plot "datos.dat"  u 1:3 w p pt 30 ps 0.5
#plot "oneDlattice.dat"  u 1:3 w p pt 20 ps 0.125
plot "ShakenLattice_spectrum.dat" u 1:3 w p pt 35 ps 0.125,"ShakenLattice_spectrum.dat" u 1:($1==0.40000000596046448?$3:1/0) w p pt 40 ps 0.125 lc rgbcolor "blue","ShakenLattice_spectrum.dat" u 1:($1==0.60000002384185791?$3:1/0) w p pt 40 ps 0.125 lc rgbcolor "black"

unset xtics
unset ytics
unset border
unset xlabel
unset ylabel 
unset label 
set origin 0.41,0.17
set size 0.5,0.38
set autoscale xy
plot 'LatticeSketch.png' binary filetype=png w rgbimage

reset
unset label
set pm3d at b map
unset surf
set origin 1.05,0
set size 0.8,1.1
set ylabel "(n,k)-state" offset -1.5,0.25
set label "n=1" at 1960,1738
set label "n=2" at 1960, 1856
set xlabel "Dressed state" offset 0,0
set xtics ("D" 1664, "D/2" 1792,  "1" 1920)
set ytics (" " 1664, " " 1792, " " 1920)
unset key
set yrange [1536+128:1792+128]
set xrange [1792+128:1536+128]
#set xrange [1536+128:1792+128]
set label "(b)" at 2000,1920 textcolor rgbcolor "black" front 
set cbtics ("0" 0, "0.5" 0.1, "1.0" 0.2)
set cbrange [0:0.2]
splot "shakenlattice_eigenvectors.dat" i 0 u 2:1:($3*$3 )matrix

unset label
set origin 2.0,0
set size 0.8,1.1
unset key
set label "(c)" at 2000,1920  textcolor rgbcolor "black" front
set ylabel "(n,k)-state" offset -1.5,0.25
set label "n=1" at 1960,1738
set label "n=2" at 1960, 1856
set xlabel "Dressed state" offset 0,0
set xtics ("D" 1664, "D/2" 1792,  "1" 1920)
set ytics (" " 1664, " " 1792, " " 1920)
splot "shakenlattice_eigenvectors.dat" i 1 u 2:1:($3*$3) matrix


unset multiplot

unset output
#set term qt 1

