set term qt 2
reset
set view 90,0
set xlabel "{/Symbol w}" offset 0,-1
set zlabel "Quasienergy"
unset key
set xrange [0.2:0.8]
set zrange [-1:1]
splot "oneDlattice.dat" i 1 w d

set term qt 1
set yrange [-0.0:0.5]
set xrange [0.37:0.8]
set xlabel "{/Symbol w}" offset 0,-1
set ylabel "Dressed energies"
set xtics 0.4,0.1
plot "datos.dat"  u 1:3 w p pt 30 ps 0.5
set terminal postscript eps color enhanced "Helvetica,18"
set output "shakenlattice.eps"
plot "oneDlattice.dat"  u 1:3 w p pt 20 ps 0.25
unset output
set term qt 1