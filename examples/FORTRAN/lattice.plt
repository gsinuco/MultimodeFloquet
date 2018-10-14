set view 90,0
set xlabel "{/Symbol w}" offset 0,-1
set zlabel "Quasienergy"
unset key
set xrange [0.2:0.8]
set zrange [-1:1]
splot "oneDlattice.dat" i 1 w d