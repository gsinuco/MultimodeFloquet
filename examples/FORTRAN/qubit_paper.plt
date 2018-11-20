reset
#set term wxt 0
set term qt 0
set xlabel "Frequency"
set ylabel "Transition probability"
set yrange [0.25:0.51]
set xrange [0.2:2.2]
set ytics 0,0.1
set xtics 0.2,0.4
unset key
plot "qubit_avgtransition.dat" u 1:3 w l lw 2 #t "Transition Prob."#, "" u 1:2 w l lw 2 lt 3 t "Remain Prob. "

#set term wxt 1
set term qt 1
reset
set pm3d at b map
unset surf
set ylabel "Time"
set xlabel "Frequency"
unset key
set yrange [0:100]
set cbtics 0,0.2
set xtics 0,0.4
set xrange [0.2:2.2]
set xtics 0.2,0.4
set title "Time oscillation probability of |up>"
splot "qubit_oscillation.dat" u 2:1:3 


reset
set terminal postscript eps color enhanced "Helvetica,14"
set output "qubit_RabiOsc.eps"
set multiplot

set origin 0.005,0.0
set size 0.91,0.47
set xlabel "Frequency"
set ylabel "Time average\nTrans. prob."
set yrange [0.25:0.53]
set xrange [0.2:2.2]
set ytics 0,0.1
set xtics 0.2,0.4
set label "(b)" at 0.25,0.5
unset key
plot "qubit_avgtransition.dat" u 1:3 w l lw 2 #t "Transition Prob."#, "" u 1:2 w l lw 2 lt 3 t "Remain Prob. "

unset label
set origin 0.0,0.325
set size 1.0,0.73
set pm3d at b map
unset surf
set ylabel "Time"
set xlabel ""
set label "(a)" at 0.25,90 front textcolor "white" 
unset key
set yrange [0:100]
set cbtics 0,0.2
set ytics 0,20
set xtics 0,0.4
set xrange [0.2:2.2]
set xtics ("" 0.2, "" 0.6, "" 1.0, "" 1.4, "" 1.8, "" 2.2)
#set title "Time oscillation probability of |up>"
splot "qubit_oscillation.dat" u 2:1:3 
unset multiplot
unset output