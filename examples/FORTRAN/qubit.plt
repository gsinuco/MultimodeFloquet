reset
#set term wxt 0
set term qt 0
set xlabel "Frequency"
set ylabel "Transition probability"
set yrange [0.27:0.8]
set ytics 0,0.1
set xtics 0,0.2
plot "qubit_avgtransition.dat" u 1:3 w l lw 2 t "Transition Prob.", "" u 1:2 w l lw 2 lt 3 t "Remain Prob. "

#set term wxt 1
set term qt 1
reset
set pm3d at b map
unset surf
set xlabel "Time"
set ylabel "Frequency"
unset key
set xrange [0:100]
set cbtics 0,0.2
set ytics 0,0.4
set yrange [0.2:2.2]
set title "Time oscillation probability of |up>"
splot "qubit_oscillation.dat" u 1:2:3 

#set term wxt 8
set term qt 8
reset
set pm3d at b map
unset surf
set xlabel "Time"
set ylabel "Frequency"
unset key
set xrange [0:100]
set cbtics 0,0.2
set ytics 0,0.4
set yrange [0.2:2.2]
set title "Time oscillation probability of |up>"
splot "qubit_oscillation_DRIVER.dat" u 1:2:3 

#set term wxt 9
set term qt 9
reset
set pm3d at b map
unset surf
set xlabel "Time"
set ylabel "Frequency"
unset key
set xrange [0:100]
set cbtics 0,0.2
set ytics 0,0.4
set yrange [0.2:2.2]
set cbrange [0:1]
set title "Time oscillation probability of |up>"
splot "spin_oscillation_DRIVER.dat" u 1:2:3 


reset
#set term wxt 2
set term qt 2
set pm3d at b map
unset surf
unset key
set xlabel "Time"
set ylabel "Frequency"
set xrange [0:1600]
set xtics 0,400
set yrange [0.035:0.1]
set ytics 0.04,0.02
set cbtics 0,0.2
set label "(a)" at 100,0.09 front textcolor "black" font "Helvetica,26"
splot "qubit_bareoscillation.dat" u 2:1:4
set terminal postscript eps color enhanced "Helvetica,18"
set output "bichromaticqubit_bare.eps"
splot "qubit_bareoscillation.dat" u 2:1:4
unset output



reset
#set term wxt 3
set term qt 3
unset label
set pm3d at b map
unset surf
unset key
set xlabel "Time"
set ylabel "Frequency"
set xrange [0:1600]
set xtics 0,400
set yrange [0.035:0.1]
set ytics 0.04,0.02
set cbtics 0,0.2
set label "(b)" at 100,0.09 front textcolor "white" font "Helvetica,26"
splot "qubit_dressedoscillation.dat" u 2:1:4
set terminal postscript eps color enhanced "Helvetica,18"
set output "bichromaticqubit_dressed.eps"
splot "qubit_dressedoscillation.dat" u 2:1:4
unset output


reset
#set term wxt 4
set term qt 4
set xlabel "Frequency"
set ylabel "Transition probability (SP)"
set yrange [0.27:0.8]
set ytics 0,0.1
set xtics 0,0.2
plot "qubit_avgtransition_SP.dat" u 1:3 w l lw 2 t "Transition Prob.", "" u 1:2 w l lw 2 lt 3 t "Remain Prob. "

#set term wxt 5
set term qt 5
reset
set pm3d at b map
unset surf
set xlabel "Time"
set ylabel "Frequency"
unset key
set xrange [0:100]
set cbtics 0,0.2
set ytics 0,0.4
set yrange [0.2:2.2]
set title "Time oscillation probability of |up> (SP)"
splot "qubit_oscillation_SP.dat" u 1:2:3 


reset
#set term wxt 6
set term qt 6
set pm3d at b map
unset surf
unset key
set xlabel "Time"
set ylabel "Frequency"
set title "Bichromatically driven qubit\nBare basis (SP)"
splot "qubit_bareoscillation_SP.dat" u 2:1:3 

reset
#set term wxt 7
set term qt 7
set pm3d at b map
unset surf
unset key
set title "Bichromatically driven qubit\nDressed basis (SP)"
set xlabel "Time"
set ylabel "Frequency"
splot "qubit_dressedoscillation_SP.dat" u 2:1:3 



reset
#set term wxt 10
set term qt 10
unset key
set xlabel "Frequency"
set title "Bichromatically driven qubit\nBare basis (SP)"
plot "qubit_bareoscillation_DRIVER.dat" u 1:3   w lp

reset
#set term wxt 11
set term qt 11
unset key
set title "Bichromatically driven qubit\nDressed basis (SP)"
set xlabel "Frequency"
plot "qubit_dressedoscillation_DRIVER.dat" u 1:3 w lp

reset
#set term wxt 12
set term qt 12
unset key
set xlabel "Frequency"
set title "Bichromatically driven 87Rb\nBare basis (SP)"
plot "Rb87_bareoscillation_DRIVER.dat" u 1:2  w lp, "" u 1:3 w l lt 3

reset
#set term wxt 13
set term qt 13
unset key
set title "Bichromatically driven 87Rb\nDressed basis (SP)"
set xlabel "Frequency"
plot "Rb87_dressedoscillation_DRIVER.dat" u 1:2 w lp, "" u 1:3 w l lt 3


