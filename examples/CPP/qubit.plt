set pm3d at b map
unset surf
set xlabel "Frequency"
set ylabel "Time"
splot "qubit_bareoscillation_SP.dat" u 1:2:3
