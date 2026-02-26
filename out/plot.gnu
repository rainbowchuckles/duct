set terminal pngcairo size 600,600
set output '/data/pls.png'
#p "cl.dat" u 1:17 title "T_T", "cl.dat" u 1:14 title "T_V", "cl.dat" u 1:16 title "p", "cl.dat" u 1:15 title "u"
set xlabel "Time, s"
set ylabel "Temperature, K"
set grid
#p "outlet.dat" u 2:12 w l title "T_T", "outlet.dat" u 2:9 w l title "T_V"
p "outlet.dat" u 2:11 w l title "T_T"

#p "cl_mole.dat" u 1:3, "cl_mole.dat" u 1:4, "cl_mole.dat" u 1:5, "cl_mole.dat" u 1:6, "cl_mole.dat" u 1:7, "cl_mole.dat" u 1:8, "cl_mole.dat" u 1:9, "cl_mole.dat" u 1:10, "cl_mole.dat" u 1:11, "cl_mole.dat" u 1:12, "cl_mole.dat" u 1:13 
