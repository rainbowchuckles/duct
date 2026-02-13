set terminal pngcairo size 600,600
set output '/data/pls.png'
#p "cl.dat" u 1:12 title "T_T", "cl.dat" u 1:9 title "T_V", "cl.dat" u 1:11 title "p", "cl.dat" u 1:10 title "u"
p "cl.dat" u 1:13, "lasta.dat" u 1:20
#p "cl_mole.dat" u 1:3, "cl_mole.dat" u 1:4, "cl_mole.dat" u 1:5, "cl_mole.dat" u 1:6, "cl_mole.dat" u 1:7, "cl_mole.dat" u 1:8, "cl_mole.dat" u 1:9, "cl_mole.dat" u 1:10, "cl_mole.dat" u 1:11, "cl_mole.dat" u 1:12, "cl_mole.dat" u 1:13 
