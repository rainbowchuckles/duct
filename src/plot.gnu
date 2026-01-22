set terminal pngcairo size 600,600

set output '/data/rho_cl.png'
p "../out/cl.dat" u 1:3 w l, "../out/cl.dat" u 1:4 w l, "../out/cl.dat" u 1:5 w l 
set output '/data/tv_cl.png'
p "../out/cl.dat" u 1:6 w l
set output '/data/u_cl.png'
p "../out/cl.dat" u 1:7 w l
set output '/data/p_cl1.png'
p "../out/cl.dat" u 1:8 w l
set output '/data/t_cl.png'
p "../out/cl.dat" u 1:9 w l
