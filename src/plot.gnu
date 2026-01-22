set terminal pngcairo size 600,600
set output '/data/3d.png'
sp "../out/flowfield.dat" u 1:2:3

set output '/data/rho_outlet.png'
p "../out/outlet.dat" u 2:3 w l
set output '/data/u_outlet.png'
p "../out/outlet.dat" u 2:4 w l
set output '/data/p_outlet.png'
p "../out/outlet.dat" u 2:5 w l

set output '/data/rho_cl.png'
p "../out/cl.dat" u 1:3 w l 
set output '/data/u_cl.png'
p "../out/cl.dat" u 1:4 w l
set output '/data/p_cl.png'
p "../out/cl.dat" u 1:5 w l
