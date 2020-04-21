run dd2D
%Time-stepping 
T = 10;
h = 5e-7;
M = 5000;%T/h;
x = xIC;
every = 1000;

%linear operator
L = genL(ks, ms, Rrho, Sc, tau, Nx, Nz);
xf = IFRK(xIC, M, h, L, Nx, Nz, km, kk, mm, every, "RK3");