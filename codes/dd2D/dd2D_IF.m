run dd2D
%Time-stepping 
T = 2000;
h = 5e-2;
M = T/h;
x = xIC;
every = 100;

%linear operator
L = genL(ks, ms, Rrho, Sc, tau, Nx, Nz);
[cx, cmat] = IFRK(xIC, M, h, L, Nx, Nz, ks, ms, every, "RK3", Rrho, Sc, tau);
