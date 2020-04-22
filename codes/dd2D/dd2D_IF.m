run dd2D
%Time-stepping 
T = 2000;
h = 5e-2;
M = T/h;
x = xIC;
every = 100;

parpool
%linear operator
%L = genL(ks, ms, Rrho, Sc, tau, Nx, Nz);
xf = IFRK(xIC, M, h, L, Nx, Nz, ks, ms, km, every, "RK3", Rrho, Sc, tau);
