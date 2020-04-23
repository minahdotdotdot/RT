run dd2D
%Time-stepping 
T = 10;
h = 5e-2;
M = T/h;
x = xIC;
every = 100;

%parpool
%linear operator
workers = 8;
L = genL(ks, ms, Rrho, Sc, tau, Nx, Nz);
xf = IFRK(xIC, M, h, Nx, Nz, ks, ms, km, every, "RK3", Rrho, Sc, tau, L, workers);
