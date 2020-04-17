% By default, every variable is in Fourier space.
% x = [Psi T S], where Psi, T, S each are columns.
% This is the default shape of our variables.
% The linear solve requires x to be flattend via (i,j) index
%[Psi(1,1) T(1,1) S(1,1) ... Psi(Nx,1) T(Nx,1) S(Nx,1)
%... Psi(1,2) T(1,2) S(1,2) ... Psi(Nx,2) T(Nx,2) S(Nx,2)
%... Psi(1,Nz) T(1,Nz) S(1,Nz) ... Psi(Nx,Nz) T(Nx,Nz) S(Nx,Nz)]'
% i.e. boxPsi = reshape(Psi, Nx, Nz)' should yield the exact physical 
%      domain with top left as the origin. 
%      reshape(boxPsi', Nz*Nx, 1) returns boxT to flattned shape.
% Function boxify assigns [Psi, T, S] = boxify(x, ...)
% Function flatten reshapes x to be appropriate for linear solve. 

%% Problem Parameters 
tau = 0.01;
Pr = 7;
Ra = 1.1;
Sc = Pr/tau;
Rrho = 1/(Ra*tau);

%% Domain
a_ratio = 2;
N = 2^7; Nx = N; Nz = a_ratio*N; NxNz = Nx*Nz;
l_o = 2*pi/( .25*(-2-Ra + Ra*sqrt(1+8/Ra)) )^(.25);
Lx = l_o; Lz = a_ratio*Lx;
x = linspace(0, Lx, N+1); x=x(1:end-1);
z = linspace(0, Lz, a_ratio*N+1); z=z(1:end-1);

%% Discretization parameters
dx = Lx/Nx;
ks = (2*pi/Lx)*[linspace(0, Nx/2, Nx/2+1)';linspace(-Nx/2+1,-1,Nx/2-1)'];
ms = (2*pi/Lz)*[linspace(0, Nz/2, Nz/2+1)';linspace(-Nz/2+1,-1,Nz/2-1)'];
[kk,mm] = meshgrid(ks,ms);
km = kk.^2 + mm.^2; km = reshape(km', NxNz, 1);

%% Initial Condition
xIC = fft2(randn(NxNz,3)/NxNz);
%[Psi, T, S] = boxify3NL(xIC);
%Psibox = boxify(Psi, Nx, Nz);
