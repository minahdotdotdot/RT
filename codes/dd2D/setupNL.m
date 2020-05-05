run dd2D

%% Modify wave numbers for dealiasing.
ks = (2*pi/dp.Lx)*[0:.75*dp.Nx -.75*dp.Nx+1:-1];
ms = (2*pi/dp.Lz)*[0:.75*dp.Nz -.75*dp.Nz+1:-1]';
[kk,mm] = meshgrid(ks,ms);
km = kk.^2 + mm.^2;

%% Zero out N/2 for odd derivatives.
ks(dp.Nx/2) = 0;
ms(dp.Nz/2) = 0;
[kk,mm] = meshgrid(ks,ms);

%% Repack into a struct, dpNL. 
% Note that everything is the same as dp except for
% km, kk, and mm. 
dpNL = struct('Nx', dp.Nx, 'Nz', dp.Nz, 'NxNz', dp.NxNz, ...
	'km', km, 'kk', kk, 'mm', mm, ...%'x',dp.x,'z',dp.z,...
	'Lx', dp.Lx, 'Lz', dp.Lz);
vars = {'dp', 'ks', 'ms', 'kk', 'mm'};
clear(vars{:});clear vars;