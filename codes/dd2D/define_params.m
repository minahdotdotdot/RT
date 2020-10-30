function [dp,pp] = define_params(N, D, Pr, a_ratio)
	%% Problem Parameters 
	tau = 0.01;
	%Pr = 1; %water: 7
	Ra = 1.1;
	Sc = Pr/tau;
	Rrho = 1/(Ra*tau); 
	pp = struct('tau', tau, 'Pr', Pr, 'Ra', Ra, 'Sc', Sc, 'Rrho', Rrho);

	%% Domain
	a_ratio = 1;
	Nx = N; Nz = a_ratio*N; NxNz = Nx*Nz;
	k_o = ( .25*(-2-Ra + Ra*sqrt(1+8/Ra)) )^(.25);     % WON'T NEED
	l_o = 2*pi/k_o;
	Lx = D*l_o; Lz = a_ratio*Lx;
	x = linspace(0, Lx, Nx+1); x=x(2:end);
	z = linspace(0, Lz, Nz+1); z=z(2:end);

	%% DEALIASED WAVE NUMBERS
	ks = (2*pi/Lx)*[0:.75*Nx -.75*Nx+1:-1];            % WON'T NEED
	ms = (2*pi/Lz)*[0:.75*Nz -.75*Nz+1:-1]';           % WON'T NEED
	[kk,mm] = meshgrid(ks,ms);
	km = kk.^2 + mm.^2;

	%% Zero out N/2 for odd derivatives.
	ks(Nx/2) = 0;                                   % WON'T NEED
	ms(Nz/2) = 0;                                   % WON'T NEED
	[kk,mm] = meshgrid(ks,ms);

	%% ALIASED WAVE NUMBERS
	dx = Lx/Nx;
	ks = (2*pi/Lx)*[0:Nx/2 -Nx/2+1:-1]';               % WON'T NEED
	ms = (2*pi/Lz)*[0:Nz/2 -Nz/2+1:-1]';               % WON'T NEED
	% Nx/2,Nz/2 are not zero since these are for setting up Laplacian.
	[kkk,mmm] = meshgrid(ks,ms);
	%km = kkk.^2 + mmm.^2; km = reshape(km', NxNz, 1); % WON'T NEED

	dp = struct('Nx', Nx, 'Nz', Nz, 'NxNz', NxNz, ...
		'ks', ks, 'ms', ms, ... %to write the linear operator. 
		'km', km, 'kk', kk, 'mm', mm, 'kkk', kkk, 'mmm', mmm, 'x', x,'z', z,...
		'Lx', Lx, 'Lz', Lz, 'l_o', l_o);
end