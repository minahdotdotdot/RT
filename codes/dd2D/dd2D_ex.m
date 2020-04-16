run dd2D
%% Set up Linear Operator
L = genL(ks, ms, Rrho, Sc, tau, Nx, Nz);

A_RK4 = [[0; .5; 0; 0] [0; 0; .5; 0] [0; 0; 0; 1.0]];
b_RK4 = [1/6; 1/3; 1/3; 1/6];
c_RK4 = [0; 1/2; 1/2; 1]; 

T = 10;
h = tau/10;
M = T/h;
x = xIC;

for ii = 1 : M
	if sum(isnan(x)) == [0 0 0]
		x = exRK(x, A_RK4, b_RK4, c_RK4, h, L, Nx, Nz, dx, km);
	else
		break
	end
end