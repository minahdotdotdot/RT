run dd2D
%% Set up Linear Operator
L = genL(pp,dp);
eigs = eigL(pp,dp);
T = 1.06;
h = 0.9/max(max(abs(eigs)))
M = T/h;
a=load('RK4.mat');
x = a.xf;
every = 50;
[xf,ES,FS] = RK4(x, M, h, every, L, dp);
save('RK4-continued.mat', 'xf');

%{
A_RK4 = [[0; .5; 0; 0] [0; 0; .5; 0] [0; 0; 0; 1.0]];
b_RK4 = [1/6; 1/3; 1/3; 1/6];
c_RK4 = [0; 1/2; 1/2; 1]; 
for ii = 1 : M
	if sum(isnan(x)) == [0 0 0]
		%x = exRK(x, A_RK4, b_RK4, c_RK4, h, L, Nx, Nz, dx, km);
		x = RK4(x, h, L, Nx, Nz, dx, km);
	else
		break
	end
end
%}
