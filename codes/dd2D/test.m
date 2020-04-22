tic
n = 200;
A = 500;
a = zeros(1,n);
parfor i = 1:n
    a(i) = max(abs(eig(rand(A))));
end
toc

tic
n = 200;
A = 500;
a = zeros(1,n);
for i = 1:n
    a(i) = max(abs(eig(rand(A))));
end
toc


%{
%test partialx and partialy
xx = linspace(0,2*pi,Nz+1)';xx=xx(1:end-1);
f = cos(xx)+1i*sin(xx);
ff = [f f];
ddy = partialy(f)/(2*dx);
zz=ddy(:,1) ./ (-sin(xx)+1i*cos(xx)) - (2*pi/Lz);
zz(1)

xx = linspace(0,2*pi,Nx+1)';xx=xx(1:end-1);
f = cos(xx)+1i*sin(xx);
ff = repmat(f', 2,1);
ddx = partialx(ff)/(2*dx);
zz=ddx(1,:)' ./ (-sin(xx)+1i*cos(xx)) - (2*pi/Lx);
zz(1)

nltend = NL(xIC, Nx, Nz, dx, km);
xvec = flatten(Psi, T, S);
z=Jacobian(Psi, T, Nx, Nz, dx);
nltend = NL(xIC, Nx, Nz, dx, km);

Psi = ones(5,1); T = 2*ones(5,1); S = 3*ones(5,1);
xvec = flatten(Psi, T, S);

x = linspace(1,15,15);
xbox = boxify(x,3,5);
x3 = flatten(ones(5,1), 2*ones(5,1), 3*ones(5,1));
%}


