%% Nonlinear tendency
function nltend = NL(x, Nx, Nz, dx, km, kk,mm)
	[Psi, T, S] = boxify3NL(x);
	nltend = [Jacobian(Psi, -km.*Psi, Nx, Nz, dx) ./km ...
	-Jacobian(Psi, T, Nx, Nz, dx)...
	-Jacobian(Psi, S, Nx, Nz, dx)];
	%Because we divide by 0 for km(1,1)=0
	nltend(1,1)=0;
end

function z = Jacobian(f,g,Nx,Nz,dx,kk,mm)
	dfdx=real(ifft2(1i*kk.*boxify(f,Nx,Nz)));
	dfdy=real(ifft2(1i*mm.*boxify(f,Nx,Nz)));
	dgdx=real(ifft2(1i*kk.*boxify(g,Nx,Nz)));
	dgdy=real(ifft2(1i*mm.*boxify(g,Nx,Nz)));
	%f=real(ifft2(boxify(f,Nx,Nz)));
	%g=real(ifft2(boxify(g,Nx,Nz)));
	%z = (partialx(f).*partialy(g) - partialy(f).*partialx(g)) /(4*dx^2);
	z = dfdx.*dgdy-dfdy.*dgdx;
	z = vectorize(fft2(z/(4*dx^2)),Nx*Nz);
end

function ddx = partialx(f)
	ddx = [f(:,2)-f(:,end) f(:,3:end)-f(:, 1:end-2) f(:,1)-f(:,end-1)];
end

function ddy = partialy(f)
	ddy = [f(2,:)-f(end,:); f(3:end,:)-f(1:end-2,:); f(1,:)-f(end-1,:)];
end