%% Nonlinear tendency
function nltend = NL(x, Nx, Nz, km, kk, mm)
	nltend = zeros(Nx*Nz, 3);
	%[Psi, T, S] = boxify3NL(x);
	%nltend = [Jacobian(Psi, -km.*Psi, Nx, Nz, kk, mm) ./km ...
	%-Jacobian(Psi, T, Nx, Nz, kk, mm)...
	%-Jacobian(Psi, S, Nx, Nz, kk, mm)];
	%Because we divide by 0 for km(1,1)=0
	%nltend(1,1)=0;
end

function z = Jacobian(f,g,Nx,Nz,kk,mm)
	dfdx=real(ifft2(1i*kk.*boxify(f,Nx,Nz)));
	dfdy=real(ifft2(1i*mm.*boxify(f,Nx,Nz)));
	dgdx=real(ifft2(1i*kk.*boxify(g,Nx,Nz)));
	dgdy=real(ifft2(1i*mm.*boxify(g,Nx,Nz)));
	%f=real(ifft2(boxify(f,Nx,Nz)));
	%g=real(ifft2(boxify(g,Nx,Nz)));
	%z = (partialx(f).*partialy(g) - partialy(f).*partialx(g)) /(4*dx^2);
	z = dfdx.*dgdy-dfdy.*dgdx;
	z = vectorize(fft2(z),Nx*Nz);
end

function ddx = partialx(f)
	ddx = [f(:,2)-f(:,end) f(:,3:end)-f(:, 1:end-2) f(:,1)-f(:,end-1)];
end

function ddy = partialy(f)
	ddy = [f(2,:)-f(end,:); f(3:end,:)-f(1:end-2,:); f(1,:)-f(end-1,:)];
end
