%% Nonlinear tendency
function nltend = NL(x, Nx, Nz, dx, km)
	[Psi, T, S] = boxify3NL(x);
	nltend = [Jacobian(Psi, -km.*Psi, Nx, Nz, dx) ./km ...
	-Jacobian(Psi, T, Nx, Nz, dx)...
	-Jacobian(Psi, S, Nx, Nz, dx)];
	%Because we divide by 0 for km(1,1)=0
	nltend(1,1)=0;
end

function z = Jacobian(f,g,Nx,Nz,dx)
	f=real(ifft2(boxify(f,Nx,Nz)));
	g=real(ifft2(boxify(g,Nx,Nz)));
	z = (partialx(f).*partialy(g) - partialy(f).*partialx(g)) /(4*dx^2);
	z = vectorize(fft2(z));
end

function ddx = partialx(f)
	ddx = [f(:,2)-f(:,end) f(:,3:end)-f(:, 1:end-2) f(:,1)-f(:,end-1)];
end

function ddy = partialy(f)
	ddy = [f(2,:)-f(end,:); f(3:end,:)-f(1:end-2,:); f(1,:)-f(end-1,:)];
end