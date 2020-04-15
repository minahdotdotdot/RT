%% Nonlinear tendency
function nltend = NL(x, Nx, Nz, dx, km)
	[Psi, T, S] = boxify3(x);
	nltend = [Jacobian(Psi, km.*Psi, Nx, Nz, dx) ./km ...
	-Jacobian(Psi, T, Nx, Nz, dx)...
	-Jacobian(Psi, S, Nx, Nz, dx)];
end

function z = Jacobian(f,g,Nx,Nz,dx)
	f=ifft(ifft(boxify(f,Nx,Nz) ,Nz,1), Nx,2);
	g=ifft(ifft(boxify(g,Nx,Nz) ,Nz,1), Nx,2);
	z = (partialx(f).*partialy(g) - partialy(f).*partialx(g)) /(4*dx^2);
	z = fft(fft(z, Nx,2), Nz,1);
	z = reshape(z', Nz*Nx,1);
end

function ddx = partialx(f)
	ddx = [f(:,2)-f(:,end) f(:,3:end)-f(:, 1:end-2) f(:,1)-f(:,end-1)];
end

function ddy = partialy(f)
	ddy = [f(2,:)-f(end,:); f(3:end,:)-f(1:end-2,:); f(1,:)-f(end-1,:)];
end