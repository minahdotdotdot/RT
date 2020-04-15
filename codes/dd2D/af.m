classdef af
	methods(Static)
		%% Reshape functions
		% These are for x = [Psi T S ]:
		function xvec = flatten(Psi, T, S)
			xvec = reshape([Psi T S]', 3*length(S), 1);
		end

		function [Psi, T, S] = boxify3(xvec)
			Psi = xvec(:,1);
			T = xvec(:,2);
			S = xvec(:,3);
		end

		% These are for each variable:
		function Q = boxify(Qvec,Nx,Nz)
			Q = reshape(Qvec, Nx, Nz)';
		end

		function Qvec = vectorize(Q, NxNz)
			Qvec = reshape(Q',NxNz,1);
		end

		%% Nonlinear tendency
		function ddx = partialx(f)
			ddx = [f(:,2)-f(:,end) f(:,3:end)-f(:, 1:end-2) f(:,1)-f(:,end-1)];
		end

		function ddy = partialy(f)
			ddy = [f(2,:)-f(end,:); f(3:end,:)-f(1:end-2,:); f(1,:)-f(end-1,:)];
		end

		function z = Jacobian(f,g,Nx,Nz,dx)
			f=ifft(ifft(boxify(f,Nx,Nz) ,Nz,1), Nx,2);
			g=ifft(ifft(boxify(g,Nx,Nz) ,Nz,1), Nx,2);
			z = (partialx(f).*partialy(g) - partialy(f).*partialx(g)) /(4*dx^2);
			z = fft(fft(z, Nx,2), Nz,1);
			z = reshape(z', Nz*Nx,1);
		end

		function nltend = NL(x, Nx, Nz, dx, km)
			[Psi, T, S] = boxify3(x);
			nltend = [Jacobian(Psi, km.*Psi, Nx, Nz, dx) ./km ...
			-Jacobian(Psi, T, Nx, Nz, dx)...
			-Jacobian(Psi, S, Nx, Nz, dx)];
		end

		%% Linear tendency
		% This creates a giant sparse operator L.
		function bigL = genL(ks, ms, Rrho, Sc, tau, Nx, Nz)
			bigL=sparse(3*Nx*Nz, 3*Nx*Nz);
			for jj = 1 : Nz
				for ii = 1 : Nx
					kk = (jj-1)*Nx+ii;
					bigL((kk-1)*3+1:3*kk, (kk-1)*3+1:3*kk)=gen3by3L(ks(ii), ms(jj), Rrho, Sc, tau);
				end
			end
		end

		% This is a 3-by-3 block of L_{k,m}.
		function lilL = gen3by3L(k, m, Rrho, Sc, tau)
			km = k^2+m^2;
			lilL = [[-Sc*km; -1i*k*Sc/(km*tau); 1i*k*Sc/(km*tau*Rrho); ] ...
			[-1i*k; -km/tau; 0] ... 
			[-1i*k; 0; -km/tau]];
		end

		%% IMEX operator
		% This creates a giant sparse operator (I-hL)^{-1}.
		function bigImexL = genImexInv(ks, ms, Rrho, Sc, tau, h, Nx, Nz)
			bigImexL=sparse(3*Nx*Nz, 3*Nx*Nz);
			for jj = 1 : Nz
				for ii = 1 : Nx
					kk = (jj-1)*Nx+ii;
					bigImexL((kk-1)*3+1:3*kk, (kk-1)*3+1:3*kk)=gen3by3ImexInv(ks(ii), ms(jj), Rrho, Sc, tau, h);
				end
			end
		end

		% This is the analytic solution to the 3-by-3 block of inv(I-hL_{k,m}).
		function lilImexL = gen3by3ImexInv(k, m, Rrho, Sc, tau, h)
			km = k^2+m^2;
			a = 1 + h*Sc*km;
			b = 1i*h*k;
			c = 1 + h*km/tau;
			d = 1i*h*k*Sc/(tau*km);
			lilImexL = [[c^2*Rrho; -b*c*Rrho; -b*c*Rrho] ...
			[-c*d*Rrho; b*d + a*c*Rrho; b*d*Rrho] ... 
			[c*d; -b*d; (a*c-b*d)*Rrho]] / (b*c*d+a*c^2*Rrho-b*c*d*Rrho);
		end

	end
end