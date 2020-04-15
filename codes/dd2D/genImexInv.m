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
