%% IMEX operator
% This creates a giant sparse operator (I-hdL)^{-1}.
function IL = genImexInv(hd, pp, dp, par)
	bigImexL=sparse(3*dp.NxNz, 3*dp.NxNz);
	if par == 0
		for jj = 1 : Nz
			for ii = 1 : Nx
				kk = (jj-1)*Nx+ii;
				IL((kk-1)*3+1:3*kk, (kk-1)*3+1:3*kk)=gen3by3ImexInv(hd, ks(ii), ms(jj), pp);
			end
		end
	else
		workers = par;
		if mod(dp.Nz, workers) == 0
            p = dp.Nz / workers;
            blocks = cell(workers,1);
            parfor w = 1 : workers
                blocks{w} = sparse(3*dp.Nx*p, 3*dp.Nx*p);
                jjj = 0;
                for jj = (w-1)*p+1:w*p
                    jjj = jjj + 1;
                    for ii = 1 : dp.Nx
                    kk = (jjj-1)*dp.Nx + ii;
                        blocks{w}((kk-1)*3+1:3*kk, (kk-1)*3+1:3*kk)=gen3by3ImexInv(hd, dp.ks(ii), dp.ms(jj), pp, dp);
                    end
                end
            end
            m = 3*dp.Nx*p;
            for w = 1 : workers
                IL((w-1)*m+1:w*m,(w-1)*m+1:w*m) = blocks{w};
            end
            if issparse(bigL) == false
                IL = sparse(IL);
            end
        else
            error('check m/workers')
        end

	end
end

% This is the analytic solution to the 3-by-3 block of inv(I-hL_{k,m}).
function iL = gen3by3ImexInv(hd, k, m, pp, dp)
	Sc = pp.Sc; tau=pp.tau; Rrho=pp.Rrho;
	km = k^2+m^2;
	if km == 0
		iL = zeros(3,3);
    else
    	if k == dp.Nx/2
    		k = 0;
    	end
		a = 1 + h*Sc*km;
		b = 1i*h*k;
		c = 1 + h*km/tau;
		d = 1i*h*k*Sc/(tau*km);
		e = 1 + h*km;
		iL = [[c*e*Rrho; -b*e*Rrho; -b*c*Rrho]...
			[-d*e*R; b*d+a*e*Rrho; b*d*Rrho]...
			[c*d; -b*d; (a*c-b*d)*Rrho]] / (a*c*e*Rrho - b*d*(c-e*Rrho));
	end
end
