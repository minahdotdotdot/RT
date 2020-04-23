function eigs = eigL(pp, dp)
	eigs = [];
    for jj = 1 : dp.Nz
        for ii = 1 : dp.Nx
        	eigs(end+1,:) = eig(gen3by3L(dp.ks(ii), dp.ms(jj), pp, dp))';
        end
    end
end
% This is a 3-by-3 block of L_{k,m}.
function lilL = gen3by3L(k, m, pp, dp)
    %Rrho, Sc, tau, Nx)
    km = k^2+m^2;
    if km == 0
        lilL = zeros(3,3);
    elseif k == dp.Nx/2
        lilL = [[-pp.Sc*km; 0; 0]...
            [0; -km/pp.tau; 0]...
            [0;0;-km]];
    else
        lilL = [[-pp.Sc*km; -1i*k; -1i*k]...
        [-1i*k*pp.Sc/(km*pp.tau); -km/pp.tau; 0]...
        [1i*k*pp.Sc/(km*pp.tau*pp.Rrho);0;-km]];
    end
end
