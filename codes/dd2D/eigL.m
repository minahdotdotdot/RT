function eigs = eigL(ks, ms, Rrho, Sc, tau, Nx, Nz)
	eigs = [];
    for jj = 1 : Nz
        for ii = 1 : Nx
        	eigs[end+1,:] = eig(gen3by3L(ks(ii), ms(jj), Rrho, Sc, tau, Nx))';
        end
    end
end
% This is a 3-by-3 block of L_{k,m}.
function lilL = gen3by3L(k, m, Rrho, Sc, tau, Nx)
    km = k^2+m^2;
    if km == 0 || k == Nx/2
        lilL = zeros(3,3);
    elseif k == Nx/3
        lilL = [[-Sc*km; 0; 0]...
            [0; -km/tau; 0]...
            [0;0;-km]];
    else
        lilL = [[-Sc*km; -1i*k; -1i*k]...
        [-1i*k*Sc/(km*tau); -km/tau; 0]...
        [1i*k*Sc/(km*tau*Rrho);0;-km]];
    end
end
