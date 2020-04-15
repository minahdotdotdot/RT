%% Linear tendency
% This creates a giant sparse operator L.
function bigL = genL(ks, ms, Rrho, Sc, tau, Nx, Nz)
    bigL=sparse(3*Nx*Nz, 3*Nx*Nz);
    for jj = 1 : Nz
        for ii = 1 : Nx
            kk = (jj-1)*Nx+ii;
            if kk > 1
                bigL((kk-1)*3+1:3*kk, (kk-1)*3+1:3*kk)=gen3by3L(ks(ii), ms(jj), Rrho, Sc, tau);
            end
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
end
