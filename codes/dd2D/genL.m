%% Linear tendency
% This creates a giant sparse operator L.

function bigL = genL(pp, dp)%ks, ms, Rrho, Sc, tau, Nx, Nz)
    bigL=sparse(3*dp.Nx*dp.Nz, 3*dp.Nx*dp.Nz);
    for jj = 1 : dp.Nz
        for ii = 1 : dp.Nx
            kk = (jj-1)*dp.Nx+ii;
            bigL((kk-1)*3+1:3*kk, (kk-1)*3+1:3*kk)=gen3by3L(dp.ks(ii), dp.ms(jj), pp, dp);
            %(ks(ii), ms(jj), Rrho, Sc, tau, Nx);
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
