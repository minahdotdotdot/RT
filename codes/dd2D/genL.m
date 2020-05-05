%% Linear tendency
% This creates a giant sparse operator L.

function bigL = genL(pp, dp, par)
    bigL=sparse(3*dp.Nx*dp.Nz, 3*dp.Nx*dp.Nz);
    if par == 0
        for jj = 1 : dp.Nz
            for ii = 1 : dp.Nx
                kk = (jj-1)*dp.Nx+ii;
                bigL((kk-1)*3+1:3*kk, (kk-1)*3+1:3*kk)=gen3by3L(dp.ks(ii), dp.ms(jj), pp, dp);
                %(ks(ii), ms(jj), Rrho, Sc, tau, Nx);
            end
        end
    else
        delete(gcp('nocreate'))
        parpool(workers)
        workers =par(1); m=3*dp.Nx*dp*Nz;
        if mod(m, workers) == 0
            p = m / workers;
            blocks = cell(workers,1);
            parfor ii = 1 : workers
                for jj = 1 : p
                    kk = (jj-1)*workers + ii;
                    if jj == 1
                        blocks{ii} = cell(p,1);
                    end
                    yy = floor(kk/Nx); xx = kk-(yy-1)dp.Nx;
                    blocks{ii}{jj}=gen3by3L(dp.ks(xx), dp.ms(yy), pp, dp);
                end
            end
            expL = sparse(m,n);
            for kk = 1 : workers*p
                jj= ceil(kk/workers); ii = kk - (jj-1)*workers;
                bigL((kk-1)*3+1:3*kk, (kk-1)*3+1:3*kk) = blocks{ii}{jj};
            end
            if issparse(bigL) == false
                expL = sparse(bigL);
            end
        else
            error('check m/workers')
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
