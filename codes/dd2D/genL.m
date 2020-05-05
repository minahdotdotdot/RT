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
        workers = par;
        %delete(gcp('nocreate'))
        %parpool(workers)
        if mod(dp.Nz, workers) == 0
            p = dp.Nz / workers;
            blocks = cell(workers,1);
            parfor w = 1 : workers
                blocks{w} = sparse(3*dp.Nx*p, 3*dp.Nx*p);
                jjj = 0;
                for jj = (w-1)*p+1:w*p
                    jjj = jjj + 1;
                    for ii = 1 : dp.Nx
                    kk = (jjj-1)*workers + ii;
                        blocks{w}((kk-1)*3+1:3*kk, (kk-1)*3+1:3*kk)=gen3by3L(dp.ks(ii), dp.ms(jj), pp, dp);
                    end
                end
            end
            m = 3*dp.Nx*p;
            for w = 1 : workers
                bigL((w-1)*m+1:w*m,(w-1)*m+1:w*m) = blocks{w};
            end
            if issparse(bigL) == false
                bigL = sparse(bigL);
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
