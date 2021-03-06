function bigexpL = expL(ch, ks, ms, Rrho, Sc, tau, Nx, Nz)
    if ch == 0
        bigexpL = speye(3*Nx*Nz);
    else
        bigexpL=sparse(3*Nx*Nz, 3*Nx*Nz);
        for kk = 1 : Nx*Nz
        	jj = ceil(kk/Nx); 
        	ii = kk - (jj-1)*Nx;
        	bigexpL((kk-1)*3+1:3*kk, (kk-1)*3+1:3*kk)=gen3by3expL(ch, ks(ii), ms(jj), Rrho, Sc, tau, Nx);
        end
    end
end
% This is a 3-by-3 block of L_{k,m}.
function lilexpL = gen3by3expL(ch, k, m, Rrho, Sc, tau, Nx)
    km = k^2+m^2;
    lilexpL = 0;
    if km == 0
        lilexpL = speye(3);
    elseif k == Nx/2
        lilexpL = sparse(expm(ch*[[-Sc*km; 0; 0],...
            [0; -km/tau; 0],...
            [0;0;-km]]));
    else
        lilexpL = sparse(expm(ch*[[-Sc*km; -1i*k; -1i*k]...
        [-1i*k*Sc/(km*tau); -km/tau; 0]...
        [1i*k*Sc/(km*tau*Rrho);0;-km]]));
    end
end
