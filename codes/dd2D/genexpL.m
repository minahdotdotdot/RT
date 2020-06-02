function expL= genexpL(L, workers, bs) % Block size has to be divisible by 3 or it will cause errors
    [m,n] = size(L); %m should be 3*Nx*Nz
    if mod(m, workers*bs) == 0
        wp = m / bs; % Number of blocks to be exponentiated.  
        blocks = cell(wp,1);
        parfor ii = 1 : wp
            [spi, spj, snz] = find(expm(L((ii-1)*bs+1: ii*bs, (ii-1)*bs+1:ii*bs)));
            blocks{ii}=[(ii-1)*bs+spi (ii-1)*bs+spj snz];
        end
        SP = zeros(9*m,3); %[rowind, colind, nzval]
        iold=1; inew=1;
        for ii = 1 : wp
            inew = iold-1 + size(blocks{ii},1);
            SP(iold:inew,:)=blocks{ii}(:,:);
            iold=inew+1;
        end
        SP = SP(1:iold-1,:);
        expL = sparse(SP(:,1), SP(:,2), SP(:,3));
    else
        error('check m/workers*bs')
    end
   
end
