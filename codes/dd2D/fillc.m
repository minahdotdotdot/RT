function cmat= fillc(cx, h, ks, ms, Rrho, Sc, tau, Nx, Nz)
	uniquec = unique(cx);
	cmat = containers.Map('KeyType', 'double', 'ValueType', 'any');
	for ii = 1 : length(uniquec)
		ch = uniquec(ii)*h;
		xx = expL(ch, ks, ms, Rrho, Sc, tau, Nx, Nz);
        %cmat(uniquec(ii)) = 1%xx;
        %expm((uniquec(ii)*h)*L);
    end
    %{
	cmat = cell(m,n);
	for ii = 1 : m
		for jj = 1 : n
			cmat{ii,jj}= c(x(ii,jj));
		end
	end
	%}
end
