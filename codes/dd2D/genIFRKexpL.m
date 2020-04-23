function RK = genIFRKexpL(name, h, L, workers, bs, pp, dp)
	A=0;b=0;cx=0;
	if name == "RK3"
		A = [[0; .5; -1],[0; 0; 2]];
		b = [1/6; 2/3; 1/6];
		cx = [0; 1/2; 1]; 
	elseif name == "RK4"
		A = [[0; .5; 0; 0],...
		[0; 0; .5; 0],...
		[0; 0; 0; 1.0]];
		b = [1/6; 1/3; 1/3; 1/6];
		cx = [0; 1/2; 1/2; 1]; 
	end
	cx = cfromx(cx);
	%delete(gcp('nocreate'))
    %parpool(workers)
	%tic
	cmat = fillc(cx, h, L, workers, bs);
	%toc
	RK = struct('A', A, 'b', b, 'cx', cx, 'cmat', cmat);
end

function c = cfromx(cx)
	s = length(cx);
	c = zeros(s+1, s);
	for ii = 1 : s
		for jj = 1 : ii-1
			c(ii,jj) = cx(ii)-cx(jj);
		end
		c(ii,ii) = cx(ii);
	end
	c(end,1:s) = 1 - cx;
end

function cmat= fillc(cx, h, L, workers, bs)
	uniquec = unique(cx);
	cmat = containers.Map('KeyType', 'double', 'ValueType', 'any');
	for ii = 1 : length(uniquec)
		cmat(uniquec(ii)) = genexpL((uniquec(ii)*h)*L, workers, bs);
		%genexpL(uniquec(ii)*h, ks, ms, Rrho, Sc, tau, Nx, Nz);
    end
end