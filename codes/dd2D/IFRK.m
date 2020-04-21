function x = IFRK(IC, M, h, bigL, Nx, Nz, km, kk, mm, every, name)
	if name == "RK3"
		A = [[0; .5; -1]...[0; 0; 2]]
		b = [1/6; 2/3; 1/6];
		x = [0; 1/2; 1]; 
	elseif name == "RK4"
		A = [[0; .5; 0; 0]...
		[0; 0; .5; 0]...
		[0; 0; 0; 1.0]];
		b = [1/6; 1/3; 1/3; 1/6];
		x = [0; 1/2; 1/2; 1]; 
	end
	x = cfromx(x);
	c = fillc(x,h,L);
	cmat = filc(c, h, L);
	for tt = 1 : M
		x = IFRK4_step(x, h, bigL, Nx, Nz, km, kk, mm, A, b, cmat);
		if mod(tt, every) == 1
            display(tt*h) 
			if ismember(1, isnan(x)) || ismember(1, isinf(x))
				break
			end
		end
	end
end

function update = IFRK_step(x, A, b, cmat, Nx, Nz, km, kk, mm)
	s = length(b);
	ks = cell(b,1); PP=0;
	% ks, x, PP are stored as vectors (flattened) within this function.
	x = flatten(x);
	ks{1} = flatten(NL(x, Nx, Nz, km, kk, mm));
	for ii = 2 : x
		PP = h*lincomIF(A(ii, 1:ii-1), cmat(ii, 1:ii-1), ks(1:ii-1);
		ks{ii} = flatten(NL(boxify3(cmat{ii,ii}*x + PP, Nx, Nz), Nx, Nz, km, kk, mm));
	end
	%new PP for update
	PP = h*lincomIF(b, c(end,:), ks)
	update = boxify3(c{end-1,end}*x+PP, Nx, Nz);
end

%lincom outputs a flat vector
function lincom = lincomIF(A, cmat, k)
	lincom = A(1)*cmat{1}*k
	for ii = 2 : length(k)
		lincom = lincom + A(ii)*cmat{ii}*k(ii)
	end
end

function cmat = fillc(x, h, L)
	[m,n] = size(x)
	uniquec = unique(x)
	c = containers.Map('KeyType', 'double', 'ValueType', 'any');
	for ii = 1 : length(uniquec)
        c(uniquec(ii)) = expm((uniquec(ii)*h)*L);
    end
	cmat = cell(m,n)
	for ii = 1 : m
		for jj = 1 : n
			cmat{ii,jj}= c(x(ii,jj))
		end
	end
end

function c = cfromx(x)
	s = length(x);
	c = zeros(s+1, s);
	for ii = 1 : s
		for jj = 1 : i-1
			c(ii,jj) = x(ii)-x(jj);
		end
		c(ii,ii) = x(ii);
	end
	c(end,1:s) = 1 - x;
end
