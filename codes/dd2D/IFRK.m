function [cx, cmat] = IFRK(x, M, h, L, Nx, Nz, ks, ms, every, name, Rrho, Sc, tau)
	[kk,mm] = meshgrid(ks,ms);
	km = kk.^2 + mm.^2; km = reshape(km', Nx*Nz, 1);
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
	cmat = fillc(cx, h, L);
	%{
	for tt = 1 : M
		x = IFRK_step(x, A, b, h, cx, cmat, Nx, Nz, km, kk, mm);
		if mod(tt, every) == 1
			%display([norm(x(:,1)) norm(x(:,2)) norm(x(:,3))])
			if ismember(1, isnan(x)) || ismember(1, isinf(x))
				%display(tt*h) 
				break
			end
		end
	end
	%}
end

function update = IFRK_step(x, A, b, h, cx, cmat, Nx, Nz, km, kk, mm)
	s = length(b);
	ks = cell(s,1); PP=0;
	% ks, x, PP are stored as vectors (flattened) within this function.
	ks{1} = flatten(NL(x, Nx, Nz, km, kk, mm));
	x = flatten(x);
	for ii = 2 : s
		PP = h*lincomIF(A(ii, 1:ii-1), cx(ii, 1:ii-1), cmat, ks(1:ii-1));
		ks{ii} = flatten(NL(boxify3(cmat{ii,ii}*x + PP, Nx, Nz), Nx, Nz, km, kk, mm));
	end
	%new PP for update
	PP = h*lincomIF(b, cx(end,:), cmat, ks);
	update = boxify3(cmat(cx(end-1,end))*x+PP, Nx, Nz);
end

%lincom outputs a flat vector
function lincom = lincomIF(A, cx, cmat, k)
	lincom = A(1)*cmat(cx(1))*k{1};
	for ii = 2 : length(k)
		lincom = lincom + A(ii)*cmat(cx(ii))*k{ii};
	end
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
