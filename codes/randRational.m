m =1000;
deg = 8;
A = randn(m,m);
%D = eig(A);
[Q,R]=qr(A); Q=Q*sparse(diag(sign(diag(R))));
D = randn(1,m) + j*randn(1,m) - 5*rand(1,m);
A = Q*sparse(diag(D))*Q';
Z = exp(D);
[r, pol, res, zer, zj, fj, wj] = aaa(Z,D,'degree', deg);
R = polzer(pol, zer, D, r);
Zm = expm(A);
Rmatrix = polzer(pol, zer, A, r);
%Zm = expm(diag(D));
%Rmatrix = polzer(pol, zer, diag(D), r);

norm(Zm-Rmatrix)

%R = polzer(pol, zer, diag(D), r);
%R = polzer(pol, zer, A, r);
%norm(Z-R)
function [R]=polzer(pol, zer, z, r)
	bary0 = r(0);                    % additive shift
	gam = bary0*prod(pol)/prod(zer); % scaling
	R = gam;
	if isvector(z)==isvector([1])
		for i = 1 : length(zer)
			R = R .* (z-zer(i));
		end
		for i = 1 : length(pol)
			R = R ./ (z-pol(i));
		end
		R = R + 1 - bary0;
	elseif isscalar(z)== isscalar(1)
		for i = 1 : length(zer)
			R = R*(z-zer(i))/(z-pol(i));
		end
		R = R + 1 - bary0;
	elseif ismatrix(z)==ismatrix(randn(1,1))
		I=eye(size(z));
		P = z - pol(1)*I;
		Z = R * (z - zer(1)*I); %gam included
		for i = 2 : length(zer)
			Z = Z*(z - zer(i)*I);
		end
		for i = 2 : length(pol)
			P = P*(z - pol(i)*I);
		end
		R = inv(P)*Z + (1-bary0)*I;
	end
end