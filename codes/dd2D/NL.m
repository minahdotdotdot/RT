%% Nonlinear tendency
function nltend = NL(x, dp)
	%Unpack into box domain form.
	[Psi, T, S] = boxify3NL(x);

	%Pad variables with zeros
	pPsi = pad_var(Psi, dp);
	pT   = pad_var(T, dp);
	pS   = pad_var(S, dp);

	%Form omega = Laplacian(Psi)
	pomega= -dp.km.*pPsi;

	%Velocities
    u = real(ifft2(-1i*dp.mm.*pPsi));
    w = real(ifft2(1i*dp.kk.*pPsi));
    
	nltend = [...
	unpad_var(Jacobian(u, w, pomega, dp) ./dp.km, dp) ...
	unpad_var(-Jacobian(u, w, pT, dp), dp)...
	unpad_var(-Jacobian(u, w, pS, dp), dp)...
	];
	
	%Because we divide by 0 for km(1,1)=0
	nltend(1,1)=0;
end

function pq =pad_var(q, dp)
	q = boxify(q, dp.Nx, dp.Nz);
	Nx = dp.Nx; Nz = dp.Nz;
	% Pad q.
	% I.e., if you set N=256 then the code uses 3*N/2=384 Fourier modes (in each
    % direction) to compute the Jacobian.
    % physical space, 3/2 grid; factor of (9/4) scales fft
    pq = zeros([1.5*dp.Nz 1.5*dp.Nx]);
    pq(1    : Nz/2+1, 1    : Nx/2+1) = (9/4)*q(1      : Nz/2+1, 1      : Nx/2+1);
    pq(1    : Nz/2+1, Nx+2 : 1.5*Nx) = (9/4)*q(1      : Nz/2+1, Nx/2+2 : Nx    );
    pq(Nz+2 : 1.5*Nz, 1    : Nx/2+1) = (9/4)*q(Nz/2+2 : Nz,     1      : Nx/2+1);
    pq(Nz+2 : 1.5*Nz, Nx+2 : 1.5*Nx) = (9/4)*q(Nz/2+2 : Nz,     Nx/2+2 : Nx    ); 
end

function z = Jacobian(u,w,pq,dp)
	xq=real(ifft2(1i*dp.kk.*pq));
	zq=real(ifft2(1i*dp.mm.*pq));
	z =fft2(u.*xq + w.*zq);
end


%unpad_var unpads AND vectorizes
function q = unpad_var(pq, dp)
	Nx = dp.Nx; Nz = dp.Nz;
	q = zeros([Nz Nx]);
    q(1      : Nz/2+1, 1      : Nx/2+1) = pq(1    : Nz/2+1, 1    : Nx/2+1);
    q(1      : Nz/2+1, Nx/2+2 : Nx)     = pq(1    : Nz/2+1, Nx+2 : 1.5*Nx);
    q(Nz/2+2 : Nz,     1      : Nx/2+1) = pq(Nz+2 : 1.5*Nz, 1    : Nx/2+1);
    q(Nz/2+2 : Nz,     Nx/2+2 : Nx)     = pq(Nz+2 : 1.5*Nz, Nx+2 : 1.5*Nx);
	q = vectorize(q, dp.NxNz);
end


%{
function ddx = partialx(f)
	ddx = [f(:,2)-f(:,end) f(:,3:end)-f(:, 1:end-2) f(:,1)-f(:,end-1)];
end

function ddy = partialy(f)
	ddy = [f(2,:)-f(end,:); f(3:end,:)-f(1:end-2,:); f(1,:)-f(end-1,:)];
end
%}