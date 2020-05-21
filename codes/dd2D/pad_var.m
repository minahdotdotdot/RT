function pq =pad_var(q, dp)
	q = boxify(q, dp.Nx, dp.Nz);
	Nx = dp.Nx; Nz = dp.Nz;
	% Pad q.
	% I.e., if you set N=256 then the code uses 3*N/2=384 Fourier modes (in each
    % direction) to compute the Jacobian.
    % physical space, 3/2 grid; factor of (9/4) scales fft
    pq = zeros([1.5*Nz 1.5*Nx]);
    pq(1    : Nz/2+1, 1    : Nx/2+1) = (9/4)*q(1      : Nz/2+1, 1      : Nx/2+1);
    pq(1    : Nz/2+1, Nx+2 : 1.5*Nx) = (9/4)*q(1      : Nz/2+1, Nx/2+2 : Nx    );
    pq(Nz+2 : 1.5*Nz, 1    : Nx/2+1) = (9/4)*q(Nz/2+2 : Nz,     1      : Nx/2+1);
    pq(Nz+2 : 1.5*Nz, Nx+2 : 1.5*Nx) = (9/4)*q(Nz/2+2 : Nz,     Nx/2+2 : Nx    ); 
end