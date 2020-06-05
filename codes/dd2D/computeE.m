function [ES, FS] = computeE(x, dp)
	[Psi, T, S] = boxify3NL(x);
    %Sbox =real(ifft2(boxify(S, dp.Nx, dp.Nz)));
    %Psi_xbox = real(ifft2(1i*dp.kkk.*boxify(Psi, dp.Nx, dp.Nz)));
	ES = .5*sum(S.*conj(S))/(dp.NxNz^2);
	FS = 1i/2*sum(dp.kkk.* boxify(Psi.*conj(S)-conj(Psi).*S,dp.Nx, dp.Nz),'all')/(dp.NxNz^2);
    %ES =.5*sum(Sbox.^2,'all')/dp.NxNz;
    %FS = .5*sum(Psi_xbox.*Sbox,'all')/dp.NxNz;
end
