function [ES, FS] = computeE(x, dp)
	[Psi, T, S] = boxify3NL(x);
	ES = .5*norm(S)^2/dp.NxNz;
	FS = 1i/2*sum(sum(dp.kk.* boxify(Psi.*conj(S)-conj(Psi).*S,dp.Nx, dp.Nz)))/dp.NxNz;
end