function [ES, FS] = computeE(x, dp)
	[Psi, T, S] = boxify3NL(x);
	ES = .5*sum(S.^2)/(dp.NxNz^2);
	FS = 1i/2*sum(sum(dp.kkk.* boxify(Psi.*conj(S)-conj(Psi).*S,dp.Nx, dp.Nz)))/(dp.NxNz^2);
end
