function xvec = flatten(Psi, T, S)
	xvec = reshape([Psi T S]', 3*length(S), 1);
end