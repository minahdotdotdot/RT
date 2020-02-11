include("setupMMT.jl")
if scheme  ∈ ["IFRK3_rat", "IFRK4_rat"]
	using MAT
	file = matopen("../data/"*scheme*"h="*string(h)*"R.mat");
	crat = read(file, "crat");
	close(file)
	if scheme == "IFRK3_rat"
		RKT = eRKTableau(IFRK3.A, IFRK3.b, crat, IFRK3.x);
	elseif scheme == "IFRK4_rat"
		RKT = eRKTableau(IFRK4.A, IFRK4.b, crat, IFRK4.x)
	end
	runMMT(RKT, M, every, IC, h, L, NLfunc, fP, k, name)
else 
	runMMT(scheme, M, every, IC, h, L, NLfunc, fP, k, name)
end

using PyPlot, LaTeXStrings
saveEnergy!(k, N, T, name, scheme=scheme, h=h, ES=true)
saveEnergy!(k, N, T, name, scheme=scheme, h=h, ES=false)

