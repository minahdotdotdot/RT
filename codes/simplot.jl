include("setupMMT.jl")
if scheme  ∈ ["IFRK3R", "IFRK4R"]
	using MAT
	file = matopen("../data/"*scheme*"h="*string(h)*"d"*string(deg)*"R.mat");
	crat = read(file, "crat");
	close(file)
	if scheme == "IFRK3R"
		RKT = eRKTableau(IFRK3.A, IFRK3.b, crat, IFRK3.x);
	elseif scheme == "IFRK4R"
		RKT = eRKTableau(IFRK4.A, IFRK4.b, crat, IFRK4.x)
	end
	runMMT(RKT, M, every, IC, h, L, NLfunc, fP, k, name)
else 
	runMMT(scheme, M, every, IC, h, L, NLfunc, fP, k, name)
end

using PyPlot, LaTeXStrings
saveEnergy!(k, N, T, name, scheme=scheme, h=h, ES=true)
saveEnergy!(k, N, T, name, scheme=scheme, h=h, ES=false)

