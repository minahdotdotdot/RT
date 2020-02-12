include("setupMMT.jl")
if scheme  ∈ ["IFRK3R", "IFRK4R"]
	using MAT
	include("IF_methods.jl");
	file = matopen("../data/"*scheme*"h="*string(h)*"d"*string(deg)*"R.mat");
	crat = read(file, "crat");
	close(file)
	#norm(IFRK3.c-crat2)/norm(IFRK3.c)
	RKT=0
	if scheme == "IFRK3R"
		RKT = eRKTableau(IFRK3.A, IFRK3.b, crat, IFRK3.x);
	elseif scheme == "IFRK4R"
		RKT = eRKTableau(IFRK4.A, IFRK4.b, crat, IFRK4.x)
	end
	#print(norm(IFRK3.c-crat))
	runMMT(RKT, M, every, IC, h, L, NLfunc, fP, k, name)
	runMMT(IFRK3, M, every, IC, h, L, NLfunc, fP, k, "IFRK3-"*string(Int(h*1000000),pad=6)*"-d"*string(deg))
else 
	runMMT(scheme, M, every, IC, h, L, NLfunc, fP, k, name)
end

using PyPlot, LaTeXStrings
saveEnergy!(k, N, T, name, scheme=scheme, h=h, ES=false, deg=deg)
saveEnergy!(k, N, T, name, scheme=scheme, h=h, ES=true)

saveEnergy!(k, N, T, "IFRK3-"*string(Int(h*1000000),pad=6)*"-d"*string(deg), scheme=scheme, h=h, ES=false)
saveEnergy!(k, N, T, "IFRK3-"*string(Int(h*1000000),pad=6)*"-d"*string(deg), scheme=scheme, h=h, ES=true)

