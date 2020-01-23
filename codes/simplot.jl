include("setupMMT.jl")
runMMT(scheme, M, every, IC, h, L, NLfunc, fP, k, name)

using PyPlot, LaTeXStrings
plotEnergy!(k, N, T, name)