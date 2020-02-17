include("setupMMT.jl")
using MAT
if scheme ∈ ["IFRK3", "IFRK3_rat", "IFRK4", "IFRK4_rat"]
	

file = matopen("../data/L.mat", "w")
write(file, "L", L)
close(file)
