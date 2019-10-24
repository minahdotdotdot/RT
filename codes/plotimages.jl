include("readwrite.jl")

using PyPlot, LaTeXStrings
#=name = "bigICRK3"
solhat = readCfile(name)
sol = ifft(solhat, 2)

T = size(sol)[1]
for i = 1 : 5: T
	fig, ax = subplots()
	plot(real.(sol[i,:]), label="real")
	plot(imag.(sol[i,:]), label="imag")
	ylim([min(minimum(real.(sol)), minimum(imag.(sol))), 
		max(maximum(real.(sol)), maximum(imag.(sol)))])
	legend()
	savefig("../plots/"*name*string(i)*".png")
	close(fig)
end=#

f200 = readdlm("../txtfiles/f200.txt"); plot(f200, label="0.2");
f100 = readdlm("../txtfiles/f100.txt"); plot(f100, label="0.1");
f050 = readdlm("../txtfiles/f050.txt"); plot(f050, label="0.05");
f025 = readdlm("../txtfiles/f025.txt"); plot(f025, label="0.025");
f010 = readdlm("../txtfiles/f010.txt"); plot(f010, label="0.010")
legend()
title("Varying Forcing")