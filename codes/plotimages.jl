include("readwrite.jl")

using PyPlot, LaTeXStrings
name = "bigICRK3"
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
end