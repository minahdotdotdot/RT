using Remez, PyPlot, LaTeXStrings

function horner(x::Array{T,1}, coef) where T<:AbstractFloat
	y = coef[end]
	for i = length(coef)-1 : -1 : 1
		y = (y.* x) .+ coef[i]
	end
	return y
end
numdeg = 12
for n = numdeg-2 : -1 : 2
	m = numdeg-n;
#for numdeg = 24 : -2 : 4
#numdeg = 16; 
dendeg = 0
Ncos, Dcos, Ecos, Xcos = ratfn_minimax(cos, (0, 2*pi), numdeg, dendeg);
Nsin, Dsin, Esin, Xsin = ratfn_minimax(sin, (0, 2*pi), numdeg, dendeg);

x = collect(range(0, stop=2*pi, length=1001))
rcosap = Array{Float64,1}(horner(x, Ncos))
rsinap = Array{Float64,1}(horner(x, Nsin))

using PyPlot


include("pade.jl")
#m = Int(numdeg/2); n = m;
x = collect(range(0, stop=2*pi, length=1001))
a_cos, b_cos = Pade(Ncos, n, m)
a_sin, b_sin = Pade(Nsin, n, m)

a_eiz, b_eiz = Pade(Ncos + im*Nsin, n, m)
peizap = Array{Complex{Float64}, 1}(horner(x, a_eiz) ./ horner(x, b_eiz))

pcosap = Array{Float64,1}(horner(x, a_cos) ./ horner(x, b_cos))
psinap = Array{Float64,1}(horner(x, a_sin) ./ horner(x, b_sin))


fig, ax = subplots()
subplot(131)
#plot(x, pcosap, label="pade of remez cos")
#plot(x, rcosap, label="remez cos")
#plot(x, cos.(x), label="actual cos")
semilogy(x, abs.(rcosap - cos.(x)), label="Remez cos err")
semilogy(x, abs.(pcosap - cos.(x)), label="pade-Remez cos err")
ylim(1e-20, 1)
legend()

subplot(132)
#plot(x, psinap, label="pade of remez sin")
#plot(x, rsinap, label="remez sin")
#plot(x, sin.(x), label="actual sin")
semilogy(x, abs.(rsinap - sin.(x)), label="Remez sin err")
semilogy(x, abs.(psinap - sin.(x)), label="pade-Remez sin err")
ylim(1e-20, 1)
legend()

subplot(133)
#plot(x, psinap, label="pade of remez sin")
#plot(x, rsinap, label="remez sin")
#plot(x, sin.(x), label="actual sin")
semilogy(x, abs.(real.(peizap) - cos.(x)), label="pade-Remez eiz real err")
semilogy(x, abs.(imag.(peizap) - sin.(x)), label="pade-Remez eiz imag err")
ylim(1e-20, 1)
legend()

suptitle(string(n)*"/"*string(m)*"Pade approximations")
savefig(string(numdeg+1)*"uneven rational_n="*string(n)*".png")
close(fig)
end
#semilogy(x, abs.(pcosap - rcosap), label="pade remez cos err")
#semilogy(x, abs.(psinap - rcosap), label="pade remez sin err")




#generate
#=
dendeg = 0;
M = 10
Ep = zeros(BigFloat, M-1,2)
coeffs = Array{Any,2}(undef, M-1,2)
for i = 2:M 
	print(i, "\n")
	Ncos, Dcos, Ecos, Xcos = ratfn_minimax(cos, (0, 2*pi), 2*i, dendeg);
	coeffs[i-1,1]=Ncos
	Ep[i-1,1] = Ecos
	Nsin, Dsin, Esin, Xsin = ratfn_minimax(sin, (0, 2*pi), 2*i, dendeg);
	Ep[i-1,2] =Esin
	coeffs[i-1,2]=Nsin
end

coeff20=coeffs[end,:]
cos20 = Array{Float64,1}(coeff20[1])
sin20 = Array{Float64,1}(coeff20[2])
eiz20 = cos20 + im*sin20 =#


#=

include("pade.jl")
x = collect(range(0, stop=2*pi, length=1001))
a_cos, b_cos = Pade(cos20, 10, 10)
a_sin, b_sin = Pade(sin20, 10, 10)


approx = horner(x, eiz20)

cosapprox = horner(x, a_cos) ./ horner(x, b_cos)
sinapprox = horner(x, a_sin) ./ horner(x, b_sin)

semilogy(x, abs.(cos.(x) - cosapprox), label="cos error")
semilogy(x, abs.(sin.(x) - sinapprox), label="sin error")
padeapprox = cosapprox + (im*sinapprox)

a, b = Pade(eiz20, 10, 10)
padeapprox2 = horner(x, a) ./ horner(x, b)

#semilogy(x, abs.(padeapprox - approx), label="pade from remez")
#semilogy(x, abs.(exp.(im*x) - approx), label="remez err")
#semilogy(x, abs.(exp.(im*x) - padeapprox), label="pade err")
xlabel("Domain")
ylabel("Error")
title("Approximation errors of "*L"exp(iz)")
legend()
=#


#=
Erat=hcat(readdlm("../txtfiles/"*"cosrat.txt")[1,:], readdlm("../txtfiles/"*"sinrat.txt")[1,:])
Epol=hcat(readdlm("../txtfiles/"*"cospol.txt")[1,:], readdlm("../txtfiles/"*"sinpol.txt")[1,:])

fig, ax = subplots()
ax.set_yscale("log")

xt = 2*collect(2:M)
xlabels = string.(xt)
xticks(xt, xlabels)
scatter(xt[E[:,1].!=0], E[:,1][E[:,1].!=0], label="Cos rat")
scatter(xt[E[:,2].!=0], E[:,2][E[:,2].!=0], label="Sin rat")
semilogy(xt[Ep[:,1].!=0], Ep[:,1][Ep[:,1].!=0], label="Cos pol")
semilogy(xt[Ep[:,2].!=0], Ep[:,2][Ep[:,2].!=0], label="Sin pol")
legend()

xticks(xt, xlabels)
legend()
ylabel("Error")
xlabel("Total Degree")
title("Error of minimax approximations")
=#

#=Ncos= Array{Float64, 1}(Ncos); 
x = collect(range(0, stop=2*pi, length=1001))
approxCos = horner(x, Ncos) #poly_eval.(Ncos, x) #
errCos = Array{Float64, 1}(cos.(x)-approxCos)

approxSin = horner(x, Nsin) #poly_eval.(Ncos, x) #
errSin = Array{Float64, 1}(sin.(x)-approxSin)

plot(x, errCos, label="cos"*string(numdegc))
plot(x, errSin, label="sin"*string(numdegs))=#


#=fig, ax=subplots(1,2)
#ax.set_yscale("log")
subplot(121)
semilogy(collect(2:M)[E[:,1].!=0], E[:,1][E[:,1].!=0], label="Cos rat")
semilogy(collect(2:M)[E[:,2].!=0], E[:,2][E[:,2].!=0], label="Sin rat")
legend()

ylim([eps(Float64)/10, 1e-01])
xt = collect(2:M)
xlabels = string.(xt)
xticks(xt, xlabels)

ylabel("Error")
xlabel("Equal Degree Rationals")
#title("Error of minimax approximations")

subplot(122)
semilogy(2*collect(2:M)[Ep[:,1].!=0], Ep[:,1][Ep[:,1].!=0], label="Cos pol")
semilogy(2*collect(2:M)[Ep[:,2].!=0], Ep[:,2][Ep[:,2].!=0], label="Sin pol")
ylim([eps(Float64)/10, 1e-01])
xt = 2*collect(2:M)
xlabels = string.(xt)
xticks(xt, xlabels)
legend()
ylabel("Error")
xlabel("Polynomials")

suptitle("Error of minimax approximations")
=#

#title("Error of Degree "*string(numdeg)*" minimax approximations")