include("readwrite.jl")
using FFTW, LinearAlgebra

struct funcparams
	α:: Float64
    β:: Float64
    λ:: Complex{Float64}
    funcparams(α, β, λ) = new(copy(α), copy(β), copy(λ))
end

@inline function Lfunc(zhat::Array{ComplexF64,1}, fp::funcparams)
	return abs.(k).^fp.α .*zhat
end

@inline function NLfunc(zhat::Array{ComplexF64,1}, fp::funcparams, k, FD)
	if fp.β != 0
		zr = ifft(abs.(k) .^(β/4) .* zhat)
    	return -im* fp.λ * abs.(k) .^(β/4) .* fft(abs.(zr).^2 .* zr) + FD .* zhat
    else
    	zr = ifft(zhat)
    	return -im* fp.λ*fft(abs.(zr).^2 .*zr) + FD .* zhat
    end
end

struct eRKTableau
	A # matrix
	b # stage weights
	c # t increments
	eRKTableau(A, b, c) = new(copy(A),copy(b), copy(c))
end

function IFRK_step(z::Array{ComplexF64,1}, h::Float64, 
	L, NLfunc::Function, nlfP::funcparams, 
	RKT::eRKTableau, ks, k, FD)
	#k = zeros(eltype(z), length(RKT.b), length(z))
	ks[1,:] = NLfunc(z, nlfP, k, FD)
	for i = 2 :length(RKT.b)
		PP=h*Transpose(ks[1:i-1,:])*RKT.A[i,1:i-1]
		ks[i,:] = exp.(-h*RKT.c[i]*L) .* NLfunc(exp.(h*RKT.c[i]*L) .*(z + PP), nlfP, k, FD)
	end
	vb = Transpose(ks)*b
	#return k, vb, L * (z+ (h*vb))
	return exp.(h*L) .* (z+ (h*vb))
end

function IFRK!(M::Int, every::Int, IC::Array{ComplexF64,1}, h::Float64, 
	L, NLfunc::Function, fP::funcparams, RKT::eRKTableau, k; name::String)
	# FFT into Fourier space
	#sqrtN = sqrt(length(IC))
	zhat = fft(IC)#/sqrtN
	newtxt!(zhat, name=name)
	ks = zeros(eltype(zhat), length(RKT.b), length(IC)) #RK stages (allocate memory)
	#L = exp.(h*L) #IF(L)
	#forcing term
	N = length(zhat)
	FD = zeros(length(k)); 
	FD[[7,8,9,-7+N,-8+N,-9+N]] .= 0.2
	#add damping term
	FD[2:end] += -196.61 * (abs.(k[2:end]).^(-8)) - 5.39* (abs.(0.001 * k[2:end]) .^ 16); FD[1] += 1
	for t = 1 : M
		zhat = IFRK_step(zhat, h, L, NLfunc, fP, RKT, ks, k, FD)
		if rem(t,every)==1
			addtxt!(zhat, name=name)
		end

	end
end

# Problem Parameters
λ = 1   #Defocusing MMT model
α = 1/2
β = 0
fP = funcparams(α, β, λ);

# Numerical Simulation Parameters
N = 2^11
x = range(0,stop=2*pi-1/N, length=Int(N));
IC = ifft(randn(ComplexF64, N)); #Array{ComplexF64,1}(sin.(x))
k = vcat(collect(0:N/2-1), collect(-N/2:-1))
L = -im*abs.(k).^α

# Set-up IFRK
A = hcat([0; .5; 0; 0], [0; 0; .5; 0], [0; 0; 0; 1.0])
b = [1/6; 1/3; 1/3; 1/6]
c = [0; 1/2; 1/2; 1]
RKT = eRKTableau(A, b, c)

# time-step, ND final time, save "every"
h = 0.005;
T = 10000
M = Int(T/h);
every = Int(M/1000) # save solution at only 1001 time locations.

name = "FDaddedRKfixed"
#IFRK!(M, every, IC, h, L, NLfunc, fP, RKT, k, name=name) 
#solhat = readCfile(name)[11:end,:]
#sol = ifft(solhat, 2)

#E = transpose(sum(abs.(solhat).^2, dims = 1)/size(solhat)[1])


using PyPlot, LaTeXStrings
kind = vcat(collect(Int(N/2)+1:N-1), collect(1:Int(N/2)))
kindnz = vcat(collect(Int(N/2)+1:N-1), collect(2:Int(N/2)))
fig, ax = subplots()
#semilogy(k[kindnz], E[kindnz], label="computed")
#semilogy(k[2:Int(end/2)], 1e-24 *(k[2:Int(end/2)]).^(-1/3), label=L"Ck^{-1/3}")
#semilogy(k[2:Int(end/2)], 1e-24 *(k[2:Int(end/2)]).^(-1/2), label=L"Ck^{-1/2}")
loglog(k[2:Int(end/2)], E[2:Int(end/2)], label="computed")
loglog(k[2:Int(end/2)], 1e-24 *(k[2:Int(end/2)]).^(-1/3), label=L"Ck^{-1/3}")
loglog(k[2:Int(end/2)], 1e-24 *(k[2:Int(end/2)]).^(-1/2), label=L"Ck^{-1/2}")
xlabel("Wave Number")
ylabel("n(k)")
legend()
title(L"t\in[110, 10000]")
#savefig(name*"ES.png")
#close(fig)

#=plot(x, abs.(sol[1,:]).^2, label="IC")
plot(x, abs.(sol[60000,:]).^2, label="middle")
plot(x, abs.(sol[end,:]).^2, label="end")
legend()
xlabel("x")
ylabel(L"|\psi(x)|^2")
title("Amplitude Squared")=#

#=
for i = 1 : 50: size(sol)[1]
	fig, ax = subplots(); 
	plot(x, real.(sol[i,:]), label="real")
	plot(x, imag.(sol[i,:]), label="imag")
	legend()
	savefig("../movie/"*name*"_"*string(i))
	close(fig)
end

colorbar(imshow(real.(sol)))
=#


