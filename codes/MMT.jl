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

@inline function NLfunc(zhat::Array{ComplexF64,1}, fp::funcparams, k)
	#forcing term
	N = length(zhat)
	F = zeros(length(k)); 
	F[[7,8,9,-7+N,-8+N,-9+N]] .= 0.2
	#damping term
	D = -196.61 * (abs.(k).^(-8)) - 5.39* (abs.(0.001 * k) .^ 16); D[1] = 1

	if fp.β != 0
		zr = ifft(abs.(k) .^(β/4) .* zhat)
    	return -im* fp.λ * abs.(k) .^(β/4) .* fft(abs.(zr).^2 .* zr)
    else
    	zr = ifft(zhat)
    	return -im* fp.λ*fft(abs.(zr).^2 .*zr) + (F+D) .* zhat
    end
end
#add dissipitation & forcing in NLfunc

struct eRKTableau
	A
	b
	eRKTableau(A,b) = new(copy(A),copy(b))
end

function IFRK_step(z::Array{ComplexF64,1}, h::Float64, 
	L, NLfunc::Function, nlfP::funcparams, 
	RKT::eRKTableau, ks)
	#k = zeros(eltype(z), length(RKT.b), length(z))
	ks[1,:] = NLfunc(z, nlfP, k)
	for i = 2 :length(RKT.b)
		PP=h*Transpose(ks[1:i-1,:])*RKT.A[i,1:i-1]
		ks[i,:] = NLfunc(z + PP, nlfP, k)
	end
	vb = Transpose(ks)*b
	#return k, vb, L * (z+ (h*vb))
	return L .* (z+ (h*vb))
end

function IFRK!(M::Int, every::Int, IC::Array{ComplexF64,1}, h::Float64, 
	L, NLfunc::Function, fP::funcparams, RKT::eRKTableau; name::String)
	# FFT into Fourier space
	#sqrtN = sqrt(length(IC))
	zhat = fft(IC)#/sqrtN
	newtxt!(zhat, name=name)
	ks = zeros(eltype(zhat), length(RKT.b), length(IC)) #RK stages (allocate memory)
	k = range(0, length(IC)-1, length=length(IC)) #wave numbers
	L = exp.(h*L) #IF(L)
	for t = 1 : M
		zhat = IFRK_step(zhat, h, L, NLfunc, fP, RKT, ks)
		if rem(t,every)==1
			addtxt!(zhat, name=name)
		end

	end
end

#Parameters
λ = 1   #Defocusing MMT model
α = 1/2
β = 0

N = 2^11
x = range(0,2*pi, length=N);
IC = ifft(randn(ComplexF64, N)); #Array{ComplexF64,1}(sin.(x))

fP = funcparams(α, β, λ);
k = vcat(collect(0:N/2), collect(-N/2+1:-1))
L = -im*abs.(k).^α

#A = hcat([0; .5])
#b = [0;1]
A = hcat([0; .5; 0; 0], [0; 0; .5; 0], [0; 0; 0; 1.0])
b = [1/6; 1/3; 1/3; 1/6]
RKT = eRKTableau(A, b)
h = 0.005;
T = 10000
M = Int(T/h);
every = Int(M/1000)

name = "FDadded"
#IFRK!(M, every, IC, h, L, NLfunc, fP, RKT, name=name)
solhat = readCfile(name) / sqrt(N)
sol = ifft(solhat, 2)

E = transpose(sum(abs.(solhat).^2, dims = 1)/size(solhat)[1])


using PyPlot, LaTeXStrings
kind = vcat(collect(Int(N/2)+1:N-1), collect(1:Int(N/2)))
loglog(k[1:Int(end/2)], E[1:Int(end/2)], label="computed")
loglog(k[1:Int(end/2)], 0.24 ./k[1:Int(end/2)], label="weak turbulence spectrum")
xlabel("Wave Number")
ylabel("n(k)")
legend()
title("Energy spectrum?")

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


