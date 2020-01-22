using MAT, PyPlot, LaTeXStrings, LinearAlgebra
#=Poles = Array{Any,1}(undef, 5)
Zeros = Array{Any,1}(undef, 5)
for i = 4 : 8
	file = matread("polzer"*string(i)*".mat");
	Poles[i-3] = file["pol"]
	Zeros[i-3] = file["zer"]
end=#

function evalbary(z::Vector{Complex{T}}, name::String) where T<:AbstractFloat
	info = matread(name*".mat")
	#fj = exp.(info["zj"])
	fj = info["fj"];
	CC = 1 ./ broadcast(-, z, transpose(info["zj"]))
	r = (CC*(info["wj"].*fj))./(CC*info["wj"])
	#=for i = 1 : length(r)
		if isnan(r[i]) == true
			if i == 1 || i == length(r)
				r[i] = (r[2]+r[end-1])/2
			elseif i > 1 && i < length(r)
				r[i] = (r[i+1]+r[i-1])/ 2
			else
				r[i] = r[i-1]
			end
		end
	end=#
	return r
end

function evalrat(z::Vector{Complex{T}}, name::String) where T<:AbstractFloat
	info = matread(name*".mat")
	D = info["pol"]
	N = info["zer"]
	denom = (z .-D[1])
	for i = 2 : length(D)
		denom .*= (z .-D[i])
	end
	num = (z .-N[1])
	for i = 2 : length(N)
		num .*= (z .-N[i])
	end
	r = num ./denom
	nr, n = evalnormal(name)
	return r*nr + (1-n)
end

function evalnormal(name::String)
	info = matread(name*".mat")
	D = info["pol"]
	N = info["zer"]
	n = evalbary(zeros(ComplexF64,2),name)[1]
	r = prod(-N)/prod(-D)
	return n/r, n
end

function evalrat(L::Matrix{Complex{T}}, name::String) where T<:AbstractFloat
	info = matread(name*".mat")
	D = info["pol"]
	N = info["zer"]
	denom = (L -D[1]*I)
	for i = 2 : length(D)
		denom *= (L -D[i]*I)
	end
	num = (L -N[1]*I)
	for i = 2 : length(N)
		num *= (L -N[i]*I)
	end
	#expL = num*inv(denom)
	expL = (D'\N')'
	nr, n = evalnormal(name)
	return expL*nr +(1-n)*I 
end

function evalrat(L::Diagonal{Complex{T}}, name::String) where T<:AbstractFloat
	return Diagonal(evalrat(diag(L), name))
end

#=
function plotcomplex!(errexpL::Array{Complex{T},1}, 
	i::Int, name::String="Matrix ") where T <: AbstractFloat
	#approxexpL = evalrat(h*L, "../codes/polzer"*string(i))
	#errexpL = trueexpL - approxexpL
	J = Int(size(errexpL)[1]/2)+1
	#norm(errexpL,2)

	fig, ax = subplots()
	for j = 2 : J
		crgb = get(rainbow, j/J)
		plot(real.(errexpL[j-1:j]), imag.(errexpL[j-1:j]), 
			c=(crgb.r, crgb.g, crgb.b))
	end
	axvline(0, color="k")
	axhline(0, color="k")
	title("Degree"*string(i)* name*"Approximation")
end

function plotcomplex!(errexpL::Diagonal{Complex{T},Array{Complex{T},1}}, 
	i::Int, name::String="Matrix ") where T <: AbstractFloat
	#approxexpL = evalrat(h*L, "../codes/polzer"*string(i))
	#errexpL = trueexpL - approxexpL
	J = Int(size(errexpL)[1]/2)+1
	#norm(errexpL,2)

	fig, ax = subplots()
	for j = 2 : J
		crgb = get(rainbow, j/J)
		plot(real.(errexpL.diag[j-1:j]), imag.(errexpL.diag[j-1:j]), 
			c=(crgb.r, crgb.g, crgb.b))
	end
	axvline(0, color="k")
	axhline(0, color="k")
	title("Degree"*string(i)* name*"Approximation")
end

function plotseparate!(errexpL, i::Int, J::Int)
	fig, ax = subplots()
	#approxexpL = evalrat(h*L, "../codes/polzer"*string(i))
	#errexpL = trueexpL - approxexpL
	plot(collect(1:J), real.(errexpL.diag[1:J]), label=string(i)*"real")
	plot(collect(1:J), imag.(errexpL.diag[1:J]), label=string(i)*"imag")
	xlabel("Wave numbers")
	legend()
	title("Degree"*string(i)* "Matrix Approximation")
end


function semilogyseparate!(errexpL, i::Int, J::Int)
	#approxexpL = evalrat(h*L, "../codes/polzer"*string(i))
	#errexpL = trueexpL - approxexpL
	semilogy(collect(1:J), abs.(real.(errexpL.diag[1:J])), label=string(i)*"real")
	semilogy(collect(1:J), abs.(imag.(errexpL.diag[1:J])), label=string(i)*"imag")
	xlabel("Wave numbers")
	legend()
end


z = collect(im*range(0, stop = 2*pi, length = 2001));
#truesol = exp.(z);
#info6, r6 = evalrat(z, "polzer6");
#Try on matrix
include("MMT.jl")
# Numerical Simulation Parameters
N = 2^11
x = range(0,stop=2*pi-1/N, length=Int(N));
IC = randn(ComplexF64, N)*sqrt(N)/1000; IC[1]=0.0; IC = ifft(IC); 
#IC = cos.(x) + im*sin.(x)
k = vcat(collect(0:N/2), collect(-N/2+1:-1)) # implies domain length is 2π
kind = vcat(collect(Int(N/2)+2:N), collect(1:Int(N/2)+1))
kindnz = vcat(collect(Int(N/2)+2:N), collect(2:Int(N/2))) # indexing w/o zero mode




# Linear operator (depends on k)
h = 0.005;
L = -im*abs.(k).^fP.α; L = Diagonal(L);
trueexpL = exp(h*L);


r = Array{Any, 1}(undef, 5);
errexpL = Array{Any, 1}(undef, 5)
for i = 4 : 8
    r[i] = evalrat(h*L, "polzer"*string(i))
    errexpL[i] = trueexpL - r[i]
end

import ColorSchemes.rainbow
for i = 4 : 8
	approxexpL = evalrat(h*L, "polzer"*string(i))
	errexpL = trueexpL - approxexpL
	J = Int(N/2+1)
	#norm(errexpL,2)

	fig, ax = subplots()
	for j = 2 : J
		crgb = get(rainbow, j/J)
		plot(real.(errexpL.diag[j-1:j]), imag.(errexpL.diag[j-1:j]), 
			c=(crgb.r, crgb.g, crgb.b))
	end
	axvline(0, color="k")
	axhline(0, color="k")
	colorbar(cmap=rainbow)
	#semilogy(k[1:Int(N/2)], abs.(real.(errexpL.diag[1:Int(N/2)])), label="real errors")
	#semilogy(k[1:Int(N/2)], abs.(imag.(errexpL.diag[1:Int(N/2)])), label="imag errors")
	#plot(k[kind], real.(errexpL.diag[kind]), label="real errors")
	#plot(k[kind], imag.(errexpL.diag[kind]), label="imag errors")
	#legend()
	title("Degree"*string(i)* "Matrix Approximation")
	savefig("complexM"*string(i)*".png")
close(fig)
end

#=
for i = 4 : 8
	name = "polzer"*string(i)
	approx = evalbary(z, name)
	approx2 = evalrat(z, name)
	E = truesol - approx
	E2 = truesol - approx2
	discr = approx - approx2
	fig, ax = subplots()
	plot(real.(E), imag.(E), label="bary error")
	plot(real.(E2), imag.(E2), label="rat error")
	plot(real.(discr), imag.(discr), label="bary-rat discrepancy")
	#ermax = max(maximum(real.(E)), maximum(imag.(E)))
	#print(ermax, "\n")
	#plot(ermax*cos.(z), ermax*sin.(z), label="errbound")
	#plot(real.(truesol), imag.(truesol), label="true")
	#plot(real.(approx), imag.(approx), label="approx")
	legend()
	xlabel("Real")
	ylabel("Imaginary")
	title("Degree "*string(i)*" Rational Approximation")
	savefig(name*".png")
	close(fig)
end
=#
#=
for i = 4 : 7
	approx = evalrat(z, Poles[i-3][:,1], Zeros[i-3][:,1])
	E = truesol - approx
	print(string(i)*": ", maximum(abs.(E)), "\n")
	semilogy(-im*z, abs.(E), label="deg"*string(i))
end
legend()
xlabel("Domain = "*L"im [0,2\pi]")
title("AAA-Lawson rational approximations")
=#
=#
