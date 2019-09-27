using LinearAlgebra

function expω(L::Int)
	c = zeros(Float64, L)
	for i = 1 : L
		c[i] = 1 / factorial(i)
	end
	return c
end

function Pade(c::Array{T,1}, n::Int, m::Int) where T<:AbstractFloat
	if m >= n
		RHS = -c[n+1:end]
		A = c[n:n+m-1]
		for i = 2 : n
			A = hcat(A, c[n-(i-1): n-i+m])
		end
		for i = 1:m-n
			A = hcat(A, vcat(zeros(i), c[1:m-i]))
		end 
		#return A, RHS
		
		b = A\RHS
		a = b[1:n]+c[1:n]
		for i = 1 : n
			a[i] += dot(b[1:i-1], c[i-1:-1:1])
		end
		return a, b 
	else #m<n
		RHS = -c[n+1:end]
		A = c[n:n+m-1]
		for i = 2 : m
			A = hcat(A, c[n-(i-1): n+m-i])
		end
		b = A\RHS
		c = vcat([1], c) #c_0=c[1]=1, c_1=c[2] and so on..
		a = c[2:n+1]       #a[1] = a_1 = c_1 (since we set a_1 = 1)
		for i = 1 : m
			a[i] += dot(b[1:i], c[i:-1:1])
		end
		return a, b
	end
end

#a = [1/2; 1/9; 1/72; 1/1008; 1/30240]
#b = [-1/2; 1/9; -1/72; 1/1008; -1/30240]
using PyPlot, LaTeXStrings

x = range(0, pi/2, length=1001)

function plotPade!(c, n, m, x)
	a, b = Pade(c, n, m)
	xx = im*x
	for i = 2 : length(c)
		xx = hcat(xx, (im*x).^i)
	end
	p = (1 .+ xx[:,1:length(a)]*a) ./(1 .+ xx[:,1:length(b)]*b)
	EEE = exp.(im*x)
	semilogy(x, abs.(angle.(p)-angle.(EEE)), label="Pade "*string(n)*"-"*string(m))
end
#=sum to 8
n = 2; m = 6
c = expω(n+m)
plotPade!(c, n, m, x)
#####################
n = 3; m = 5
c = expω(n+m)
plotPade!(c, n, m, x)
#####################
n = 4; m = 4
c = expω(n+m)
plotPade!(c, n, m, x)
#####################

n = 5; m = 3
c = expω(n+m)
plotPade!(c, n, m, x)
#####################
n = 6; m = 2
c = expω(n+m)
plotPade!(c, n, m, x)
#####################
=#

#=
n = 2; m = 2
c = expω(n+m)
plotPade!(c, n, m, x)
#####################
n = 3; m = 3
c = expω(n+m)
plotPade!(c, n, m, x)
#####################
n = 4; m = 4
c = expω(n+m)
plotPade!(c, n, m, x)
#####################

n = 5; m = 5
c = expω(n+m)
plotPade!(c, n, m, x)
#####################
n = 6; m = 6
c = expω(n+m)
plotPade!(c, n, m, x)
#####################
n = 7; m = 7
c = expω(n+m)
plotPade!(c, n, m, x)
#####################
n = 8; m = 8
c = expω(n+m)
plotPade!(c, n, m, x)
#####################

xt = [0, .125*pi, .25*pi, .375*pi, .5*pi]#, 3*pi, 4*pi]
axvline.(xt[2:end], color="k")
xlabels = ("0", L"\frac{\pi}{8}",  L"\frac{\pi}{4}", L"\frac{3\pi}{8}", L"\frac{\pi}{2}")#, L"3\pi", L"4\pi")
xticks(xt, xlabels)
legend()
title("Pade approx errors of exp(ix)")
=#