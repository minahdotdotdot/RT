include("MMT.jl")
include("al_approx.jl")

function cfromx(x::Vector{T}) where T<:AbstractFloat
	pushfirst!(x, T(0))
	s = length(x)
	c = zeros(T, s+1, s)
	for i = 1 : s
		for j = 1 : i-1
			c[i,j] = x[i]-x[j]
		end
		c[i,i] = x[i]
	end
	c[end,1:s] = 1 .- x
	return c
end


function fillc(x, h, L)
	uniquec = unique(x)
	uniquet = Dict{Float64, typeof(L)}();
	for u in uniquec
		push!(uniquet, u => exp.((h*u)*L))#+1e-16*randn(length(L)))
	end
	tmp = Matrix{typeof(L)}(undef,size(x))
	for i = 1 : size(x)[1]
		for j = 1 : size(x)[2]
			tmp[i,j] = uniquet[x[i,j]]
		end
	end
	return tmp
end

function fillc_rat(c, h, L, deg::Int)
	name = "../data/MMT_deg"*string(deg)
	tmp = Array{typeof(L),1}(undef,0)
	for i = 1 : length(c)
		push!(tmp, evalrat(-h*c[i]*L, name))
		push!(tmp, evalrat(h*c[i]*L, name))
	end
	push!(tmp, evalrat(h*L, name))
	return tmp
end

#RK2

#RK3
A = hcat([0; .5; -1],[0; 0; 2]);
b = [1/6; 2/3; 1/6];
x = [1/2; 1]; 
x = cfromx(x)
c =fillc(x,h,L);
#crat = fillc_rat(x, h, L, 6);
IFRK3 = eRKTableau(A, b, c, x);
#IFRK3_rat = eRKTableau(A, b, crat)

#RK4
A = hcat([0; .5; 0; 0], [0; 0; .5; 0], [0; 0; 0; 1.0]); #3by4
b = [1/6; 1/3; 1/3; 1/6];
x = [1/2; 1/2; 1]; 
x = cfromx(x)
c =fillc(x,h,L);
#crat = fillc_rat(x, h, L, 6);
IFRK4 = eRKTableau(A, b, c, x); 
#IFRK4_rat = eRKTableau(A, b, crat);  =#

IFdict = Dict( 
	"IFRK3" => IFRK3, 
	"IFRK4" => IFRK4
	);
