using LinearAlgebra

include("RT.jl")

struct eRKTableau
	A
	b
	eRKTableau(A,b) = new(copy(A),copy(b))
end

function IFRK(z::Array{ComplexF64,1}, h::Float64, L, NLfunc::Function, nlfP::NLfuncparams, RKT::eRKTableau)
	k = zeros(eltype(z), length(RKT.b), length(z))
	k[1,:] = NLfunc(z, nlfP)
	for i = 2 :length(RKT.b)
		PP=h*Transpose(k[1:i-1,:])*RKT.A[i,1:i-1]
		k[i,:] = NLfunc(z + PP, nlfP)
	end
	vb = Transpose(k)*b
	return k, vb, L * (z+ (h*vb))
end

#=function IFRK(z::Array{ComplexF64,1}, h::Float64, L, NLfunc::Function, nlfP::NLfuncparams, RKT::eRKTableau)
	k = Array{ComplexF64,2}(undef, length(RKT.b), length(z))
	k[1,:] = NLfunc(z, nlfP)
	vn = b[1]*k[1,:]
	for i = 2 :length(RKT.b)
		PP= zeros(ComplexF64,3)
		for j = 1 : i-1
			PP += RKT.A[i,j]*k[j,:]
		end
		k[i,:] = NLfunc(z + PP, nlfP)
		vn += b[i]*k[i,:]
	end
	return k, vb, L * (z+ (h*vb))
end=#

h = 0.3
ω = [-1,3,-2]
L = Diagonal(exp.(im*h*ω))

ϵ = 0.01;
C = floatRT(5)
nlf = NLfuncparams(ϵ, C)

A = hcat([0; .5; 0; 0], [0; 0; .5; 0], [0; 0; 0; 1.0])
b = [1/6; 1/3; 1/3; 1/6]
RKT = eRKTableau(A, b)

IC = onUnitCircle(3);
kk, vv, gg = @time(IFRK(IC, h, L, NLfunc, nlf, RKT))
kk2, vv2, gg2= @time(IFRK4(h, IC; ω=ω, ϵ=ϵ, C=C))
