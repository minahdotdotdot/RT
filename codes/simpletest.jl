include("MMT.jl")

IC=[im*pi]
L = [0.0];k=0;fP=0
function NLfunc(zhat::Array{ComplexF64,1}, fP::funcparams, k)
    return 3*im*zhat#
end

function NLfunc(zhat::Array{ComplexF64,2}, fP::funcparams, k)
    return 3*im*zhat
end
λ = 1;  #Defocusing MMT model
α = 1/2;
β = 0;
F = 0.05;
D = [2.51e-57, 16];
fP = funcparams(α, β, λ, F, D);
# time-step, ND final time, save "every"
h = 1e-3;
T = 2
M = Int(T/h);
#T = floor(Int,M*h)
every = Int(M/1000) # save solution at only 101 time locations.

include("IF_methods.jl")
name = "IFRK3"
RKT=RK3;
IFRK!(M, every, IC, h, L, NLfunc, fP, RKT, k, name=name);
print(name," done.\n")

name="IFRK4"
RKT=RK4
IFRK!(M, every, IC, h, L, NLfunc, fP, RKT, k, name=name) 
print(name," done.\n")

include("ETD_methods.jl")
name="ETDRK3"
RKT=ETDRK3
ETDRK!(M, every, IC, h, L, NLfunc, fP, RKT, k, name=name)
print(name," done.\n")

name="ETDRK4"
RKT=ETDRK4
ETDRK!(M, every, IC, h, L, NLfunc, fP, RKT, k, name=name)
print(name," done.\n")

name="ETDRK4B"
RKT=ETDRK4B
ETDRK!(M, every, IC, h, L, NLfunc, fP, RKT, k, name=name)
print(name," done.\n")

Names = ["IFRK3", "IFRK4", "ETDRK3", "ETDRK4", "ETDRK4B"]
solutions = Array{Any,1}(undef, length(Names))

truesol = IC[1]/1000*exp.(3*im*range(0, stop=T, length=Int(M/every+1)))
	#hcat(IC[1]/1000*exp.(L[1]*range(0, stop=T, length=Int(M/every+1))),
	#IC[2]/1000*exp.(L[2]*range(0, stop=T, length=Int(M/every+1)))
	#)
for i = 5:5
	name=Names[i];
	solutions[i]=readCfile(name);
	#plot(-im*solutions[i], label=name)
	plot(abs.(truesol[2:end]-solutions[i][2:end]), label=name)
	#plot(abs.(truesol[2:end,:]-solutions[i][2:end,:])[:,1], label=name*"1")
	#plot(abs.(truesol[2:end,:]-solutions[i][2:end,:])[:,2], label=name*"2")
end
#plot(-im*truesol, label="truesol")
title("Testing NonLinear Operator: Absolute Errors")
legend()