include("MMT.jl")

L = [3*im];k=0;fP=0
function NLfunc(zhat::Array{ComplexF64,1}, fP::funcparams, k)
    return 0*zhat#3*im*zhat#
end

function NLfunc(zhat::Array{ComplexF64,2}, fP::funcparams, k)
    return 0*zhat#3*im*zhat
end

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


truesol = IC[1]/1000*exp.(L[1]*range(0, stop=T, length=Int(M/every+1)))
	#hcat(IC[1]/1000*exp.(L[1]*range(0, stop=T, length=Int(M/every+1))),
	#IC[2]/1000*exp.(L[2]*range(0, stop=T, length=Int(M/every+1)))
	#)
for i = 1:5
	name=Names[i];
	solutions[i]=readCfile(name);
	#plot(-im*solutions[i], label=name)
	plot(abs.(truesol[2:end]-solutions[i][2:end]), label=name)
	#plot(abs.(truesol[2:end,:]-solutions[i][2:end,:])[:,1], label=name*"1")
	#plot(abs.(truesol[2:end,:]-solutions[i][2:end,:])[:,2], label=name*"2")
end
#plot(-im*truesol, label="truesol")
title("Testing Linear Operator: Absolute Errors")
legend()