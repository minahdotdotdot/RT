include("MMT.jl")
# Problem Parameters
λ = 1   #Defocusing MMT model
α = 1/2
β = 0
F = 0.2
D = [7.51e-25, 7]
fP = funcparams(α, β, λ, F, D)

# Numerical Simulation Parameters
N = 2^13
x = range(0,stop=2*pi-1/N, length=Int(N));
IC = randn(ComplexF64, N)*sqrt(N)/1000; IC[1]=0.0; IC = ifft(IC); 
#IC = cos.(x) + im*sin.(x)
k = vcat(collect(0:N/2), collect(-N/2+1:-1)) # implies domain length is 2π
kind = vcat(collect(Int(N/2)+2:N), collect(1:Int(N/2)+1))
kindnz = vcat(collect(Int(N/2)+2:N), collect(2:Int(N/2))) # indexing w/o zero mode

# Linear operator (depends on k)
L = -im*abs.(k).^fP.α; #L[Int(N/4+2):Int(3*N/4)].=0
L[[6+1, 7+1, 8+1, 9+1, -6+(N+1), -7+(N+1), -8+(N+1), -9+(N+1)]] .= fP.F;
L[2:end] += -196.61 * (abs.(k[2:end]).^(-8)) - fP.D[1]* (abs.(k[2:end]) .^ fP.D[2]); 
L[1]= -200.0;

# Set-up IFRK
A = hcat([0; .5; 0; 0], [0; 0; .5; 0], [0; 0; 0; 1.0])
b = [1/6; 1/3; 1/3; 1/6]
c = [0; 1/2; 1/2; 1];
RK4 = eRKTableau(A, b, c);
A = hcat([0; .5; -1],[0; 0; 2])
b = [1/6; 2/3; 1/6]
c = [0; 1/2; 1];
RK3 = eRKTableau(A, b, c)
RKT = RK4

# time-step, ND final time, save "every"
h = 0.015;
T = 150
M = Int(T/h);
every = 5# Int(M/1000) # save solution at only 1001 time locations.

name = "A"
IFRK!(M, every, IC, h, L, NLfunc, fP, RKT, k, name=name) 
solhat = readCfile(name)#[11:end,:]
#sol = ifft(solhat, 2)
E = transpose(sum(abs.(solhat).^2, dims = 1)/size(solhat)[1])/N^2


using PyPlot, LaTeXStrings
fig, ax = subplots()
#semilogy(k[kindnz], E[kindnz], label="computed")
#semilogy(k[2:Int(end/2)], 1e-24 *(k[2:Int(end/2)]).^(-1/3), label=L"Ck^{-1/3}")
#semilogy(k[2:Int(end/2)], 1e-24 *(k[2:Int(end/2)]).^(-1/2), label=L"Ck^{-1/2}")
loglog(k[2:Int(end/2)], E[2:Int(end/2)], label="computed")
loglog(k[2:Int(end/2)], .24 *(k[2:Int(end/2)]).^(-2), label=L"Ck^{-2}")
loglog(k[2:Int(end/2)], .24 *(k[2:Int(end/2)]).^(-1), label=L"Ck^{-1}")
xlabel("Wave Number")
ylabel("n(k)")
legend()
title(L"t\in[11, 1000]")
savefig(name*"ES.png")
close(fig)



