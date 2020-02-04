include("MMT.jl")
scheme="IFRK3"
# time-step, ND final time, save "every"
h=0.01
name=scheme*"-"*string(Int(h*1000000),pad=6)
T =10000
M = ceil(Int, T/h);
#T = floor(Int,M*h)
every = floor(Int, M/1000) # save solution at only 10 time locations.

# Problem Parameters
λ = 1;  #Defocusing MMT model
α = 1/2;
β = 0;
F = 0.05;
D = [2.51e-57, 16];
fP = funcparams(α, β, λ, F, D);

# Numerical Simulation Parameters
N = 2^13;
k = vcat(collect(0:N/2), collect(-N/2+1:-1)); # implies domain length is 2π
kind = vcat(collect(Int(N/2)+2:N), collect(1:Int(N/2)+1));
kindnz = vcat(collect(Int(N/2)+2:N), collect(2:Int(N/2))); # indexing w/o zero mode

#IC
#IC = cos.(range(0,2*pi,length=N)) + im*sin.(range(0,2*pi,length=N))
IC = randn(ComplexF64, N)*sqrt(N)/1000; IC[1]=0.0; IC = ifft(IC);

# Linear operator (depends on k)
L = -im*abs.(k).^fP.α; #L[Int(N/4+2):Int(3*N/4)].=0
L[[6+1, 7+1, 8+1, 9+1, -6+(N+1), -7+(N+1), -8+(N+1), -9+(N+1)]] .+= fP.F;
L[2:end] += -196.61 * (abs.(k[2:end]).^(-8)) - fP.D[1]* (abs.(k[2:end]) .^ fP.D[2]); 
L[1]= -200.0;

if scheme ∈ ["IFRK3_rat", "IFRK4_rat"]
	using MAT
	include("IF_methods.jl")
	file = matopen("../data/Lhc_"*name*".mat","w")
	write(file, "scheme", scheme)
	write(file, "h", h)
	close(file)
	if scheme == "IFRK3_rat"
		file = matopen("../data/"*scheme*"h="*string(h)*".mat", "w")
		write(file, "L", L)
		write(file, "x", IFRK3.x)
		write(file, "crat", IFRK3.c)
		close(file)
	else
		file = matopen("../data/Lhc.mat", "w")
		write(file, "L", L)
		write(file, "h", h)
		write(file, "x", IFRK4.x)
		write(file, "crat", IFRK4.c)
		close(file)
	end

end

function runMMT(method::String, 
	M::Int, every::Int, IC, h, L, NLfunc, fP, k, name, cont::Bool=false)
	if method ∈ ["ETDRK2", "ETDRK3", "ETDRK4", "ETDRK4B"]
		include("ETD_methods.jl")
		ETDRK!(M, every, IC, h, L, NLfunc, fP, ETDdict[method], k, name=name, cont=cont)
	elseif method ∈ ["IFRK3", "IFRK4"]
		include("IF_methods.jl")
		IFRK!(M, every, IC, h, L, NLfunc, fP, IFdict[method], k, name=name, cont=cont)
	elseif method ∈ ["ARK3", "ARK4"]
		include("IMEX_methods.jl")
		IMEXRK!(M, every, IC, h, L, NLfunc, fP, IMEXdict[method], k, name=name, cont=cont)
	else
		error("method must be ETD, IF, or IMEX.")
	end
end

function runMMT(method::eRKTableau, 
	M::Int, every::Int, IC, h, L, NLfunc, fP, k, name, cont::Bool=false)
	IFRK!(M, every, IC, h, L, NLfunc, fP, method, k; name=name, cont=cont)
end

function saveEnergy!(k, N, T, name::String; scheme::String, h, ES::Bool=true)
	solhat = readCfile(name)[11:end,:]
	#sol = ifft(solhat, 2)
	E = k .* transpose(sum(abs.(solhat).^2, dims = 1)/size(solhat)[1])/N^2;
	if ES == false
		fig, ax = subplots()
		loglog(k[2:Int(end/2)], E[2:Int(end/2)], label=L"k\times"*"computed")
		xlabel("Wave Number")
		ylabel("n(k)")
		legend()
		title("h="*string(h)*", "*scheme)
		savefig(scheme*"-"*string(Int(h*1000000),pad=6)*"ES.png")
		close(fig)
	else
		newtxt!(E[2:Int(end/2)], name)
	end
end
