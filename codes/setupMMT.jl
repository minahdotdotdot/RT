 include("MMT.jl")
#=scheme="IFRK3R"; deg = 8;
# time-step, ND final time, save "every"
h=0.025
name=scheme*"-"*string(Int(h*1000000),pad=6)*"-d"*string(deg)
=#
scheme="IFRK3"; deg = 0;
h=0.0001
name=scheme*"-"*string(Int(h*1_000_000),pad=6);
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

import Base: +

function +(c::Matrix{Vector{Complex{T}}}, x::T) where T<:AbstractFloat
    m,n = size(c);
    for i = 1 : m
        for j = 1 : n
            c[i,j] = c[i,j] .+ x
        end
    end
    return c
end

function +(c::Matrix{Vector{Complex{T}}}, x::Vector{Complex{T}}) where T<:AbstractFloat
    m,n = size(c);
    for i = 1 : m
        for j = 1 : n
            c[i,j] = c[i,j] + x
        end
    end
    return c
end

function runMMT(method::String, 
    M::Int, every::Int, IC, h, L, NLfunc, fP, k, name, cont::Bool=false; deg::Int=0)
    if method ∈ ["ETDRK2", "ETDRK3", "ETDRK4", "ETDRK4B"]
        include("ETD_methods.jl")
        ETDRK!(M, every, IC, h, L, NLfunc, fP, ETDdict[method], k, name=name, cont=cont)
    elseif method ∈ ["IFRK3", "IFRK4"]
        include("IF_methods.jl")
        #x = 5e-14*rand(length(L)) .* exp.(im*2*pi*rand(length(L)))
        #x = (5e-14*(sign.(randn(length(L))).*rand(eltype(L),length(L))) .+5e-15)
        #x = 5e-15 .+ 5e-14*exp.(im * range(0, 2pi, length=length(L)+1)[1:end-1])
        RKT = eRKTableau(IFdict[method].A, 
            IFdict[method].b, 
            IFdict[method].c,#+ x,
            IFdict[method].x)
        IFRK!(M, every, IC, h, L, NLfunc, fP, RKT, k, name=name, cont=cont)
    elseif method ∈ ["ARK3", "ARK4"]
        include("IMEX_methods.jl")
        IMEXRK!(M, every, IC, h, L, NLfunc, fP, IMEXdict[method], k, name=name, cont=cont)
    elseif method ∈ ["IFRK3R", "IFRK4R"]
        include("IF_methods.jl");
        file = matopen("../data/"*scheme*"h="*string(h)*"d"*string(deg)*"R.mat");
        crat = read(file, "crat");
        close(file)
        
        if method == "IFRK3R"
            #err = IFRK3.c-crat;
            crat = convertcrat(Vector{ComplexF64}, crat, IFRK3.c)
            RKT = eRKTableau(IFRK3.A, IFRK3.b, crat, IFRK3.x)
        else
            crat = convertcrat(Vector{ComplexF64}, crat, IFRK4.c)
            RKT = eRKTableau(IFRK4.A, IFRK4.b, crat, IFRK4.x)
        end
        IFRK!(M, every, IC, h, L, NLfunc, fP, RKT, k, name=name, cont=cont)
    end
end

function convertcrat(dt::DataType, a::Matrix{Any}, c)
    tmp = copy(c)
    for i = 1 : size(c)[1]
        for j = 1:size(c)[2]
            tmp[i,j] = dt(a[i,j][:,1])
        end
    end
    return tmp
end


#=
function runMMT(method::eRKTableau, 
    M::Int, every::Int, IC, h, L, NLfunc, fP, k, name, cont::Bool=false)
    IFRK!(M, every, IC, h, L, NLfunc, fP, method, k, name=name, cont=cont)
end
=#

function saveEnergy!(k, N, T, name::String; 
    scheme::String, h, ES::Bool=true, deg::Int=0)
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
        if deg == 0
            savefig(name*"-ES.png")
                #scheme*"-"*string(Int(h*1000000),pad=6)*"ES.png")
            close(fig)
        else
            savefig(name*"-ES.png")
            #scheme*"-"*string(Int(h*1000000),pad=6)*"d"*string(deg)*"ES.png")
            close(fig)
        end
    else
        newtxt!(E[2:Int(end/2)], name=name*"avg")
    end
end


#=
if scheme ∈ ["IFRK3R", "IFRK4R"]
    include("IF_methods.jl")
    file = matopen("../data/Lhc_"*scheme*"h="*string(h)*"d"*string(deg)*".mat","w")
    write(file, "scheme", scheme)
    write(file, "h", h)
    write(file, "deg", deg)
    close(file)
    if scheme == "IFRK3R"
        file = matopen("../data/"*scheme*"h="*string(h)*"d"*string(deg)*".mat", "w")
        write(file, "L", L)
        write(file, "x", IFRK3.x)
        write(file, "crat", IFRK3.c)
        close(file)
    else
        file = matopen("../data/"*scheme*"h="*string(h)*"d"*string(deg)*".mat", "w")
        write(file, "L", L)
        write(file, "x", IFRK4.x)
        write(file, "crat", IFRK4.c)
        close(file)
    end
end
=#
