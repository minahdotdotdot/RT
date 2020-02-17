include("MMT.jl")
scheme="IFRK3R"; deg = 8;
# time-step, ND final time, save "every"
h=0.01
name=scheme*"-"*string(Int(h*1000000),pad=6)*"-d"*string(deg)
T =10.0#000
M = ceil(Int, T/h);
#T = floor(Int,M*h)
every = 50;#floor(Int, M/1000) # save solution at only 10 time locations.

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
        newtxt!(E[2:Int(end/2)], name=name)
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
hdict = Dict{String, Float64}();
mdict = Dict{String, String}();
cdict = Dict{String, Symbol}();
mds = ["IFRK3" "ETDRK3" "ARK3" "ARK4"]
cs = [:red; :orange; :green; :blue]
hs = [0.06 0.05 0.025 0.01 0.005 0.0025 0.001 0.0005 0.00025 0.0001 0.000075 0.00005]
for (i,m) in enumerate(mds)
    for h in hs
        push!(hdict, m*"-"*string(Int(h*1000000),pad=6) => h)
        push!(mdict, m*"-"*string(Int(h*1000000),pad=6) => m)
        push!(cdict, m*"-"*string(Int(h*1000000),pad=6) => cs[i])
    end
end

ddict=Dict{String, Int64}();
cs2 = [:purple, :pink, :grey]
hs2 = [0.06 0.05 0.025 0.01 0.005 0.0025 0.001 0.0005 0.00025 0.0001 0.000075 0.00005]
ds = [4 6 8]
for h in hs
    for (i,d) in enumerate(ds)
        push!(hdict, "IFRK3R-"*string(Int(h*1000000),pad=6)*"-d"*string(d) => h)
        push!(mdict, "IFRK3R-"*string(Int(h*1000000),pad=6)*"-d"*string(d) => "IFRK3R")
        push!(cdict, "IFRK3R-"*string(Int(h*1000000),pad=6)*"-d"*string(d) => cs2[i])
        push!(ddict, "IFRK3R-"*string(Int(h*1000000),pad=6)*"-d"*string(d) => ds[i])        
    end
end

IFlist = ["IFRK3-060000", "IFRK3-050000","IFRK3-025000", "IFRK3-010000",
"IFRK3-005000","IFRK3-002500", "IFRK3-001000", "IFRK3-000500",
"IFRK3-000250", "IFRK3-000100"];
ETDlist = ["ETDRK3-002500","ETDRK3-001000", "ETDRK3-000500", "ETDRK3-000250", 
"ETDRK3-000100"];
ARK3list = ["ARK3-025000", "ARK3-010000", "ARK3-005000","ARK3-002500", 
"ARK3-001000", "ARK3-001000", "ARK3-000500", "ARK3-000250",
"ARK3-000100", "ARK3-000075", "ARK3-000050"];
ARK4list = ["ARK4-060000", "ARK4-050000", "ARK4-025000", "ARK4-010000", 
"ARK4-005000", "ARK4-002500", "ARK4-001000", "ARK4-000500", 
"ARK4-000250"];

listn = vcat(IFlist, ETDlist, ARK3list, ARK4list);

function saveerror(k, listn::Vector{String}, tru::String, nt::T, rel::Bool=true) where T<:Real
    k = k[2:Int(end/2)];
    truth= log.(readdlm("../txtfiles/"*tru*".txt")' ./ k);
    edict = Dict{String, Float64}();
    for (i, na) in enumerate(listn)
        E = log.(readdlm("../txtfiles/"*na*".txt")' ./ k);
        if rel == true
            push!(edict, na => norm(truth-E,nt)/norm(truth,nt))
        else
            push!(edict, na => norm(truth-E,nt))
        end
    end
    return edict
end
#=IFdict = saveerror(k, IFlist, "IFRK3-001000",2);
ETDdict = saveerror(k, ETDlist, "IFRK3-001000",2);
ARK3dict = saveerror(k, ARK3list, "IFRK3-001000",2);
ARK4dict = saveerror(k, ARK4list, "IFRK3-001000",2); =#

function plotEnergy!(k, name::Vector{String}, m,n,hdict, mdict;#, m::Int, n::Int
    tru::String="IFRK3-025000")
    k = k[2:Int(end/2)]
    truth=readdlm("../txtfiles/"*tru*".txt")' ./ k
    β = truth[65]*64;
    fig, ax = subplots(m,n)
    loglog(k, β ./k, label="1/k+C",c=:black)
    for (i,na) in enumerate(name)
        #subplot(m*100+n*10+i)
        E = readdlm("../txtfiles/"*na*".txt")' ./ k
        err = norm(truth-E,2)/norm(truth,2)
        loglog(k, abs.(β ./k - E))
        #loglog(k, abs.(truth-readdlm("../txtfiles/"*na*".txt")') ./ truth)#, label=na)
        if i > (m-1)*n
            xlabel("Wave Number")
        end
        if mod(i,n)==1
            ylabel("n(k)")
        end
        legend()
        title(na)
    end
    suptitle("Integrating Factor Methods")#mdict(na)*hdict(na))
end

function plotErrvH!(k, name::Vector{String}, hdict, mdict;#, m::Int, n::Int
    tru::String="IFRK3-025000")  
    k = k[2:Int(end/2)]
    truth=log.(readdlm("../txtfiles/"*tru*".txt")' ./ k)
    β = truth[65]*64;
    fig, ax = subplots()#(m,n)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ymin=-1
    ymax=-1
    for (i,na) in enumerate(name)
        #subplot(m*100+n*10+i)
        E = log.(readdlm("../txtfiles/"*na*".txt")' ./ k)
        err = norm(truth-E,2)/norm(truth,2)
        if mdict[na] == "IFRK3R"
            if ddict[na] != 8
                scatter(hdict[na], err, c=cdict[na], s=(ddict[na]-4)*20+20)
            end
        else
            scatter(hdict[na], err, c=cdict[na])
        end
        if i == 1
            ymin = err
            ymax == err
        else
            if ymin > err && na != tru
                ymin = err
            end
            if ymax < err
                ymax = err
            end
        end
        #print(na*": ymin is "*string(ymin)*", and ymax is "*string(ymax)*".\n")
    end 
    for (i,m) in enumerate(mds)
        scatter([], [], c=cs[i], label=m)
    end
    for (i,d) in enumerate(ds)
        scatter([], [], c=cs2[i], label="deg="*string(d), s=(d-4)*20+20)
    end
    ylim(0.5*ymin, ymax*1.5)
    #ylim(5e-3, 10^(0.5))
    xlabel("h: time-step size")
    ylabel("Relative error (2-norm)")
    title("Error in log(Energy Spectrum)")
    legend()
    ax.set_axisbelow(true)
    ax.grid(true)
    # Turn on the minor TICKS, which are required for the minor GRID
    ax.minorticks_on()

    # Customize the major grid
    ax.grid(which="major", linestyle="-", linewidth=0.5, color=:red)
    # Customize the minor grid
    ax.grid(which="minor", linestyle=":", linewidth=0.5, 
        color=:black, alpha=0.7)
end
