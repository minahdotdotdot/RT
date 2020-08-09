using PyPlot, DelimitedFiles, LinearAlgebra, DataFrames, GLM, Printf

# Numerical Simulation Parameters
N = 2^13;
k = vcat(collect(0:N/2), collect(-N/2+1:-1)); # implies domain length is 2π
kind = vcat(collect(Int(N/2)+2:N), collect(1:Int(N/2)+1));
kindnz = vcat(collect(Int(N/2)+2:N), collect(2:Int(N/2))); # indexing w/o zero mode

hdict = Dict{String, Float64}();
mdict = Dict{String, String}();
cdict = Dict{String, Symbol}();
mds = ["IFRK3" "ETDRK3" "ARK3" "ARK4"]
cs = [:red; :orange; :green; :blue]
hs = [0.06; 0.05; 0.025; 0.01;
0.005; 0.0025; 0.001; 0.0005;
0.00025; 0.0001; 0.000075; 0.00005]
for (i,m) in enumerate(mds)
    for h in hs
        push!(hdict, m*"-"*string(Int(h*1000000),pad=6) => h)
        push!(mdict, m*"-"*string(Int(h*1000000),pad=6) => m)
        push!(cdict, m*"-"*string(Int(h*1000000),pad=6) => cs[i])
    end
end

ddict=Dict{String, Int64}();
cs2 = [:purple, :pink, :grey]
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
"ARK4-000250", "ARK4-000100"];
IFr4list = ["IFRK3R-060000-d4", "IFRK3R-050000-d4", "IFRK3R-025000-d4", "IFRK3R-010000-d4",
"IFRK3R-005000-d4", "IFRK3R-002500-d4", "IFRK3R-001000-d4", "IFRK3R-000500-d4",
"IFRK3R-000250-d4"];
IFr6list = ["IFRK3R-060000-d6", "IFRK3R-050000-d6", "IFRK3R-025000-d6", "IFRK3R-010000-d6",
"IFRK3R-005000-d6", "IFRK3R-002500-d6", "IFRK3R-001000-d6", "IFRK3R-000500-d6",
"IFRK3R-000250-d6"];
IFr8list = ["IFRK3R-060000-d8", "IFRK3R-050000-d8", "IFRK3R-025000-d8", "IFRK3R-010000-d8",
"IFRK3R-005000-d8", "IFRK3R-002500-d8", "IFRK3R-001000-d8", "IFRK3R-000500-d8",
"IFRK3R-000250-d8"];

listn = vcat(ETDlist, ARK3list, ARK4list, IFr8list, IFr6list, IFlist[1:end-1], IFr4list);

function saveerror(k, listn::Vector{String}, tru::String, nt::T, rel::Bool=true;
    errtype::String="Dict", hdict=0,logB::Bool=true) where T<:Real
    k = k[2:Int(end/2)];
    truth= log.(readdlm("../txtfiles/"*tru*".txt")' ./ k);
    if errtype=="Dict"
        edict = Dict{String, Float64}();
        for (i, na) in enumerate(listn)
            E = log.(readdlm("../txtfiles/"*na*".txt")' ./ k);
            if rel == true
                if logB == false    
                    push!(edict, na => norm(truth-E,nt)/norm(truth,nt))
                else
                    push!(edict, na => log(norm(truth-E,nt)/norm(truth,nt)))
                end
            else
                if logB == false
                    push!(edict, na => norm(truth-E,nt))
                else
                    push!(edict, na => log(norm(truth-E,nt)))
                end
            end
        end
    elseif errtype=="Matrix"
        err = Matrix{Float64}(undef, length(listn),2);
        for (i, na) in enumerate(listn)
            E = log.(readdlm("../txtfiles/"*na*".txt")' ./ k);
            err[i,1] = hdict[na]
            if rel==true
                err[i,2] = norm(truth-E,nt)/norm(truth,nt)
            else
                err[i,2] = norm(truth-E,nt)
            end
        end
        p = sortperm(err[:,1]);
        err = err[p,:]
        if logB == false
            return err
        else
            return log.(err)
        end
    end
end
#=IFdict = saveerror(k, IFlist, "IFRK3-000100",2);
ETDdict = saveerror(k, ETDlist, "IFRK3-000100",2);
ARK3dict = saveerror(k, ARK3list, "IFRK3-000100",2);
ARK4dict = saveerror(k, ARK4list, "IFRK3-000100",2);
IFr4dict = saveerror(k, IFr4list, "IFRK3-000100",2);
IFr6dict = saveerror(k, IFr6list, "IFRK3-000100",2);
IFr8dict = saveerror(k, IFr8list, "IFRK3-000100",2);=#
IFerr = saveerror(k, IFlist, "IFRK3-000100",2, true, errtype="Matrix", hdict=hdict);
ETDerr = saveerror(k, ETDlist, "IFRK3-000100",2, true, errtype="Matrix", hdict=hdict);
A3err = saveerror(k, ARK3list, "IFRK3-000100",2, true, errtype="Matrix", hdict=hdict);
A4err = saveerror(k, ARK4list, "IFRK3-000100",2, true, errtype="Matrix", hdict=hdict);
IFr4err = saveerror(k, IFr4list, "IFRK3-000100",2, true, errtype="Matrix", hdict=hdict);
IFr6err = saveerror(k, IFr6list, "IFRK3-000100",2, true, errtype="Matrix", hdict=hdict);
IFr8err = saveerror(k, IFr8list, "IFRK3-000100",2, true, errtype="Matrix", hdict=hdict);

function plotErrRationals!(k, hs::Vector{T}; nt::Float64=2.0) where T<:AbstractFloat
    fig, ax = subplots()
    ax.set_xscale("log")
    ax.set_yscale("log")
    k2 = k[2:Int(end/2)]
    truth = log.(readdlm("../txtfiles/IFRK3-000100.txt")' ./ k2)
    for h in hs
        #truth = log.(readdlm("../txtfiles/IFRK3-"*string(Int(h*1000000),pad=6)*".txt")' ./ k2)
        for i = 3 : -1 : 1
            E = log.(readdlm("../txtfiles/IFRK3R-"*string(Int(h*1000000),pad=6)*"-d"*string(2*i+2)*".txt")' ./ k2)
            err = norm(truth-E,nt)/norm(truth,nt)
            scatter(h, err, c=cs2[i], s=40*i-20)
        end
    end
    for i = 1 : 3
        scatter([], [], c=cs2[i], s=40*i-20, label="deg="*string(2*i+2))
    end
    xlabel("h: time-step size")
    ylabel("Relative error ("*string(nt)*"-norm)")
    title("Error in log(ES) of rational approx against IFRK3, h=0.0001")
    legend()
    xlim(0.5*minimum(hs), 1.5*maximum(hs))
    ylim(2e-3, 1e-2)
    ax.set_axisbelow(true)
    ax.grid(true)
    # Turn on the minor TICKS, which are required for the minor GRID
    ax.minorticks_on()

    # Customize the major grid
    ax.grid(which="major", axis="both", linestyle="-", linewidth=0.5, color=:red)
    # Customize the minor grid
    ax.grid(which="minor", axis="both", linestyle=":", linewidth=0.5, 
        color=:black, alpha=0.7)

end
fitdict=Dict{String,Float64}();
function plotErrvH!(k, name::Vector{String}, hdict, mdict;#, m::Int, n::Int
    tru::String="IFRK3-000100", regfit::Bool=false, minV=0, nt::Float64=2.0)  
    k2 = k[2:Int(end/2)]
    truth=log.(readdlm("../txtfiles/"*tru*".txt")' ./ k2)
    #fig, ax = subplots()#(m,n); #ymin=-1; ymax=-1
    for (i,na) in enumerate(name)
        E = log.(readdlm("../txtfiles/"*na*".txt")' ./ k2)
        err = norm(truth-E,nt)/norm(truth,nt)
        if mdict[na] == "IFRK3R"
            scatter(hdict[na], err, c=cdict[na], s=(ddict[na]-4)*20+20)
        else
            scatter(hdict[na], err, c=cdict[na])
        end
        #=if i == 1
            ymin = err
            ymax == err
        else
            if ymin > err && na != tru
                ymin = err
            end
            if ymax < err
                ymax = err
            end
        end=#
    end 
    if regfit == true
        err = saveerror(k, name, tru,nt, true, errtype="Matrix", hdict=hdict);
        data = DataFrame(logh=err[:,1], logerr=err[:,2]);
        ols=lm(@formula(logh ~ logerr), data);
        push!(fitdict, mdict[name[1]] => ols.model.pp.beta0[1])
        errrange = range(minV, stop=maximum(err), length=1001)
        fitline = exp.(predict(ols, DataFrame(logerr=errrange)))
        #axvline(fitline[1], c=cdict[name[1]])
        plot(fitline, exp.(errrange), c=cdict[name[1]], linewidth=0.75) 
        scatter(fitline[1], exp.(errrange)[1], c=cdict[name[1]], marker="x") 
    end
end

function plotErrvH2!(nt::Float64=2.)
    fig,ax = subplots()
    err = saveerror(k, IFlist[1:end-1], "IFRK3-000100", Inf, true, errtype="Matrix", hdict=hdict);
    minV=minimum(err[:,2]);
    #print(minV)
    ax.set_xscale("log")
    ax.set_yscale("log")
    axhline(y=exp(minV),c=:red, label="Minimum Error")
    plotErrvH!(k,ETDlist, hdict, mdict, regfit=true, minV=minV, nt=Inf)
    plotErrvH!(k,ARK3list, hdict, mdict, regfit=true, minV=minV, nt=Inf)
    plotErrvH!(k,ARK4list, hdict, mdict, regfit=true, minV=minV, nt=Inf)
    plotErrvH!(k,IFr8list, hdict, mdict, regfit=false, nt=Inf)
    plotErrvH!(k,IFr6list, hdict, mdict, regfit=false, nt=Inf)
    plotErrvH!(k,IFlist[1:end-1], hdict, mdict, regfit=false, nt=Inf)
    plotErrvH!(k,IFr4list, hdict, mdict, regfit=false, nt=Inf)
    for (i,m) in enumerate(mds)
        if m == "IFRK3"
            scatter([], [], c=cs[i], label=m)
        else
            scatter([], [], c=cs[i], label=m*", "*@sprintf("m≈%.2f",fitdict[m]))
        end
    end
    for (i,d) in enumerate(ds)
        scatter([], [], c=cs2[i], label="deg="*string(d), s=(d-4)*20+20)
    end
    plot([], [], c=:black, marker="x", label="fitted")
    xlabel("h: time-step size")
    ylabel("Relative error ("*string(nt)*"-norm)")
    title("Error in log(Energy Spectrum)")
    legend()
    ax.set_axisbelow(true)
    ax.grid(true)
    #ax.grid(which="major", axis="both", linestyle="--")
    #ax.grid(which="major", axis="both", linestyle="--")
    # Turn on the minor TICKS, which are required for the minor GRID
    ax.minorticks_on()

    # Customize the major grid
    ax.grid(which="major", axis="both", linestyle="-", linewidth=0.5, color=:red)
    # Customize the minor grid
    ax.grid(which="minor", axis="both", linestyle=":", linewidth=0.5, 
        color=:black, alpha=0.7)
    axvline.(vcat((2:2:10)*1e-8, (2:2:10)*1e-7, (2:2:10)*1e-6, (2:2:10)*1e-5, 
        (2:2:10)*1e-4, (2:2:10)*1e-3, (2:2:10)*1e-2, (2:2:10)*1e-9), 
    linestyle=":", linewidth=0.5, color=:black, alpha=0.7)
    ylim(1e-3, 3e1)
    xlim(7e-10, 0.1)
end
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

