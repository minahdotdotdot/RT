using LinearAlgebra, PyPlot, Random, DelimitedFiles, Printf

#generating IC on unit circle
function onUnitCircle(n::Int=3)
    z = Array{ComplexF64,1}(undef, n);
    for i = 1 : n
        x = 2*(rand()-0.5)
        if abs(x) > 1
            error("x must be in [-1,1].")
        else
            y = rand()
            if y >=0.5
                z[i] = x + im * sqrt(1-x^2) 
            else
                z[i] = x - im * sqrt(1-x^2) 
            end
        end
    end
    return z
end

#Finds 3 integers that add to 0, where N is the largest possible magnitude.
function intRT(N::Int64)
    ω = randperm(N)[1:3]
    x = rand(3)
    for i = 1 : 3  
        x[i]>0.5 ? ω[i]=-ω[i] : ω[i]=ω[i]
    end
    x = rand();
    x<=1/3 ? ω[1] = - (ω[2]+ω[3]) :
    x>2/3 ? ω[3] = - (ω[1] + ω[2]) : ω[2] = - (ω[1] + ω[3]) 
    return ω
end 

#Finds 3 floats that add to 0.
function floatRT(N::Int64)
    #pick  3 floats from uniform distribution(-N,N)
    C = (rand(3) .- .5)*2*N
    x = rand();
    # change the first float  if x <=1/3    to meet C1+C2+C3=0
    # change the second float if 1/3<=x<2/3 to meet C1+C2+C3=0
    # change the third float  if x>2/3      to meet C1+C2+C3 = 0
    x<=1/3 ? C[1] = -(C[2]+C[3]) : #
    x>2/3 ? C[3] = -(C[1] + C[2]) : C[2] = -(C[1] + C[3])
    return C
end 

@inline function Lfunc(z, ω, h)
    return exp.(im*h*ω) .* z
end

@inline function NLfunc(z, ϵ, C)
    N = deepcopy(z)
    for i = 1 : 3
        N[i] = ϵ*C[i]*prod(conj.(z[1:end .!=i]))
    end
    return N
end

#tendency for a resonant triad with ω = [ω₁, ω₂, ω₃] frequencies and ϵ=slow time scale << 1.
@inline function tendRT(z::Array{T,1}, tend::Array{T,1}; ω, ϵ, C) where T<:ComplexF64
    return im*ω .* z + NLfunc(z, ϵ, C)
end

#Runge-Kutta 4(explicit)
@inline function RK4(h::Float64, z::Array{T,1}, tend::Array{T,1};ω, ϵ, C) where T<:ComplexF64
    yn = deepcopy(z);
    tendRT(yn, tend, ω=ω, ϵ=ϵ, C=C);
    k = h * tend; #k1
    yn += 1/6*k;
    tendRT(z + .5*k, tend, ω=ω, ϵ=ϵ, C=C);
    k = h * tend; #k2
    yn += 1/3*k;
    tendRT(z + .5*k, tend, ω=ω, ϵ=ϵ, C=C);
    k = h * tend; #k3
    yn += 1/3*k;
    tendRT(z + k, tend, ω=ω, ϵ=ϵ, C=C);
    k = h * tend; #k4
    return yn + (1/6)*k
end

@inline function IFE(h::Float64, z::Array{T,1}, update::Array{T,1};ω, ϵ, C) where T<:ComplexF64
    return exp.(im*h*ω) .* (z + NLfunc(z, ϵ, C))
end


@inline function phi(A::T) where T<:Complex
    return expm1(A)/A
end

@inline function ETD1(h::Float64, z::Array{T,1}, update::Array{T,1};ω, ϵ, C) where T<:ComplexF64
    return Lfunc(z, ω, h) + h*NLfunc(z, ϵ, C) .* phi.(im*h*ω)
end

@inline function EUimex(h::Float64, z::Array{T,1}, RHS::Array{T,1};ω, ϵ, C) where T<:ComplexF64
    RHS = z + h*NLfunc(z, ϵ, C)
    return RHS ./ ( 1 .- im*h*ω)
end

@inline function CNimex(h::Float64, z::Array{T,1}, RHS::Array{T,1};ω, ϵ, C) where T<:ComplexF64
    RHS = (1 .+ im*ω*h/2) .*z + h*NLfunc(z, ϵ, C)
    return RHS ./ (1 .- im*h/2*ω)
end

@inline function IFRK2(h::Float64, z::Array{T,1}, RHS::Array{T,1};ω, ϵ, C) where T<:ComplexF64
    yn = deepcopy(z)
    k = NLfunc(z, ϵ, C)       # k1
    yn += h/2*k
    k = NLfunc(z + h*k, ϵ, C) # k2
    return exp.(im*h*ω) .* (yn .+ h/2*k)
end


@inline function IFRK4(h::Float64, z::Array{T,1}, RHS::Array{T,1};ω, ϵ, C) where T<:ComplexF64
    yn = deepcopy(z)
    k = NLfunc(z, ϵ, C)         # k1
    yn += h/6*k
    k = NLfunc(z + h/2*k, ϵ, C) # k2
    yn += h/3*k
    k = NLfunc(z + h/2*k, ϵ, C) # k3
    yn += h/3*k
    k = NLfunc(z + h*k, ϵ, C)   # k4
    return  exp.(im*h*ω) .* (yn + h/6*k)
end

@inline function newtxt!(zAmp::Array{T,1}; name::String="zAmp") where T<:Float64
    writedlm("../txtfiles/"*name*".txt", zAmp')
end

@inline function addtxt!(zAmp::Array{T,1}; name::String="zAmp") where T<:Float64
    open("../txtfiles/"*name*".txt", "a") do io
        writedlm(io, zAmp')
    end
end

@inline function newtxt!(zAmp::Array{T,1}; name::String="zAmp") where T<:Complex{Float64}
    writedlm("../txtfiles/"*name*"_Re.txt", real.(zAmp'))
    writedlm("../txtfiles/"*name*"_Im.txt", imag.(zAmp'))
end

@inline function addtxt!(zAmp::Array{T,1}; name::String="zAmp") where T<:Complex{Float64}
    open("../txtfiles/"*name*"_Re.txt", "a") do io
        writedlm(io, real.(zAmp'))
    end
    open("../txtfiles/"*name*"_Im.txt", "a") do io
        writedlm(io, imag.(zAmp'))
    end
end

@inline function readCfile(name::String)
    return readdlm("../txtfiles/"*name*"_Re.txt") + im * readdlm("../txtfiles/"*name*"_Im.txt")
end

function RT_amp(N::Int, h::Float64, every::Int, IC::Array{ComplexF64,1}; 
    ω, ϵ, C, stepper::Function=RK4, name::String="zAmp")
    z = deepcopy(IC)
    tend = deepcopy(IC)
    newtxt!(abs.(z), name=name)
    if every != 1
        for i=2:N+1
            z = stepper(h, z, tend, ω=ω, ϵ=ϵ, C=C)
            if rem(i, every) ==1
                #@printf("%d\t",div(i, every))
                addtxt!(abs.(z), name=name)
            end
        end
    else
        for i=2:N+1
            z = stepper(h, z, tend, ω=ω, ϵ=ϵ, C=C)
            addtxt!(abs.(z), name=name)
        end
    end
end

function RT(N::Int, h::Float64, every::Int, z::Array{ComplexF64,1}; 
    ω, ϵ, C, stepper::Function=RK4_step, name::String="z")
    tend = deepcopy(z)
    newtxt!(z, name=name)
    for i=2:N+1
        z = stepper(h, z, tend, ω=ω, ϵ=ϵ, C=C)
        print(z)
        if rem(i, every) ==1
            addtxt!(z, name=name)
        end
    end
end

function plot_RT(name::String="z")
    return 0
end
