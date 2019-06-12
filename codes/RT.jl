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

function floatRT(N::Int64)
    #pick from uniform distribution(-N,N)
    C = (rand(3) .- .5)*2*N
    x = rand();
    x<=1/3 ? C[1] = -(C[2]+C[3]) :
    x>2/3 ? C[3] = -(C[1] + C[2]) : C[2] = -(C[1] + C[3])
    return C
end 

#tendency for a resonant triad with ω = [ω₁, ω₂, ω₃] frequencies and ϵ=slow time scale << 1.
@inline function tendRT(z::Array{T,1}, tend::Array{T,1}; ω, ϵ, C) where T<:ComplexF64
    for i = 1 : 3
        tend[i] = im*ω[i]*z[i] + ϵ*C[i]*prod(conj.(z[1:end .!=i]))
    end
    return tend
end

#Runge-Kutta 4(explicit)
@inline function RK4_step(h::Float64, z::Array{T,1}, tend::Array{T,1};ω, ϵ, C) where T<:ComplexF64
    yn = deepcopy(z);
    tendRT(yn, tend, ω=ω, ϵ=ϵ, C=C);
    k = h * tend; #k1
    yn += 1/6*k;
    tendRT(yn + .5*k, tend, ω=ω, ϵ=ϵ, C=C);
    k = h * tend; #k2
    yn += 1/3*k;
    tendRT(yn + .5*k, tend, ω=ω, ϵ=ϵ, C=C);
    k = h * tend; #k3
    yn += 1/3*k;
    tendRT(yn + k, tend, ω=ω, ϵ=ϵ, C=C);
    k = h * tend; #k4
    return yn + (1/6)*k
end

@inline function IFE_step(h::Float64, z::Array{T,1}, update::Array{T,1};ω, ϵ, C) where T<:ComplexF64
    for i = 1 : 3
        update[i] = exp(im*ω[i]*h) * (z[i] + ϵ*C[i]*h*prod(conj.(z[1:end .!=i])))
    end
    return update
end

@inline function phi(A::T) where T<:Complex
    return expm1(A)/A
end

@inline function ETD1_step(h::Float64, z::Array{T,1}, update::Array{T,1};ω, ϵ, C) where T<:ComplexF64
    for i = 1 : 3
        update[i] = exp(im*ω[i]*h)*z[i] + ϵ*C[i]*h*phi(im*ω[i]*h)*prod(conj.(z[1:end .!=i]))
    end
    return update
end

@inline function EUimex_step(h::Float64, z::Array{T,1}, RHS::Array{T,1};ω, ϵ, C) where T<:ComplexF64
    for i = 1 : 3
        RHS[i] = z[i] + h*ϵ*C[i]*prod(conj.(z[1:end .!=i]))
    end
    for i = 1 : 3
        z[i] = RHS[i]/(1-im*h*ω[i])
    end
    return z
end

@inline function CNimex_step(h::Float64, z::Array{T,1}, RHS::Array{T,1};ω, ϵ, C) where T<:ComplexF64
    for i = 1 : 3
        RHS[i] = (1+im*ω[i]*h/2)*z[i] + h*ϵ*C[i]*prod(conj.(z[1:end .!=i]))
    end
    for i = 1 : 3
        z[i] = RHS[i]/(1-im*h*ω[i]/2)
    end
    return z
end

@inline function newtxt!(zAmp::Array{T,1}; name::String="zAmp") where T<:Float64
    writedlm("./"*name*".txt", zAmp')
end

@inline function addtxt!(zAmp::Array{T,1}; name::String="zAmp") where T<:Float64
    open("./"*name*".txt", "a") do io
        writedlm(io, zAmp')
    end
end

function RT_amp(N::Int, h::Float64, every::Int, IC::Array{ComplexF64,1}; 
    ω, ϵ, C, stepper::Function=RK4_step, name::String="zAmp")
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
    newtxt!(z[1], name=name*string(1))
    newtxt!(z[2], name=name*string(2))
    newtxt!(z[3], name=name*string(3))
    for i=2:N+1
        z = stepper(h, z, tend, ω=ω, ϵ=ϵ, C=C)
        if rem(i, every) ==1
            #@printf("%d\t",div(i, every))
            addtxt!(z[1], name=name*string(1))
            addtxt!(z[1], name=name*string(1))
            addtxt!(z[1], name=name*string(1))
        end
    end
end

function plot_RT(name::String="z")
    return 0
end
