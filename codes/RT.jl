using LinearAlgebra, PyPlot, Random, DelimitedFiles, Printf, SparseArrays
include("readwrite.jl")
include("IMEX_methods.jl")
#generating IC on unit circle
function onUnitCircle(n::Int=3)
    z = Vector{ComplexF64}(undef, n);
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
    # change the third float  if x>2/3      to meet C1+C2+C3=0
    x<=1/3 ? C[1] = -(C[2]+C[3]) : #
    x>2/3 ? C[3] = -(C[1] + C[2]) : C[2] = -(C[1] + C[3])
    return C
end 

@inline function Lfunc(z::Vector{ComplexF64}, ω::Vector{Int}, h::Float64)
    return exp.(im*h*ω) .* z
end

struct NLfuncparams
    ϵ::Float64
    C::Array{Float64,1}
    NLfuncparams(ϵ,C) = new(copy(ϵ), copy(C))
end

@inline function NLfunc(z::Array{ComplexF64,1}, nlfp::NLfuncparams)
    N = deepcopy(z)
    for i = 1 : 3
        N[i] = prod(conj.(z[1:end .!=i]))
    end
    return nlfp.ϵ * nlfp.C .* N
end

@inline function NLfunc(z, ϵ, C)
    N = deepcopy(z)
    for i = 1 : 3
        N[i] = ϵ*C[i]*prod(conj.(z[1:end .!=i]))
    end
    return N
end


#tendency for a resonant triad with ω = [ω₁, ω₂, ω₃] frequencies and ϵ=slow time scale << 1.
@inline function tendRT!(z::Array{T,1}; ω, ϵ, C) where T<:ComplexF64
    return im*ω .* z + NLfunc(z, ϵ, C)
end

#Runge-Kutta 4(explicit)
@inline function RK4RT(h::Float64, z::Array{T,1}, tend::Array{T,1}; ω, ϵ, C) where T<:ComplexF64
    yn = deepcopy(z)
    tend = tendRT!(yn, ω=ω, ϵ=ϵ, C=C);
    k = h*tend; #k1
    yn += 1/6*k;
    tend = tendRT!(z + (.5*k), ω=ω, ϵ=ϵ, C=C);
    k = h*tend; #k2
    yn += 1/3*k;
    tend = tendRT!(z + (.5*k), ω=ω, ϵ=ϵ, C=C);
    k = h*tend; #k3
    yn += 1/3*k;
    tend = tendRT!(z + k, ω=ω, ϵ=ϵ, C=C); #tend=k4
    return yn + (1/6)*h*tend
end

@inline function IFE(h::Float64, z::Array{T,1}, update::Array{T,1};ω, ϵ, C) where T<:ComplexF64
    return exp.(im*h*ω) .* (z + h*NLfunc(z, ϵ, C))
end


@inline function phi1(A::T) where T<:Complex
    return expm1(A)/A
end
####### ETD methods ###########
using ExponentialUtilities

@inline function ETD1(h::Float64, z::Array{T,1}, update::Array{T,1};ω, ϵ, C) where T<:ComplexF64
    return Lfunc(z, ω, h) + h*NLfunc(z, ϵ, C) .* phi1.(im*h*ω)
end

@inline function ETDRK2(h::Float64, z::Array{T,1}, update::Array{T,1};ω, ϵ, C) where T<:ComplexF64
    phis = phi(Diagonal(h*im*ω),2);    #phis[i] is the (i-1)^th phi function.
    yn = zeros(ComplexF64,3)
    k  = NLfunc(z, ϵ, C)                      # k1
    yn += (phis[2]-phis[3]) *k
    k = NLfunc(phis[1]*z + h*phis[2]*k, ϵ, C) # k2
    yn += phis[3]*k
    return phis[1] * z + (h*yn)
end

@inline function ETDRK3(h::Float64, z::Array{T,1}, update::Array{T,1};ω, ϵ, C) where T<:ComplexF64
    #phis=[1=>exp(hL), 2=>phi_1(hL), 3=>phi_2(hL), 4=>phi_3(hL), 5=>exp(.5hL), 6=>phi_1(.5hL)]
    phis = [phi(Diagonal(h*im*ω),3); phi(Diagonal(h/2*im*ω),1)] 
    yn = zeros(ComplexF64,3)
    k13  = NLfunc(z, ϵ, C)                          # k1
    yn += (phis[2]-3phis[3]+4phis[4]) *k13
    k2 = NLfunc(phis[5]*z + h*.5*phis[6]*k13, ϵ, C) # k2
    yn += (4*phis[3]-8*phis[4])*k2 
    k13 = NLfunc(phis[1]*z+ h*(-phis[2]*k13 + 2*phis[2]*k2), ϵ, C) # k3
    yn += (-phis[3]+4*phis[4])*k13
    return phis[1]* z + (h*yn)
end

@inline function ETDRK4(h::Float64, z::Array{T,1}, update::Array{T,1};ω, ϵ, C) where T<:ComplexF64
    #phis=[1=>exp(hL), 2=>phi_1(hL), 3=>phi_2(hL), 4=>phi_3(hL), 5=>exp(.5hL), 6=>phi_1(.5hL)]
    phis = [phi(Diagonal(h*im*ω),3); phi(Diagonal(h/2*im*ω),1)] 
    yn = zeros(ComplexF64,3)
    k14  = NLfunc(z, ϵ, C)                          # k1
    yn += (phis[2]-3phis[3]+4phis[4]) *k14
    k23 = NLfunc(phis[5]*z + h*.5*phis[6]*k14, ϵ, C) # k2
    yn += (2*phis[3]-4*phis[4])*k23 
    k23 = NLfunc(phis[5]*z + h*.5*phis[6]*k23, ϵ, C) # k3
    yn += (2*phis[3]-4*phis[4])*k23 
    k14 = NLfunc(phis[1]*z+h*((.5*phis[6]*(phis[5]-sparse(Matrix(I,3,3))))*k14 + phis[6]*k23), ϵ, C) # k4
    yn += (-phis[3]+4*phis[4])*k14
    return phis[1]* z + (h*yn)
end

@inline function ETDRK4B(h::Float64, z::Array{T,1}, update::Array{T,1};ω, ϵ, C) where T<:ComplexF64
    #phis=[1=>exp(hL), 2=>phi_1(hL), 3=>phi_2(hL), 4=>phi_3(hL), 5=>exp(.5hL), 6=>phi_1(.5hL), 7=>phi_2(.5hL)]
    phis = [phi(Diagonal(h*im*ω),3); phi(Diagonal(h/2*im*ω),2)] 
    yn = zeros(ComplexF64,3)
    k14  = NLfunc(z, ϵ, C)                          # k1
    yn += (phis[2]-3phis[3]+4phis[4]) *k14
    k23 = NLfunc(phis[5]*z + h*.5*phis[6]*k14, ϵ, C) # k2
    yn += (2*phis[3]-4*phis[4])*k23 
    k23 = NLfunc(phis[5]*z + h*((.5*phis[6]-phis[7])*k14+.5*phis[6]*k23), ϵ, C) # k3
    yn += (2*phis[3]-4*phis[4])*k23 
    k14 = NLfunc(phis[1]*z+h*((phis[2]-2*phis[3])*k14 + 2*phis[3]*k23), ϵ, C) # k4
    yn += (-phis[3]+4*phis[4])*k14
    return phis[1]* z + (h*yn)
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
    yn = zeros(ComplexF64,3)
    k = NLfunc(z, ϵ, C)               # k1
    yn += .5*k
    k = NLfunc(z + h*k, ϵ, C)         # k2
    yn += .5*k
    return exp.(im*h*ω) .* (z + h*yn)
end

@inline function IFRK3(h::Float64, z::Array{T,1}, RHS::Array{T,1}; ω, ϵ, C) where T<:ComplexF64
    yn = zeros(ComplexF64,3)
    k13 = NLfunc(z, ϵ, C)               # k1
    yn += k13/6
    k2 = NLfunc(z + .5*h*k13, ϵ, C)     # k2
    yn += 2. *k2/3
    k13 = NLfunc(z +h*(-k13+2*k2), ϵ, C)  # k3
    yn += k13/6
    return exp.(im*h*ω) .* (z + h*yn)
end


@inline function IFRK4(h::Float64, z::Array{T,1}, RHS::Array{T,1}; ω, ϵ, C) where T<:ComplexF64
    yn = zeros(ComplexF64,3)
    k = NLfunc(z, ϵ, C)               # k1
    yn += k/6
    k = NLfunc(z + .5*h*k, ϵ, C)      # k2
    yn += k/3
    k = NLfunc(z + .5*h*k, ϵ, C)      # k3
    yn += k/3
    k = NLfunc(z + h*k, ϵ, C)         # k4
    yn += k/6
    return exp.(im*h*ω) .* (z + h*yn)
end

#####IMEX methods
# @inline function lincom(Aerow, stages)
#     tmp = Aerow[1] * stages[1,:]
#     for i = 2 : length(Aerow)
#         tmp += Aerow[i] * stages[i,:]
#     end
#     return tmp
# end

@inline function lincomN(Airow, stages, nlfunc::Function; ϵ, C)
    tmp = Airow[1]*nlfunc(stages[1,:], ϵ, C)
    for i = 2 : length(Airow)
        tmp += Airow[i]*nlfunc(stages[i,:], ϵ, C)
    end
    return tmp
end

@inline function ARK3RT(h::Float64, z::Array{T,1}, RHS::Array{T,1}; ω, ϵ, C) where T<:ComplexF64
    invd = 1 ./(1 .- (h*ARK3.d*im*ω))
    #stages
    ks = zeros(ComplexF64,length(ARK3.b),3)
    ks[1,:] = z;
    #L.* zhat + NLfunc(zhat, fP, k); #first stage
    for i = 2 :length(ARK3.b)
        ks[i,:] = invd .* (z + h*(
            im* (ω.* lincom(ARK3.Ae[i-1,1:i-1], ks[1:i-1,:]))
            + lincomN(ARK3.Ai[i-1,1:i-1], ks[1:i-1,:], NLfunc, ϵ=ϵ, C=C)
            )
        )
    end
    return z + h*(im*(ω.*lincom(ARK3.b, ks)) + lincomN(ARK3.b, ks, NLfunc, ϵ=ϵ, C=C))
end

@inline function ARK4RT(h::Float64, z::Array{T,1}, RHS::Array{T,1}; ω, ϵ, C) where T<:ComplexF64
    invd = 1 ./(1 .- (h*ARK4.d*im*ω))
    #stages
    ks = zeros(ComplexF64,length(ARK4.b),3)
    ks[1,:] = z;
    #L.* zhat + NLfunc(zhat, fP, k); #first stage
    for i = 2 :length(ARK4.b)
        ks[i,:] = invd .* (z + h*(
            im* (ω.* lincom(ARK4.Ae[i-1,1:i-1], ks[1:i-1,:]))
            + lincomN(ARK4.Ai[i-1,1:i-1], ks[1:i-1,:], NLfunc, ϵ=ϵ, C=C)
            )
        )
    end
    return z + h*(im*(ω.*lincom(ARK4.b, ks)) + lincomN(ARK4.b, ks, NLfunc, ϵ=ϵ, C=C))
end





@inline function no0rem(x::Int, y::Int)
    r = rem(x,y)
    if r == 0
        return y
    else
        return r
    end
end

@inline function IFab1(state::Array{T,1}, ϵ, C, h, ω,
    tendlist::Array{Array{T,1},1}, i::Int; init::Bool=false) where T<:ComplexF64
    update =  state + h * tendlist[1] 
    ind = 2
    if init==false
        ind = no0rem(i,1)
    end
    tendlist[ind] = NLfunc(state, ϵ, C)
    return update, tendlist
end

@inline function IFab2(state::Array{T,1}, ϵ, C, h, ω,
    tendlist::Array{Array{T,1},1}, i::Int,
    xind::Array{Array{Int,1},1}=[[2,1],[1,2]]; init::Bool=false) where T<:ComplexF64
    update =  state + h * (
        3/2 * tendlist[xind[no0rem(i,2)][1]] 
        -1/2 * tendlist[xind[no0rem(i,2)][2]] 
        )
    ind = 3
    if init==false
        ind = no0rem(i,2)
    end
    tendlist[ind] = NLfunc(state, ϵ, C)
    return update, tendlist
end

@inline function IFab3(state::Array{T,1}, ϵ, C, h, ω,
    tendlist::Array{Array{T,1},1}, i::Int,
    xind::Array{Array{Int,1},1}=[[3,2,1],[1,3,2],[2,1,3]]; init::Bool=false) where T<:ComplexF64
    update =  state + h * (
        23/12 * tendlist[xind[no0rem(i,3)][1]] 
        -16/12 * tendlist[xind[no0rem(i,3)][2]] 
        +5/12 * tendlist[xind[no0rem(i,3)][3]]
        )
    ind = 4;
    if init ==false
        ind = no0rem(i,3)
    end
    tendlist[ind] = NLfunc(state, ϵ, C)
    return update, tendlist
end

@inline function IFab4(state::Array{T,1}, ϵ, C, h, ω,
    tendlist::Array{Array{T,1},1}, i::Int,
    xind::Array{Array{Int,1},1}=[[4,3,2,1],[1,4,3,2],[2,1,4,3],[3,2,1,4]]) where T<:ComplexF64
    #i>5 is the true index number.
    #This function will calculate the state at t_i. 
    #state is the state at t_{i-1}
    update =  state + h * (
        55/24 * tendlist[xind[no0rem(i,4)][1]] 
        -59/24 * tendlist[xind[no0rem(i,4)][2]] 
        +37/24 * tendlist[xind[no0rem(i,4)][3]] 
        -9/24 * tendlist[xind[no0rem(i,4)][4]] 
        )
    tendlist[no0rem(i,4)] = NLfunc(state, ϵ, C)
    return update, tendlist
end


function IFab_amp(N::Int, h::Float64, every::Int, IC::Array{T,1}; 
    ω, ϵ, C, ab::Int, name::String="zAmp",
    msfunc::Array{Function,1}=[IFab1, IFab2, IFab3, IFab4]) where T<:ComplexF64
    stepper = msfunc[ab];
    tendlist = Array{Array{T,1},1}(undef, ab)

    state = deepcopy(IC)
    newtxt!(abs.(state), name=name)
    tendlist[1] = NLfunc(state, ϵ, C)
    for i = 1 : ab-1
        state, tendlist = msfunc[i](state, ϵ, C, h, ω, tendlist, i, init=true)
        if rem(i, every) ==1
            addtxt!(abs.(exp.(im*h*i*ω) .* state), name=name)
        end
    end
    for i = ab: N
        state, tendlist = stepper(state, ϵ, C, h, ω, tendlist, i)
        if rem(i, every)==1
            addtxt!(abs.(exp.(im*h*i*ω) .* state), name=name)
        end
    end
end


function IFab(N::Int, h::Float64, every::Int, IC::Array{T,1}; 
    ω, ϵ, C, ab::Int, name::String="zAmp",
    msfunc::Array{Function,1}=[IFab1, IFab2, IFab3, IFab4]) where T<:ComplexF64
    stepper = msfunc[ab];
    tendlist = Array{Array{T,1},1}(undef, ab)

    state = deepcopy(IC)
    newtxt!(state, name=name)
    tendlist[1] = NLfunc(state, ϵ, C)
    for i = 1 : ab-1
        state, tendlist = msfunc[i](state, ϵ, C, h, ω, tendlist, i, init=true)
        if rem(i, every) ==1
            addtxt!(exp.(im*h*i*ω) .* state, name=name)
        end
    end
    for i = ab: N
        state, tendlist = stepper(state, ϵ, C, h, ω, tendlist, i)
        if rem(i, every)==1
            addtxt!(exp.(im*h*i*ω) .* state, name=name)
        end
    end
end

function RT_amp(N::Int, h::Float64, every::Int, IC::Array{ComplexF64,1}; 
    ω, ϵ, C, stepper::Function=RK4RT, name::String="zAmp")
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
    ω, ϵ, C, stepper::Function=RK4RT, name::String="z")
    tend = deepcopy(z)
    newtxt!(z, name=name)
    for i=2:N+1
        z = stepper(h, z, tend, ω=ω, ϵ=ϵ, C=C)
        #print(z)
        if rem(i, every) ==1 || every ==1
            addtxt!(z, name=name)
        end
    end
end

function plot_RT(name::String="z")
    return 0
end
