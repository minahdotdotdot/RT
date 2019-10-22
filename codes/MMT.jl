include("readwrite.jl")
using FFTW, LinearAlgebra

struct funcparams
    α:: Float64
    β:: Float64
    λ:: Complex{Float64}
    funcparams(α, β, λ) = new(copy(α), copy(β), copy(λ))
end

@inline function Lfunc(zhat::Array{ComplexF64,1}, fp::funcparams)
    return abs.(k).^fp.α .*zhat
end

@inline function NLfunc(zhat::Array{ComplexF64,1}, fp::funcparams, k, FD)
    if fp.β != 0
        zr = ifft(abs.(k) .^(β/4) .* zhat)
        return -im* fp.λ * abs.(k) .^(β/4) .* fft(abs.(zr).^2 .* zr) + D .* zhat + F
    else
        zr = ifft(zhat)
        tmp = -im* fp.λ*fft(abs.(zr).^2 .*zr)
        tmp[(N/4+2):(3*N/4)] .= 0
        return  tmp + FD .* zhat
    end
end

struct eRKTableau
    A # matrix
    b # stage weights
    c # t increments
    eRKTableau(A, b, c) = new(copy(A),copy(b), copy(c))
end

function IFRK_step(z::Array{ComplexF64,1}, h::Float64, 
    L, NLfunc::Function, nlfP::funcparams, 
    RKT::eRKTableau, ks, k, FD)
    #k = zeros(eltype(z), length(RKT.b), length(z))
    ks[1,:] = NLfunc(z, nlfP, k, FD)
    for i = 2 :length(RKT.b)
        PP=h*Transpose(ks[1:i-1,:])*RKT.A[i,1:i-1]
        ks[i,:] = exp.(-h*RKT.c[i]*L) .* NLfunc(exp.(h*RKT.c[i]*L) .*(z + PP), nlfP, k, FD)
    end
    vb = Transpose(ks)*b
    #return k, vb, L * (z+ (h*vb))
    return exp.(h*L) .* (z+ (h*vb))
end

function IFRK!(M::Int, every::Int, IC::Array{ComplexF64,1}, h::Float64, 
    L, NLfunc::Function, fP::funcparams, RKT::eRKTableau, k; name::String)
    # FFT into Fourier space
    zhat = fft(IC)*0.001
    newtxt!(zhat, name=name)
    ks = zeros(eltype(zhat), length(RKT.b), length(IC)) #RK stages (allocate memory)
    #forcing term
    N = length(zhat)
    FD = zeros(length(k)); 
    FD[[6+1, 7+1, 8+1, 9+1, -6+(N+1), -7+(N+1), -8+(N+1), -9+(N+1)]] .= 0.2
    #add damping term
    FD[2:end] += -196.61 * (abs.(k[2:end]).^(-8)) - 53.9* (abs.(0.001 * k[2:end]) .^ 16); 
    FD[1]= -200.0
    for t = 1 : M
        zhat = IFRK_step(zhat, h, L, NLfunc, fP, RKT, ks, k, FD)
        if rem(t,every)==1
            if isnan(abs(zhat))==true || isinf(abs(zhat))==true
                break
            end
            addtxt!(zhat, name=name)
        end
    end
end

# Problem Parameters
λ = 1   #Defocusing MMT model
α = 1/2
β = 0
fP = funcparams(α, β, λ)
