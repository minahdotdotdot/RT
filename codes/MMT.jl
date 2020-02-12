include("readwrite.jl")
using FFTW, LinearAlgebra, Printf

####### MMT Problem Set-up #######
struct funcparams
    α:: Float64
    β:: Float64
    λ:: Complex{Float64}
    F:: Float64
    D:: Array{Float64,1}
    funcparams(α, β, λ, F, D) = new(copy(α), copy(β), copy(λ), copy(F), copy(D))
end

@inline function NLfunc(zhat::Array{ComplexF64,1}, fP::funcparams, k)
    if fP.β != 0
        zr = ifft(abs.(k) .^(fP.β/4) .* zhat)
        return -im* fP.λ * abs.(k) .^(fP.β/4) .* fft(abs.(zr).^2 .* zr)
    else
        zr = ifft(zhat);
        tmp = -im* fP.λ*fft(abs.(zr).^2 .*zr) #tmp[Int(N/4+2):Int(3*N/4)] .= 0
        return  tmp
    end
end

@inline function NLfunc(zhat::Array{ComplexF64,2}, fP::funcparams, k)
    tmp = copy(zhat)
    if fP.β != 0
        zr = ifft(abs.(k) .^(fP.β/4) .* zhat, 2)
        return -im* fP.λ * abs.(k) .^(fP.β/4) .* fft(abs.(zr).^2 .* zr, 2)
    else
        zr = ifft(zhat, 2);
        tmp = -im* fP.λ*fft(abs.(zr).^2 .*zr, 2) #tmp[Int(N/4+2):Int(3*N/4)] .= 0
        return  tmp
    end
end


####### IFRK Set-up #######
struct eRKTableau
    A # matrix        :: (s-1)-by-s matrix
    b # stage weights :: s-length vector
    c # for i > j, e^((c_i - c_j)hL)
    x # 
    eRKTableau(A, b, c, x) = new(copy(A),copy(b), copy(c), copy(x))
end
@inline function lincomIF(A, ks)
    tmp = A[1] .* ks[1,:]
    for i = 2 : length(A)
        tmp += A[i] .* ks[i,:]
    end
    return tmp
end

function IFRK_step(zhat::Array{ComplexF64,1}, h::Float64, 
    L, NLfunc::Function, fP::funcparams, RKT::eRKTableau, ks, k)
    ks[1,:] = c[1,1] .* NLfunc(zhat, fP, k); #c[1,1] should be just ones.
    for i = 2 :length(RKT.b)
        PP = h* lincomIF(RKT.A[i,1:i-1] .* RKT.c[i,1:i-1], ks[1:i-1,:]); # sum
        ks[i,:] = NLfunc(PP+RKT.c[i,i].*zhat, fP, k)
        #PP = h* lincom(RKT.A[i,1:i-1], ks[1:i-1,:])
        #ks[i,:] = RKT.c[2*(i-2)+1] .* NLfunc(RKT.c[2*(i-2)+2].*(zhat+PP), fP, k)
    end
    PP = h* lincomIF(RKT.b .* RKT.c[end,:], ks)
    return PP + RKT.c[end-1,end].*zhat
    #return RKT.c[end] .* (zhat + h*lincom(RKT.b, ks))
end

function IFRK!(M::Int, every::Int, IC::Array{ComplexF64,1}, h::Float64, 
    L, NLfunc::Function, fP::funcparams, RKT::eRKTableau, k; name::String, cont::Bool=false)
    # FFT into Fourier space
    zhat = fft(IC)*0.001; N = length(zhat); #zhat[Int(N/4+2):Int(3*N/4)].= 0.0;
    if cont==false
        newtxt!(zhat, name=name); kmax = sqrt(maximum(abs.(k)));
    end
    ks = zeros(eltype(zhat), length(RKT.b), length(IC)); #RK stages (allocate memory)
    #forcing term
    N = length(zhat);
    for t = 1 : M
        zhat = IFRK_step(zhat, h, L, NLfunc, fP, RKT, ks, k)[:,1]
        if rem(t,every)==1 || every==1
            if any(isnan,zhat)  || any(isinf,zhat)
                error("Blowup!!! at ND time="*string(t*h))
            end
            addtxt!(zhat, name=name)
        end
    end
end


####### ETDRK Problem Set-up #######
using ExponentialUtilities
struct ETDRKTableau
    A #matrix :: (s-1)-by-s matrix
    b #stages :: s-length vector
    c #times  :: s-length vector
    ETDRKTableau(A, b, c) = new(copy(A),copy(b), copy(c))
end

function buildAbc(A, b, c, h, L)
    s = size(A)[2]
    L = Diagonal(L)
    for i = 1 : size(A)[1]
        phis = phi(c[i]*h*L, s)[2:end]
        for j = 1 : s
            A[i,j] = sum(A[i,j].*phis)
        end
    end
    phis = phi(h*L,s)[2:end]
    for j = 1 : length(b)
        b[j] = sum(b[j].*phis)
    end
    for j = 1 : length(c)
        c[j] = exp(c[j]*h*L)
    end
    return Array{typeof(A[1,1]),2}(A), Array{typeof(b[1,1]),1}(b), Array{typeof(c[1]),1}(c)
end

@inline function lincom(A, ks)
    tmp = A[1] * ks[1,:]
    for i = 2 : length(A)
        tmp += A[i] * ks[i,:]
    end
    return tmp
end

function ETDRK_step(zhat::Array{ComplexF64,1}, h::Float64, 
    L, NLfunc::Function, fP::funcparams, RKT::ETDRKTableau, ks, k, exphL)
    ks[1,:] = NLfunc(zhat, fP, k); #first stage
    for i = 2 :length(RKT.b)
        ks[i,:] = NLfunc(c[i-1]*zhat+h*lincom(RKT.A[i-1,1:i-1], ks[1:i-1,:]), fP,k)
    end
    return (exphL.*zhat) + h*lincom(RKT.b, ks)
end

function ETDRK!(M::Int, every::Int, IC::Array{ComplexF64,1}, h::Float64, 
    L, NLfunc::Function, fP::funcparams, RKT::ETDRKTableau, k; name::String, cont::Bool=false)
    # FFT into Fourier space
    zhat = fft(IC)*0.001; N = length(zhat); #zhat[Int(N/4+2):Int(3*N/4)].= 0.0;
    if cont==false
        newtxt!(zhat, name=name); kmax = sqrt(maximum(abs.(k)));
    end
    ks = zeros(eltype(zhat), length(RKT.b), length(IC)); #RK stages (allocate memory)
    #forcing term
    N = length(zhat);
    #fname=name*"f"*string(Int(fP.F*1000), pad=3)
    #newtxt!(maximum(abs.(ifft(zhat)))^2/kmax, name=fname)
    exphL = exp.(h*L)
    for t = 1 : M
        zhat = ETDRK_step(zhat, h, L, NLfunc, fP, RKT, ks, k, exphL)
        if rem(t,every)==1 || every==1
            if any(isnan,zhat)  || any(isinf,zhat)
                error("Blowup!!! at ND time="*string(t*h))#break
            end
            addtxt!(zhat, name=name)
            #addtxt!(maximum(abs.(ifft(zhat)))^2/kmax, name=fname)
        end
    end
end

####### IMEX Problem Set-up #######
struct IMEXTableau
    Ae #:: s-by-s matrix (explicit)
    Ai #:: s-by-s matrix (implicit)
    b #stages :: s-length vector
    c #times  :: s-length vector
    d #diagonal entry of A^(implicit)
    IMEXTableau(Ae, Ai, b, c, d) = new(copy(Ae), copy(Ai), copy(b), copy(c), d)
end
@inline function lincomN(A, ks, fp, k, nlfunc::Function=NLfunc)
    tmp = A[1]*nlfunc(ks[1,:], fp, k)
    for i = 2 : length(A)
        tmp += A[i] * nlfunc(ks[i,:], fp, k)
    end
    return tmp
end

function IMEXRK_step(zhat::Array{ComplexF64,1}, h::Float64, 
    L, NLfunc::Function, fP::funcparams, RKT::IMEXTableau, ks, k, invd)
    ks[1,:] = copy(zhat)
    #L.* zhat + NLfunc(zhat, fP, k); #first stage
    for i = 2 :length(RKT.b)
        ks[i,:] = invd .* (zhat + h*(
            L .* lincom(RKT.Ae[i-1,1:i-1], ks[1:i-1,:])
            + lincomN(RKT.Ai[i-1,1:i-1], ks[1:i-1,:], fP, k, NLfunc)
            )
        )
    end
    return zhat + h*( L .* lincom(RKT.b, ks) + lincomN(RKT.b, ks, fP, k, NLfunc))
end

function IMEXRK!(M::Int, every::Int, IC::Array{ComplexF64,1}, h::Float64, 
    L, NLfunc::Function, fP::funcparams, RKT::IMEXTableau, k; name::String, cont::Bool=false)
    # FFT into Fourier space
    zhat = fft(IC)*0.001; N = length(zhat); #zhat[Int(N/4+2):Int(3*N/4)].= 0.0;
    if cont==false
        newtxt!(zhat, name=name); kmax = sqrt(maximum(abs.(k)));
    end
    ks = zeros(eltype(zhat), length(RKT.b), length(IC)); #RK stages (allocate memory)
    #forcing term
    N = length(zhat);
    invd = 1 ./(1 .- h*RKT.d*L)
    for t = 1 : M
        zhat = IMEXRK_step(zhat, h, L, NLfunc, fP, RKT, ks, k, invd)
        if rem(t,every)==1 || every==1
            if any(isnan,zhat)  || any(isinf,zhat)
                error("Blowup!!! at ND time="*string(t*h))#break
            end
            addtxt!(zhat, name=name)
        end
    end
end
