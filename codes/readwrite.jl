using DelimitedFiles
@inline function newtxt!(zAmp::Array{T,1}; name::String="zAmp") where T<:Float64
    writedlm("../txtfiles/"*name*".txt", transpose(zAmp))
end

@inline function addtxt!(zAmp::Array{T,1}; name::String="zAmp") where T<:Float64
    open("../txtfiles/"*name*".txt", "a") do io
        writedlm(io, transpose(zAmp))
    end
end

@inline function newtxt!(zAmp::T; name::String="zAmp") where T<:Float64
    writedlm("../txtfiles/"*name*".txt", zAmp)
end

@inline function addtxt!(zAmp::T; name::String="zAmp") where T<:Float64
    open("../txtfiles/"*name*".txt", "a") do io
        writedlm(io, transpose(zAmp))
    end
end

@inline function newtxt!(zAmp::Array{T,1}; name::String="zAmp") where T<:Complex{Float64}
    writedlm("../txtfiles/"*name*"_Re.txt", transpose(real.(zAmp)))
    writedlm("../txtfiles/"*name*"_Im.txt", transpose(imag.(zAmp)))
end

@inline function addtxt!(zAmp::Array{T,1}; name::String="zAmp") where T<:Complex{Float64}
    open("../txtfiles/"*name*"_Re.txt", "a") do io
        writedlm(io, transpose(real.(zAmp)))
    end
    open("../txtfiles/"*name*"_Im.txt", "a") do io
        writedlm(io, transpose(imag.(zAmp)))
    end
end

@inline function readCfile(name::String)
    return readdlm("../txtfiles/"*name*"_Re.txt") + im * readdlm("../txtfiles/"*name*"_Im.txt")
end


using FFTW, PyPlot
