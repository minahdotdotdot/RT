using DelimitedFiles
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