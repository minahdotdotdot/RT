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

function intRT()
	ω = rand(1:20,3)
	x = rand(3)
	for i = 1 : 3
    	if x[i]>0.5
    		ω[i] = -ω[i]
    	end
    end
    x = rand();
    if x <= 1/3;
    	ω[1] = - (ω[2]+ω[3])
    elseif x > 1/3 && x <=2/3
    	ω[2] = - (ω[1] + ω[3])
    else
    	ω[3] = - (ω[1] + ω[2])
    end
    return ω
end 
#tendency for a resonant triad with ω = [ω₁, ω₂, ω₃] frequencies and ϵ=slow time scale << 1.
@inline function tendRT(z::Array{T,1}, tend::Array{T,1}; ω, ϵ, C) where T<:ComplexF64
	tend[1] = im*ω[1]*z[1] + ϵ*C[1]*conj(z[2])*conj(z[3])
	tend[2] = im*ω[2]*z[2] + ϵ*C[2]*conj(z[3])*conj(z[1])
	tend[3] = im*ω[3]*z[3] + ϵ*C[3]*conj(z[1])*conj(z[2])
	return tend
end

@inline function RK4_step(h::Float64, z::Array{T,1}, tend::Array{T,1}; ω, ϵ, C) where T<:ComplexF64
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

@inline function newtxt!(zAmp::Array{T,1}) where T<:Float64
	writedlm("./zAmp.txt", zAmp')
end

@inline function addtxt!(zAmp::Array{T,1}) where T<:Float64
    open("./zAmp.txt", "a") do io
    	writedlm(io, zAmp')
    end
end

function RT(N::Int, h::Float64, every::Int, z::Array{ComplexF64,1}; 
	ω, ϵ, C, stepper::Function=RK4_step)
	zAmp = Array{Float64,2}(undef,N+1,3)
	tend = deepcopy(z)
	newtxt!(abs.(z))
	for i=1:N
		z = stepper(h, z, tend,  ω=ω, ϵ=ϵ, C=C)
		if rem(i, every) ==1
			@printf("%d\t",div(i, every))
			addtxt!(abs.(z))
		end
	end
end


Random.seed!(4)
ϵ = 0.01;
ω = [1, 2, -3]
C = [-2, 1, 1]
tend = Array{ComplexF64,1}(undef, 3);

T = 1000; h = 1e-1; every = Int(1e1);
z = onUnitCircle(3);
RT(Int(ceil(T/h)), h, every, z, ω=ω, ϵ=ϵ, C=C);
plot(readdlm("./zAmp.txt"))
savefig("h1.png")
