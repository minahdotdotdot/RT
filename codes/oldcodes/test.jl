include("RT.jl")
using PyPlot
seed=9
Random.seed!(seed);

ω = [-1, 3, -2]
C = floatRT(5)
IC = onUnitCircle(3);
#RK4_step, EUimex_step
hs = collect(-4:0)

#### "true" solve by RK4 h=1e-5.
name="true_seed"*string(seed); h=exp10(-6)
N = Int(ceil(T/h)); 
every = Int(ceil(N/L)) #only save 1001 values total
RT_amp(N, h, every, IC, ω=ω, ϵ=ϵ, C=C, stepper=RK4_step, name=name);


#### Methods we are comparing.
steppers = [IFE_step, ETD1_step, CNimex_step]
ss = length(steppers)
for i in hs, j = 1 : ss
	name=string(steppers[j])*string(i)
	h = exp10(i); N = Int(ceil(T/h)); 
	if N > L
		every = Int(ceil(N/L))
	else
		every = 1
	end
	RT_amp(N, h, every, IC, ω=ω, ϵ=ϵ, C=C, stepper=steppers[j], name=name);		
end



###PLOT errors for each of the 3 waves separately. 
colors = ["r","b","g"]
for k = 1 : 3
	for i in hs
		h = exp10(i); 
		tsol = readdlm("true_seed"*string(seed)*".txt")
		subplot(100+length(hs)*10+i+1-hs[1])
		for j = 1 : ss
			name=string(steppers[j])*string(i)
			y = readdlm(name*".txt")[:,k]
			if length(y)==length(x)
				semilogy(x, abs.(tsol[:,k]-y), c=colors[j], label=name)
			elseif length(y) < length(x)
				nn = Int((length(x)-1)/(length(y)-1))
				semilogy(x[1:nn:end], abs.(tsol[1:nn:end,k]-y), c=colors[j], label=name)
			else
				"True solution cannot have lower resolution."
			end
		end
		legend()
	end
	savefig("seed"*string(seed)*"-errz"*string(k)*".png")
	close()
end

####PLOT all methods for each of the 3 waves separately.
steppers = [IFE_step, ETD1_step, CNimex_step]
ss = length(steppers)
for k = 1 : 3
	for i in hs
		tsol = readdlm("true_seed"*string(seed)*".txt")
		subplot(100+length(hs)*10+i+1-hs[1])
		for j = 1 : ss-1
			name=string(steppers[j])*string(i)
			y = readdlm(name*".txt")[:,k]
			if length(y)==length(x)
				plot(x,y, c=colors[j], label=name)
			elseif length(y) < length(x)
				nn = Int((length(x)-1)/(length(y)-1))
				plot(x[1:nn:end], y, c=colors[j], label=name)
			else
				"True solution cannot have lower resolution."
			end
		end
		plot(x,tsol[:,k], c="k", label="true")
		legend()
	end
	savefig("seed"*string(seed)*"-z"*string(k)*".png")
	close()
end

####PLOT each method separately.
for i in hs
	for j = 1 : ss
		name=string(steppers[j])*string(i)
		y = readdlm(name*".txt")
			if size(y)[1]==length(x)
				plot(x,y, c=colors[j], label=name)
			elseif size(y)[1] < length(x)
				nn = Int((length(x)-1)/(size(y)[1]-1))
				plot(x[1:nn:end], y, c=colors[j], label=name)
			else
				"True solution cannot have lower resolution."
			end
		legend()
		savefig("everything"*string(i)*string(steppers[j])*".png")
		close()
	end
end

plot(x, readdlm("true_seed"*string(seed)*".txt"), c="k", label="true")
legend()
savefig("everything_seed"*string(seed)*".png")
close()

