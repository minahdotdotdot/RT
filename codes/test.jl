include("RT.jl")
using PyPlot
seed=12
Random.seed!(seed);
#T = 1000; 
#ϵ = 0.001;
h = exp10(-3); 

ω = intRT(3)
C = floatRT(5)
IC = onUnitCircle(3);
#RK4_step, EUimex_step

for i in R
	name="true"*string(i)
	T = exp10(i); ϵ=exp10(-i); h=exp10(-8)
	N = Int(ceil(T/h)); 
	every = Int(ceil(N/100)) #only save 101 values total
	RT_amp(N, h, every, IC, ω=ω, ϵ=ϵ, C=C, stepper=RK4_step, name=name);	
end


steppers = [IFE_step, ETD1_step, CNimex_step]
colors = ["r","b","g"]
ss = length(steppers)
R = collect(2:3) 
for i in R
	for j = 1 : ss
		name=string(steppers[j])*string(i)
		T = exp10(i); ϵ=exp10(-i); 
		N = Int(exp10(i+3)); #Int(ceil(T/h)); 
		every =Int(exp10(i+1)); #Int(ceil(N/100))
		RT_amp(N, h, every, IC, ω=ω, ϵ=ϵ, C=C, stepper=steppers[j], name=name);		
	end
end

#Plot amplitudes for each of the variables separately. 

for k = 1 : 3
	for i in R
		tsol = readdlm("true"*string(i)*".txt")
		print(100+length(R)*10+i-1, "\n")
		subplot(100+length(R)*10+i-1)
		for j = 1 : ss
			name=string(steppers[j])*string(i)
			semilogy(abs.(tsol[:,k]-readdlm(name*".txt")[:,k]),c=colors[j], label=name)
		end
		legend()
	end
	savefig("seed"*string(seed)*"-errz"*string(k)*".png")
	close()
end

steppers = [IFE_step, ETD1_step, CNimex_step,:true]
colors = ["r","b","g","k"]
ss = length(steppers)
for k = 1 : 3
	for i in R
		tsol = readdlm("true"*string(i)*".txt")
		print(100+length(R)*10+i-1, "\n")
		subplot(100+length(R)*10+i-1)
		for j = 1 : ss
			name=string(steppers[j])*string(i)
			plot(readdlm(name*".txt")[:,k],c=colors[j], label=name)
		end
		legend()
	end
	savefig("seed"*string(seed)*"-z"*string(k)*".png")
	close()
end
