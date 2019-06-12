include("RT.jl")
using PyPlot
seed=13
Random.seed!(seed);
h = exp10(-2); #fixed 
ω = intRT(3)
C = floatRT(5)
IC = onUnitCircle(3);
#RK4_step, EUimex_step
R = collect(2:3); #10 .^(R) gives final time, 10 .^(-R) gives ϵ.

for i in R
	name="true"*string(i)
	T = 5*exp10(i); ϵ=exp10(-i); h=exp10(-5)
	N = Int(ceil(T/h)); 
	every = Int(ceil(N/100)) #only save 101 values total
	RT_amp(N, h, every, IC, ω=ω, ϵ=ϵ, C=C, stepper=RK4_step, name=name);	
end

steppers = [IFE_step, ETD1_step, CNimex_step]
colors = ["r","b","g"]
ss = length(steppers)
for i in R, j = 1 : ss
	name=string(steppers[j])*string(i)
	T = 5*exp10(i); ϵ=exp10(-i); 
	N = Int(ceil(T/h)); 
	every = Int(ceil(N/100))
	RT_amp(N, h, every, IC, ω=ω, ϵ=ϵ, C=C, stepper=steppers[j], name=name);		
end

#Plot amplitudes for each of the variables separately. 

for k = 1 : 3
	for i in R
		T = 5*exp10(i); x = range(0, T, length=101)
		tsol = readdlm("true"*string(i)*".txt")
		subplot(100+length(R)*10+i-1)
		for j = 1 : ss
			name=string(steppers[j])*string(i)
			semilogy(x, abs.(tsol[:,k]-readdlm(name*".txt")[:,k]),
				c=colors[j], label=name)
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
		T = 5*exp10(i); x = range(0, T, length=101)
		tsol = readdlm("true"*string(i)*".txt")
		print(100+length(R)*10+i-1, "\n")
		subplot(100+length(R)*10+i-1)
		for j = 1 : ss
			name=string(steppers[j])*string(i)
			plot(x, readdlm(name*".txt")[:,k],c=colors[j], label=name)
		end
		legend()
	end
	savefig("seed"*string(seed)*"-z"*string(k)*".png")
	close()
end

#plot 3 waves for each method separately.
for i in R
	T = exp10(i); x = range(0, T, length=101)
	for j = 1 : ss
		name=string(steppers[j])*string(i)
		plot(x, readdlm(name*".txt"),c=colors[j], label=name)
		legend()
		savefig("everything"*string(i)*string(steppers[j])*".png")
		close()
	end
end
