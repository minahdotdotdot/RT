include("RT.jl")
using PyPlot
seed=7
Random.seed!(seed);
ϵ = 0.01; T=1000; L=1000; x=range(0,1000,length=L+1);
ω = [-1, 3, -2]
C = floatRT(5)
IC = onUnitCircle(3);
#RK4, EUimex
hs = [0.02, 0.025, 0.05]
test="_B"
#### "true" solve by RK4 h=1e-6.
trname="true_seed"*string(seed)*"B"; h=exp10(-6)
N = Int(ceil(T/h)); 
every = Int(ceil(N/L)) #only save 1001 values total
#=
RT_amp(N, h, every, IC, ω=ω, ϵ=ϵ, C=C, stepper=RK4, name=trname);
=#

#### Methods we are comparing.
steppers = [IFE, ETD1, CNimex]
ss = length(steppers)
#=
for i = 1 : length(hs), j = 1 : ss
	txtname=string(steppers[j])*test*string(i)
	h = hs[i]; N = Int(ceil(T/h)); 
	if N > L
		every = Int(ceil(N/L))
	else
		every = 1
	end
	RT_amp(N, h, every, IC, ω=ω, ϵ=ϵ, C=C, stepper=steppers[j], name=txtname);		
end
=#

###PLOT errors for each of the 3 waves separately. 
colors = ["r","b","g"]
#=
for k = 1 : 3
	for i =1: length(hs)
		tsol = readdlm("../txtfiles/"*trname*".txt")
		subplot(100+length(hs)*10+i)
		title("h="*string(hs[i]))
		for j = 1 : ss
			txtname="../txtfiles/"*string(steppers[j])*test*string(i)
			y = readdlm(txtname*".txt")[:,k]
			if length(y)==length(x)
				semilogy(x, abs.(tsol[:,k]-y), c=colors[j], label=string(steppers[j]))
				ylim(top=10.0)
			elseif length(y) < length(x)
				nn = Int((length(x)-1)/(length(y)-1))
				semilogy(x[1:nn:end], abs.(tsol[1:nn:end,k]-y), c=colors[j], label=string(steppers[j]))
				ylim(top=10.0)
			else
				"True solution cannot have lower resolution."
			end
		end
		legend()
	end
	savefig("../plots/seed"*string(seed)*"-errz"*string(k)*".png")
	close()
end
=#

####PLOT each wave separately.
for k = 1 : 3
	for i =1:length(hs)
		tsol = readdlm("../txtfiles/"*trname*".txt")
		subplot(100+length(hs)*10+i)
		title("h="*string(hs[i]))
		for j = 1 : ss
			txtname="../txtfiles/"*string(steppers[j])*test*string(i)
			y = readdlm(txtname*".txt")[:,k]
			if length(y)==length(x)
				plot(x,y, c=colors[j], label=string(steppers[j]))
			elseif length(y) < length(x)
				nn = Int((length(x)-1)/(length(y)-1))
				plot(x[1:nn:end], y, c=colors[j], label=string(steppers[j]))
			else
				"True solution cannot have lower resolution."
			end
		end
		plot(x,tsol[:,k], c="k", label="true")
		legend()
	end
	savefig("../plots/seed"*string(seed)*test*"-z"*string(k)*".png")
	close()
end



####PLOT each method separately.
for j = 1 : ss
	for i =1:length(hs)
		subplot(100*length(hs)+10+i)
		title(string(hs[i]))
		txtname="../txtfiles/"*string(steppers[j])*test*string(i)
		y = readdlm(txtname*".txt")
		if size(y)[1]==length(x)
			for k = 1 : 3
				plot(x,y[:,k], label=k)
			end
		elseif size(y)[1] < length(x)
			nn = Int((length(x)-1)/(size(y)[1]-1))
			for k = 1 : 3
				plot(x[1:nn:end], y[:, k], label=k)
			end
		else
			"True solution cannot have lower resolution."
		end
		legend()
	end
	savefig("../plots/"*test*string(steppers[j])*".png")
	close()
end

plot(x, readdlm("../txtfiles/"*trname*".txt"), c="k", label="true")
legend()
savefig("../plots/"*test*"everything_seed"*string(seed)*".png")
close()
