include("../codes/RT.jl")

#Set-up
seed=1205
Random.seed!(seed);
ϵ = 0.01;             # Nonlinear scale
ω = [-1, 3, -2];      # "slow" wave numbers
C = floatRT(5);       # Energy conserving constants
IC = onUnitCircle(3)  # Initial condition
T=1000;               # Final time
L=1000;               # Number of state/amplitudes to save

#Time-steppers
#IF_steppers = [IFE, IFRK2, IFRK3, IFRK4];
#ETD_steppers = [ETD1, ETDRK2, ETDRK3, ETDRK4, ETDRK4B]
#steppers=[IF_steppers;ETD_steppers];
steppers1=[IFE, ETD1, CNimex, EUimex, ImplicitEuler]
steppers2=[IFRK3, ETDRK3, ARK3RT, PDIRK44]
#ss = length(steppers)

#True solution
tsol = readdlm("../txtfiles/true_seed"*string(seed)*".txt")[:,3]

#Other solutions
sol1 = Matrix{Vector{Float64}}(undef, length(steppers1),2)
for i = 1 : length(steppers1)
    sol1[i,1] = readdlm("../txtfiles/"*string(steppers1[i])*"h=0.001seed"*string(seed)*".txt")[:,3]
    sol1[i,2] = readdlm("../txtfiles/"*string(steppers1[i])*"h=0.05seed"*string(seed)*".txt")[:,3]
end
sol2 = Matrix{Vector{Float64}}(undef,length(steppers2),2);
for i = 1 : length(steppers2)
    sol2[i,1] = readdlm("../txtfiles/"*string(steppers2[i])*"h=1.0seed"*string(seed)*".txt")[:,3]
    sol2[i,2] = readdlm("../txtfiles/"*string(steppers2[i])*"h=4.0seed"*string(seed)*".txt")[:,3]
end

#Plot
hs = [0.001, 1.0, 0.05, 4.0];
fig, axs = plt.subplots(nrows=2, ncols=2, sharex=true, sharey=true)
for i = 1:4
	#axs[i].axhline(0.1, linestyle="dashed", c=:grey, label="0.1 relative error")
	sol=0
	if i == 1 
		sol=view(sol1,:,1)
	elseif i == 2
		sol=view(sol2,:,1)
	elseif i ==3
		sol=view(sol1,:,2)
	elseif i ==4
		sol=view(sol2,:,2)
	end
	steppers=0;ind=1:1:T+1;
	if i ==1||i==3
		steppers=steppers1
	else
		steppers=steppers2
		ind = 1:Int(hs[i]):T+1
	end
	for j = 1 : length(steppers)-1
    	axs[i].semilogy(ind, abs.((tsol[ind]-sol[j])./tsol[ind]), label=string(steppers[j]))
    end
    j = length(steppers)
    axs[i].semilogy(ind, abs.((tsol[ind]-sol[j])./tsol[ind]), linestyle="dashed", label=string(steppers[j]))
    axs[i].set_ylim(1e-7, 1.2)
    axs[i].set_title("h="*string(hs[i]))
    #axs[i].set_title("Relative errors of |z₃| with h="*string(hs[i]))
    if i >2
   		axs[i].legend(loc="lower right")
   	end
end
#for i = 1 :2
# 	axs[i].axhline(0.1, linestyle="dashed", c=:grey, label="0.1 relative error")
# 	for j = 1 : ss
#     	axs[i].semilogy(ind,abs.((tsol[ind]-sol2[j,i-2])./tsol[ind]), label=string(steppers[j]))
#     end
#     axs[i].set_ylim(1e-7, 1.)
#     axs[i].set_title("Relative errors of |z₃| with h="*string(hs[i]))
#     axs[i].legend(loc="bottom right")
# end
tight_layout(pad=1.0)