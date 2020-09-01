using DifferentialEquations, PyPlot, DelimitedFiles
include("RT.jl")

struct RTparams
	ω::Vector{Float64}
    ϵ::Float64
    C::Vector{Float64}
    RTparams(ω,ϵ,C) = new(copy(ω), copy(ϵ), copy(C))
end

function RTtend!(z, p, t)
	dz = deepcopy(z);
	dz[1] = p.C[1]*(z[2]*z[3]-z[5]*z[6]);
	dz[4] = p.C[1]*(-z[2]*z[6]-z[3]*z[5]);
	dz[2] = p.C[2]*(z[1]*z[3]-z[4]*z[6]);
	dz[5] = p.C[2]*(-z[1]*z[6]-z[3]*z[4]);
	dz[3] = p.C[3]*(z[1]*z[2]-z[4]*z[5]);
	dz[6] = p.C[3]*(-z[1]*z[5]-z[2]*z[4]);
	#dz = p.ϵ*dz + im*p.ω.*z
	#dz = p.ϵ*dz + vcat(-p.ω.*z[4:6],p.ω.*z[1:3]);
	return p.ϵ*dz + vcat(-p.ω.*z[4:6],p.ω.*z[1:3]);
end

function reshapetomat(u::Vector{Vector{Float64}})
	newu = zeros(ComplexF64,length(u),3)
	newmag = zeros(Float64, length(u),3)
	for i = 1 : length(u)
		newu[i,1] = u[i][1]+im*u[i][4];
		newu[i,2] = u[i][2]+im*u[i][5];
		newu[i,3] = u[i][3]+im*u[i][6];
		newmag[i,1] = abs(newu[i,1]);
		newmag[i,2] = abs(newu[i,2]);
		newmag[i,3] = abs(newu[i,3]);
	end
	return newu, newmag
end

seed=1205
Random.seed!(seed);
ϵ = 0.01;             # Nonlinear scale
ω = [-1, 3, -2];      # "slow" wave numbers
C = floatRT(5);       # Energy conserving constants
IC = onUnitCircle(3)  # Initial condition
ICr = vcat(real.(IC), imag.(IC));
p = RTparams(ω,ϵ,C);
tsol = readdlm("../txtfiles/true_seed"*string(seed)*".txt");

T = 1000.;
tspan = (0., T);
dz = deepcopy(ICr); #p=3.;
prob = ODEProblem(RTtend!, ICr, tspan, p);
#prob = ODEProblem(easytend2!, ICr, tspan, p);



# Choose algorithm.
alg = ImplicitEuler();
#alg = ImplicitMidpoint();
#alg = RadauIIA3();
#alg = RadauIIA5();
#alg = PDIRK44();

# Set time step-size.
dt = 0.001;#5e-5

strideI = Int(1/dt);
ind = 1:strideI:strideI*1000+1;
#sol=solve(prob,alg,dt=dt);
#truesol = 5*exp.(p*sol.t);
#scatter(1:length(sol.t)-1, sol.t[2:end]-sol.t[1:end-1] .- 0.05)
#complexsol, absvalsol = reshapetomat(sol.u);
sol = Vector{Vector{Float64}}(undef,Int(1000/dt)+1);
sol[1] = ICr;
for i = 2 : length(sol)
	tspan = (0, dt);
	prob = ODEProblem(RTtend!, sol[i-1], tspan, p);
	sol[i] = solve(prob, alg, dt=dt).u[2];
end
complexsol, absvalsol = reshapetomat(sol);
writedlm("../txtfiles/ImplicitEulerh="*string(dt)*"seed"*string(seed)*".txt", absvalsol[ind,:])
#=
fig,ax=subplots()
subplot(121)
#semilogy(1:1000, abs.(tsol[2:end,1]-absvalsol[ind,1]), label="z1")
#semilogy(1:1000, abs.(tsol[2:end,2]-absvalsol[ind,2]), label="z2")
semilogy(1:1000, abs.(tsol[2:end,3]-absvalsol[ind,3]), label="z3")
#plot(sol.t[1:strideI:end], absvalsol[1:strideI:end,1],label="z1")
#plot(sol.t[1:strideI:end], absvalsol[1:strideI:end,2],label="z2")
#plot(sol.t[1:strideI:end], absvalsol[1:strideI:end,3],label="z3")
title(@sprintf("h=%.3f",dt))
legend()
=#
dt = 0.05;#5e-5
strideI = Int(1/dt);
ind = 1:strideI:strideI*1000+1;
#sol=solve(prob,alg,dt=dt);
#truesol = 5*exp.(p*sol.t);
#scatter(1:length(sol.t)-1, sol.t[2:end]-sol.t[1:end-1] .- 0.05)
#complexsol, absvalsol = reshapetomat(sol.u);
sol = Vector{Vector{Float64}}(undef,Int(1000/dt)+1);
sol[1] = ICr;
for i = 2 : length(sol)
	tspan = (0, dt);
	prob = ODEProblem(RTtend!, sol[i-1], tspan, p);
	sol[i] = solve(prob, alg, dt=dt).u[2];
end
complexsol, absvalsol = reshapetomat(sol);
writedlm("../txtfiles/ImplicitEulerh="*string(dt)*"seed"*string(seed)*".txt", absvalsol[ind,:])
#=
subplot(122)
#semilogy(1:1000, abs.(tsol[2:end,1]-absvalsol[ind,1]), label="z1")
#semilogy(1:1000, abs.(tsol[2:end,2]-absvalsol[ind,2]), label="z2")
semilogy(1:1000, abs.(tsol[2:end,3]-absvalsol[ind,3]), label="z3")
#plot(sol.t[1:strideI:end], absvalsol[1:strideI:end,1],label="z1")
#plot(sol.t[1:strideI:end], absvalsol[1:strideI:end,2],label="z2")
#plot(sol.t[1:strideI:end], absvalsol[1:strideI:end,3],label="z3")
title(@sprintf("h=%.3f",dt))
legend()
suptitle("Implicit Euler")
=#

alg = PDIRK44();
# Set time step-size.
dt = 1.
strideI = Int(1/dt);
ind = 1:strideI:strideI*1000+1;
sol=solve(prob,alg,dt=dt);
#truesol = 5*exp.(p*sol.t);
#scatter(1:length(sol.t)-1, sol.t[2:end]-sol.t[1:end-1] .- 0.05)
complexsol, absvalsol = reshapetomat(sol.u);
writedlm("../txtfiles/PDIRK44h="*string(dt)*"seed"*string(seed)*".txt", absvalsol[ind,:])
#=fig,ax=subplots()
subplot(121)
#semilogy(sol.t[ind], abs.(tsol[2:end,1]-absvalsol[ind,1]), label="z1")writedlm("../txtfiles/ImplicitEulerh="*string(dt)*"seed"*string(seed)*".txt"
writedlm("../txtfiles/ImplicitEulerh="*string(dt)*"seed"*string(seed)*".txt"
writedlm("../txtfiles/PDIRK44h="*string(dt)*"seed"*string(seed)*".txt"
writedlm("../txtfiles/PDIRK44h="*string(dt)*"seed"*string(seed)*".txt"
#semilogy(sol.t[ind], abs.(tsol[2:end,2]-absvalsol[ind,2]), label="z2")
semilogy(sol.t[ind], abs.(tsol[2:end,3]-absvalsol[ind,3]), label="z3")
ylim(1e-1,1.5e0)
#plot(sol.t[1:strideI:end], absvalsol[1:strideI:end,1],label="z1")
#plot(sol.t[1:strideI:end], absvalsol[1:strideI:end,2],label="z2")
#plot(sol.t[1:strideI:end], absvalsol[1:strideI:end,3],label="z3")
title(@sprintf("h=%.2f",dt))
#legend()
=#
dt = 4.;#5e-5
sol=solve(prob,alg,dt=dt);
#truesol = 5*exp.(p*sol.t);
#scatter(1:length(sol.t)-1, sol.t[2:end]-sol.t[1:end-1] .- 0.05)
complexsol, absvalsol = reshapetomat(sol.u);
writedlm("../txtfiles/PDIRK44h="*string(dt)*"seed"*string(seed)*".txt", absvalsol)
#=
dt = Int(dt);
subplot(122)
ylim(1e-1,1.5e0)
#semilogy(sol.t[2:end], abs.(tsol[dt:dt:end,1]-absvalsol[2:end,1]), label="z1")
#semilogy(sol.t[2:end], abs.(tsol[dt:dt:end,2]-absvalsol[2:end,2]), label="z2")
semilogy(sol.t[2:end], abs.(tsol[dt:dt:end,3]-absvalsol[2:end,3]), label="z3")
#plot(sol.t[1:strideI:end], absvalsol[1:strideI:end,1],label="z1")
#plot(sol.t[1:strideI:end], absvalsol[1:strideI:end,2],label="z2")
#plot(sol.t[1:strideI:end], absvalsol[1:strideI:end,3],label="z3")
title(@sprintf("h=%.2f",dt))
#legend()
suptitle("PDIRK44: 2 processor 4th order diagonally non-adaptive implicit method")
=#