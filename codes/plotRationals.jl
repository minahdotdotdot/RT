include("readwrite.jl")
include("setupMMT.jl")
using PyPlot, LaTeXStrings, MAT

#=
file = matopen("../data/R025.mat");
R025 = read(file, "R025");
close(file)
file = matopen("../data/R250.mat");
R250= read(file, "R250")
close(file)

trueexp1 = exp.(0.025*L);
trueexp2 = exp.(0.25*L);
fig, ax = subplots();
#ax.set_xscale("log")
ax.set_yscale("log")
for deg = 4:12
	err1 = norm(trueexp1 - R025[:,deg-3],2);
	err2 = norm(trueexp2 - R250[:,deg-3],2);
	scatter(deg, err1, c=:blue)
	scatter(deg, err2, c=:red)
	#loglog(abs.(real.(err)), abs.(imag.(err)), label="deg="*string(deg))
end
scatter([], [], c=:blue, label="h=0.025")
scatter([], [], c=:red, label="h=0.250")
legend()
title("Errors for rational approximations of exp(hL)")
xlabel("Degree of rational approximation")
ylabel("2-norm Error")

=#

#=name = "bigICRK3"
solhat = readCfile(name)
sol = ifft(solhat, 2)

T = size(sol)[1]
for i = 1 : 5: T
	fig, ax = subplots()
	plot(real.(sol[i,:]), label="real")
	plot(imag.(sol[i,:]), label="imag")
	ylim([min(minimum(real.(sol)), minimum(imag.(sol))), 
		max(maximum(real.(sol)), maximum(imag.(sol)))])
	legend()
	savefig("../plots/"*name*string(i)*".png")
	close(fig)
end=#

#=f200 = readdlm("../txtfiles/f200.txt"); plot(f200, label="0.2");
f100 = readdlm("../txtfiles/f100.txt"); plot(f100, label="0.1");
f050 = readdlm("../txtfiles/f050.txt"); plot(f050, label="0.05");
f025 = readdlm("../txtfiles/f025.txt"); plot(f025, label="0.025");
f010 = readdlm("../txtfiles/f010.txt"); plot(f010, label="0.010")
legend()
xlabel("Time/10")
title("Varying Forcing")
=#

H=[0.04; 0.025; 0.01]
colors =[:blue, :red, :green]
#trueexp1 = exp.(h*L);
#trueexp2 = exp.(0.5*h*L);
fig, ax = subplots();
ax.tick_params(labelsize="large")

#ax.set_xscale("log")
ax.set_yscale("log")

file=matopen("../data/R2.mat")
R2=read(file,"R2");
close(file)
hL = [];
for i = 1 : length(H)
	global hL = vcat(hL, H[i]*L, 0.5*H[i]*L)
end
hL = vcat(hL, 0)
trueexp = exp.(hL);

for deg = 4 : 12
	err = norm(trueexp - R2[:,deg-3],2)/norm(trueexp,2);
	scatter(deg, err, c=:black, s=100, marker="x")
end
scatter([], [], c=:black, s=100, marker="x", label="D")
R2=0

file=matopen("../data/R.mat")
R=read(file,"R");
close(file)

for i = 1 : length(H)
	h = H[i];
	hL = vcat(h*L, 0.5*h*L, 0);
	trueexp = exp.(hL);
	for deg = 4:12
		err = maximum(abs.(trueexp - R[i,:,deg-3]))/maximum(abs.(trueexp))
		#err = norm(trueexp - R[i,:,deg-3],2)/norm(trueexp,2);
		scatter(deg, err, c=colors[i])
		#scatter(deg, err2, c=:red)
		#loglog(abs.(real.(err)), abs.(imag.(err)), label="deg="*string(deg))
	end
	scatter([], [], c=colors[i], label="h="*string(h))
end

#scatter([], [], c=:blue, label="h="*string(h))
#scatter([], [], c=:red, label="h="*string(.5*h))
legend(fontsize=13)
title("Errors for rational approximations on "*L"D_h"*"'s and D", fontsize=14)
xlabel("Degree of rational approximation", fontsize=14)
#ylabel("2-norm Error")
ylabel("Max-norm Error", fontsize=14)
ylim(1e-15, 1e-4)


#=
fig,ax = subplots();
ax.tick_params(labelsize="large")
for i = 1 : length(H)
	h = H[i]
	hL = vcat(h*L, 0.5*h*L, 0);
	loglog(abs.(real.(hL)), abs.(imag.(hL)),c=colors[i], label="h="*string(h))
end
legend(fontsize=13)
title(L"D_h"*": h-dependent Domains", fontsize=14)
xlabel("Absolute real values", fontsize=14)
ylabel("Absolute imaginary values", fontsize=14)
=#
