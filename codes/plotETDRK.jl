using PyPlot, LaTeXStrings, ExponentialUtilities, LinearAlgebra

function plotETDerror!(h::Vector{T}, c::Vector{Vector{T}}, b::Vector{Vector{Vector{S}}}, ω::S,
 s::Vector{Int}, colors::Vector{Symbol}, labels=Vector{String};plotfunc::Function=plot) where S<:Union{T,Complex{T}} where T<:AbstractFloat
	#fig, ax = subplots()
	#ax.set_yscale("log")
	for i = length(s) : -1 : 1
		err = zeros(length(h))
		for j = 1 : length(b[i])
			err = err + b[i][j].*exp.((ω*(c[i][j]-1))*h)
		end
		err = 1 .- err#abs.(err)
		if i == 3
			plotfunc(h, err, c=colors[i], label=labels[i], linestyle="--")
		elseif i == 4
			plotfunc(h, err, c=colors[i], label=labels[i], linewidth=5)
		else
			plotfunc(h, err, c=colors[i], label=labels[i])
		end
	end
	#plotfunc(ones(length(h)), zeros(length(h)), c=:black, label="truth")
	xlabel(L"h\omega", fontsize=13)
	#ylabel(L"Arg\left(1-\sum_{j=1}^sb_j(ih\omega)e^(c_j-1)ih\omega}\right)")
	ylabel(L"1-\sum_{j=1}^sb_j(ih\omega)e^{i(c_j-1)h\omega}",fontsize=13)
	legend()
	#ylim(5e-15,2)
	#xt = [0, pi/2, pi, 3/2*pi, 2*pi]
	#xlabels = ("0", L"\frac{\pi}{2}", L"\pi", L"\frac{3\pi}{2}", L"2\pi")
	#xticks(xt, xlabels)
	title("Interaction cofficient errors for ETDRK schemes, "*L"\omega=1")
end

#h = Vector(range(0, stop=2*pi, length=1001));
h = Vector(exp10.(range(-16, stop=log10(2*pi), length=1001)))
ω = 1.0*im
n = 4
b = Vector{Vector{Vector{typeof(ω)}}}(undef, n)
c = Vector{Vector{Float64}}(undef, n)
s = Vector{Int}(undef, n)
colors = Vector{Symbol}(undef, n)
labels= Vector{String}(undef, n)
L =Diagonal(h*ω);
phis = diag.(phi(L, 4)[2:4]); #phi functions 1 to 3.

i = 1
b[i]=[phis[1]]
c[i]= [0.0]
s[i] = length(b);
colors[i]=:red
labels[i]="ETD1"

i = 2
b[i]=[phis[1]-phis[2], phis[2]]
c[i]= [0.0; 1]
s[i] = length(b);
colors[i]=:green
labels[i]="ETD2RK"

i = 3
b[i]=[phis[1]-3*phis[2]+4*phis[3], 
4*phis[2]-8*phis[3], 
-phis[2]+4*phis[3]]
c[i]= [0.0; 1/2; 1]
s[i] = length(b);
colors[i]=:blue
labels[i]="ETD3RK"

i = 4
b[i]= [phis[1]-3*phis[2]+4*phis[3], 
2*phis[2]-4*phis[3], 
2*phis[2]-4*phis[3], 
-phis[2]+4*phis[3]]
c[i]=[0; 1/2; 1/2; 1]
s[i] = length(b);
colors[i]=:yellow
labels[i]="ETD4RK"

#=
i = 5
b[i]=[1/6; 1/3; 1/3; 1/6]
c[i]= [0.0; 1/2; 1/2; 1]
s[i] = length(b);
colors[i]=:blue
labels[i]="Classic RK4"

i = 6
b[i]=[1/8; 3/8; 3/8; 1/8]
c[i]= [0.0; 1/3; 2/3; 1]
s[i] = length(b);
colors[i]=:purple
labels[i]="3/8 Rule RK4"=#

plotETDerror!(h, c, b, ω, s, colors, labels, plotfunc=plot)



