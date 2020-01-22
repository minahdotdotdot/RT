include("MMT.jl")
include("al_approx.jl")


function fillc(c, h, L)
	tmp = Array{typeof(L),1}(undef,0)
	for i = 1 : length(c)
		push!(tmp, exp.(-h*c[i]*L))
		push!(tmp, exp.(h*c[i]*L))
	end
	push!(tmp,exp.(h*L))
	return tmp
end

function fillc_rat(c, h, L, deg::Int)
	name = "../data/MMT_deg"*string(deg)
	tmp = Array{typeof(L),1}(undef,0)
	for i = 1 : length(c)
		push!(tmp, evalrat(-h*c[i]*L, name))
		push!(tmp, evalrat(h*c[i]*L), name))
	end
	push!(tmp, evalrat(h*L, name))
	return tmp
end

#RK2

#RK3
A = hcat([0; .5; -1],[0; 0; 2]);
b = [1/6; 2/3; 1/6];
x = [1/2; 1]; 
c=fillc(x,h,L);
crat = fillc_rat(x, h, L, 6);
RK3 = eRKTableau(A, b, c);
RK3_rat = eRKTableau(A, b, crat)

#RK4
A = hcat([0; .5; 0; 0], [0; 0; .5; 0], [0; 0; 0; 1.0]); #3by4
b = [1/6; 1/3; 1/3; 1/6];
x = [1/2; 1/2; 1]; 
c =fillc(x,h,L);
rat = fillc_rat(x, h, L, 6);
RK4 = eRKTableau(A, b, c); 
RK4_rat = eRKTableau(A, b, crat); 
