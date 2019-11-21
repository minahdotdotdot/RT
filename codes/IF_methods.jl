include("MMT.jl")

function fillc(c, h, L)
	tmp = Array{typeof(L),1}(undef,0)
	for i = 1 : length(c)
		push!(tmp, exp.(-h*c[i]*L))
		push!(tmp, exp.(h*c[i]*L))
	end
	push!(tmp,exp.(h*L))
	return tmp
end

#RK2

#RK3
A = hcat([0; .5; -1],[0; 0; 2]);
b = [1/6; 2/3; 1/6];
c = [1/2; 1]; 
c=fillc(c,h,L)
RK3 = eRKTableau(A, b, c);

#RK4
A = hcat([0; .5; 0; 0], [0; 0; .5; 0], [0; 0; 0; 1.0]); #3by4
b = [1/6; 1/3; 1/3; 1/6];
c = [1/2; 1/2; 1]; 
c=fillc(c,h,L)
RK4 = eRKTableau(A, b, c); 
