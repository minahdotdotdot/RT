include("MMT.jl")

#RK2


#RK3
A = hcat([0; .5; -1],[0; 0; 2]);
b = [1/6; 2/3; 1/6];
c = [0; 1/2; 1];
RK3 = eRKTableau(A, b, c);

#RK4
A = hcat([0; .5; 0; 0], [0; 0; .5; 0], [0; 0; 0; 1.0]); #3by4
b = [1/6; 1/3; 1/3; 1/6];
c = [0; 1/2; 1/2; 1];
RK4 = eRKTableau(A, b, c);

