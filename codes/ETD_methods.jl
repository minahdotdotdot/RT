include("MMT.jl")
#This file must be ran after h and L have already been defined.
#=
Initially, I fill the entries of A matrix and b vector 
with coefficients of the phi functions indexed at 1.
Then, the function buildAbc, which is defined in MMT.jl, builds 
the true operators so that you can apply the RK table directly. 
For example, A[1,1] = [1,0] indicates:
A_{11} = 1 x ϕ₁(c₁hL) + 0 x ϕ₂(c₁hL), 
and b[2]=[0,4,-8] indicates: 
b_{2} = 0 x ϕ₁(hL) + 4 x ϕ₂(hL) - 8 x ϕ₃(hL). 
The i^{th} row of A gets cᵢ multiplied to the argument, hL, 
and the elements of the b vector always have the standard argument, hL.
Finally, the c vector will store {exp(c[i]hL)}i=1^{s} for an s-stage method.
=#

#ETDRK2
A=fill(zeros(2), 1,2)
A[1,1]=[1,0];A=Array{Any,2}(A)
c=[1];c=Array{Any,1}(c)
b=[[1.0,-1],[0,1]];b=Array{Any,1}(b)
A,b,c=buildAbc(A,b,c,h,L)
ETDRK2=ETDRKTableau(A,b,c)

#ETDRK3
A=fill(zeros(3),2,3);
A[1,1]=[1/2,0,0];
A[2,1:2]=[[-1,0,0],[2,0,0]];A=Array{Any,2}(A)
c=[1/2,1];c=Array{Any,1}(c)
b=[[1.0,-3,4],[0,4,-8],[0,-1,4]];b=Array{Any,1}(b)
A,b,c = buildAbc(A,b,c,h,L);
ETDRK3=ETDRKTableau(A,b,c)

#ETDRK4
A=fill(zeros(4),3,4);
A[1,1]=[1/2,0,0,0];
A[2,2]=[1/2,0,0,0];
A[3,[1,3]]=[[1/2,0,0,0],[1,0,0,0]];A=Array{Any,2}(A)
c=[1/2,1/2,1];c=Array{Any,1}(c)
b=[[1.0,-3,4,0],[0,2,-4,0],[0,2,-4,0],[0,-1,4,0]];b=Array{Any,1}(b);
A,b,c=buildAbc(A,b,c,h,L); 
A[3,[1,3]] = [A[1,1]*(exp(.5*h*Diagonal(L))-I),2*A[1,1]]
ETDRK4=ETDRKTableau(A,b,c)

#ETDRK4-B
A=fill(zeros(4),3,4);
A[1,1]=[1/2,0,0,0];
A[2,1:2]=[[1/2,-1,0,0],[0,1,0,0]];
A[3,[1,3]]=[[1,-2,0,0],[0,2,0,0]];A=Array{Any,2}(A)
c=[1/2,1/2,1];c=Array{Any,1}(c)
b=[[1.0,-3,4,0],[0,2,-4,0],[0,2,-4,0],[0,-1,4,0]];b=Array{Any,1}(b);
A,b,c=buildAbc(A,b,c,h,L);
ETDRK4B=ETDRKTableau(A,b,c)

ETDdict = Dict(
	"ETDRK2"=> ETDRK2, 
	"ETDRK3" => ETDRK3, 
	"ETDRK4" => ETDRK4, 
	"ETDRK4B" => ETDRK4B
	);
