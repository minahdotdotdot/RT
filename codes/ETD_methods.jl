include("MMT.jl")

#ETDRK2
A=fill(zeros(2), 1,2)
A[1,1]=[1,0];A=Array{Any,2}(A)
c=[1]
b=[[1.0,-1],[0,1]];b=Array{Any,1}(b)
A,b=buildAandb(A,b,c,h,L)
ETDRK2=ETDRKTableau(A,b,c)

#ETDRK3
A=fill(zeros(3),2,3);
A[1,1]=[1/2,0,0];
A[2,1:2]=[[-1,0,0],[2,0,0]];A=Array{Any,2}(A)
c=[1/2,1];
b=[[1.0,-3,4],[0,4,-8],[0,-1,4]];b=Array{Any,1}(b)
A,b = buildAandb(A,b,c,h,L);
ETDRK3=ETDRKTableau(A,b,c)

#ETDRK4
A=fill(zeros(4),3,4);
A[1,1]=[1/2,0,0,0];
A[2,2]=[1/2,0,0,0];
A[3,[1,3]]=[[1/2,0,0,0],[1,0,0,0]];A=Array{Any,2}(A)
c=[1/2,1/2,1];
b=[[1.0,-3,4,0],[0,2,-4,0],[0,2,-4,0],[0,-1,4,0]];b=Array{Any,1}(b);
A,b=buildAandb(A,b,c,h,L); A[3,1] *= (exp(h/2*Diagonal(L)) -I);
ETDRK4=ETDRKTableau(A,b,c)

#ETDRK4-B
A=fill(zeros(4),3,4);
A[1,1]=[1/2,0,0,0];
A[2,1:2]=[[1/2,-1,0,0],[0,1,0,0]];
A[3,[1,3]]=[[1,-2,0,0],[0,2,0,0]];A=Array{Any,2}(A)
c=[1/2,1/2,1];
b=[[1.0,-3,4,0],[0,2,-4,0],[0,2,-4,0],[0,-1,4,0]];b=Array{Any,1}(b);
A,b=buildAandb(A,b,c,h,L);
ETDRK4B=ETDRKTableau(A,b,c)