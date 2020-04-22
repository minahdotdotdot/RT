 function x = function ETDRK(x, M, h, L, Nx, Nz, km, kk, mm, every, name)
 	A=0;b=0;c=0;
 	if name == "RK2"
 		A=fill(zeros(2), 1,2)
		A[1,1]=[1,0];A=Array{Any,2}(A)
		c=[1];c=Array{Any,1}(c)
		b=[[1.0,-1],[0,1]];b=Array{Any,1}(b)
 	elseif name == "RK3"
 		A=fill(zeros(3),2,3);
		A[1,1]=[1/2,0,0];
		A[2,1:2]=[[-1,0,0],[2,0,0]];A=Array{Any,2}(A)
		c=[1/2,1];c=Array{Any,1}(c)
		b=[[1.0,-3,4],[0,4,-8],[0,-1,4]];b=Array{Any,1}(b)
 	elseif name == "RK4"
 		A=fill(zeros(4),3,4);
		A[1,1]=[1/2,0,0,0];
		A[2,2]=[1/2,0,0,0];
		A[3,[1,3]]=[[1/2,0,0,0],[1,0,0,0]];A=Array{Any,2}(A)
		c=[1/2,1/2,1];c=Array{Any,1}(c)
		b=[[1.0,-3,4,0],[0,2,-4,0],[0,2,-4,0],[0,-1,4,0]];b=Array{Any,1}(b);
 	end
 end

 function [A,b,c] = buildABC(A, b, c, h, L)
 end