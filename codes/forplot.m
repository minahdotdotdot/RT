file=load('../data/L.mat');
L=file.L;
H = [0.04 0.025 0.01];
hL=[];
for i=1:length(H)
	h=H(i);
	hL = [hL; h*L; .5*h*L];
end
hL = [hL; 0];
R2 = zeros(length(hL),9);
Z=exp(hL);
for deg = 4 : 12
	[r, pol, res, zer, zj, fj, wj] = aaa(Z,hL,'degree', deg);
	R2(:,deg-3)=r(hL);
end
save('../data/R2.mat','R2')
    
%{
H = [0.04 0.025 0.01];
R = zeros(length(H),length(L)*2+1,9);
for i=1:length(H)
	h=H(i);
	hL = [h*L; .5*h*L; 0];
	Z=exp(hL);
	for deg = 4 : 12
	    [r, pol, res, zer, zj, fj, wj] = aaa(Z,hL,'degree', deg);
	    R(i,:,deg-3)=r(hL);
	end
end
save('../data/R.mat','R')
%}