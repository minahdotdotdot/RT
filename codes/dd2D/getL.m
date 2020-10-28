function L = getL(N,D,Pr,dp,pp,p)
	filename=sprintf('../../data/D%2dN%02dP%1dL.mat',D,log2(N),Pr);
	if(exist(filename,'file'))
		disp('L file already exists.');
		data=load(filename);
		L=data.L;
	else
		disp('L file does not exist, generating it now.');
		L = genL(pp,dp,p(3));
		pname = '../../data/'+p(1)+p(2)+'params.mat';
		save(pname,'L');
	end
end