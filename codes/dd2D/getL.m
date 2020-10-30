function L = getL(N,D,Pr,dp,pp,dname,p)
	filename=sprintf('../../data/D=%sL.mat',dname);
	if(exist(filename,'file'))
		disp('L file already exists.');
		data=load(filename);
		L=data.L;
	else
		disp('L file does not exist, generating it now.');
		L = genL(pp,dp,p(2));
		pname = '../../data/'+p(1)+p(2)+'params.mat';
		save(pname,'L');
	end
end