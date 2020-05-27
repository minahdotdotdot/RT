run dd2D % to get h
%Time-stepping 
%h = 5e-2;
dname=join(['D32N',sprintf('%02d',log2(dp.Nx)),'P1'],'');%dname='D32N07P1'
bname = join(['A3',sprintf('%03d',h*1e3)],'');
%% Set-up Linear operator
%par = 0;tic
%L=genL(pp,dp,par);toc
workers = 4;
par = workers;tic
delete(gcp('nocreate'))
parpool(workers)
L = genL(pp, dp, par);toc
pname = join(['../../data/',bname,dname,'params.mat'],'');
save(pname,'L');

%% Set-up Additive-Runge_kutta scheme.
run ARKtableaus
arks = ARK3;
clear ARK4
save(pname,'arks','-append');

%% Set-up IMEX operator: invd = (I - h*arks.d*L)^{-1}
hd = h*arks.d;
invd = genImexInv(h, hd, pp, dp, par);
save(pname,'invd','-append');

run setupNL
save(pname,'dpNL','-append');



