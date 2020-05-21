run dd2D
%Time-stepping 
h = 5e-2;
bname = join(['ARK3',sprintf('%02d',h*1e2),'D32N07P1'],'');

%% Set-up Linear operator
%par = 0;tic
%L=genL(pp,dp,par);toc
workers = 4;
par = workers;tic
delete(gcp('nocreate'))
parpool(workers)
L = genL(pp, dp, par);toc
save(join(['../../data/',bname,'.mat'],''),'L');

%% Set-up Additive-Runge_kutta scheme.
run ARKtableaus
arks = ARK3
clear ARK4
save(join(['../../data/',bname,'.mat'],''),'arks','-append');

%% Set-up IMEX operator: invd = (I - h*arks.d*L)^{-1}
hd = h*arks.d;
invd = genImexInv(h, hd, pp, dp, par);
save(join(['../../data/',bname,'.mat'],''),'invd','-append');

run setupNL
save(join(['../../data/',bname,'.mat'],''),'dpNL','-append');



