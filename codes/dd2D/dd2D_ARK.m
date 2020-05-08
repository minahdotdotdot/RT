%Time-stepping 
T = 1500;
h = 5e-2;
M = T/h;
every = 100;
bname = join(['ARK3',sprintf('%02d',h*1e2),'D32N09S1'],'');
data = load(join(['../../data/',bname,'.mat'],''));

L = data.L; invd = data.invd; arks = data.arks;
clear data;
%% Initial Condition
run dd2D
xIC = randn(dp.NxNz,3)/(dp.NxNz);


[xf, ES, FS] = ARK(xIC, M, h, every, dp, arks, L, invd);
save('A305D32N09S1.mat','xf', 'ES', 'FS');
