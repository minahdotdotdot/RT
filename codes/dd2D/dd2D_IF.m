%Time-stepping 
T = 1500;
h = 5e-2;
M = T/h;
every = 1000;
bname = join(['IFRK3',sprintf('%02d',h*1e2),'D32N10'],'');
data = load(join(['../../data/',bname,'.mat'],''));
RK = data.RK; dp = data.dp;
clear data;
%% Initial Condition
xIC = randn(data.dp.NxNz,3)/(data.dp.NxNz);


[xf, ES, FS] = IFRK(xIC, M, h, every, dp, RK);
save('I305D32N10.mat','xf', 'ES', 'FS');
