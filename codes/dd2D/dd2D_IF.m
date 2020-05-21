%Time-stepping 
T = 1200;
h = 5e-1;
M = T/h;
every = 100;
bname = join(['IFRK3',sprintf('%02d',h*1e2),'D32N07'],'');
data = load(join(['../../data/',bname,'.mat'],''));
RK = data.RK; dp = data.dpNL;
clear data;
%% Initial Condition
xIC = randn(dp.NxNz,3)/sqrt(dp.NxNz);
%data = load('./I305D32N07.mat');
%xIC = data.xf;
%clear data


[xf, ES, FS] = IFRK(xIC, M, h, every, dp, RK);
save('I350D32N07.mat','xf', 'ES', 'FS');
