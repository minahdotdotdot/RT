%% setup IFRK
run dd2D_IFsetup.m

%% Time-stepping 
T = 4800;
%h = 5e-2;dname='D32N07';bname = join(['IFRK3',sprintf('%02d',h*1e2),dname],'');
M = T/h;
every = 100; 
data = load(join(['../../data/',bname,'.mat'],''));
RK = data.RK; dp = data.dpNL;
clear data;
%% Initial Condition
xIC = randn(dp.NxNz,3);
%data = load('./I305D32N07.mat');
%xIC = data.xf;
%clear data

[xf, ES, FS] = IFRK(xIC, M, h, every, dp, RK);
save(join(['../../data/I3',sprintf('%02d',h*1e2),dname,'.mat'],''),'xf', 'ES', 'FS');
