%Time-stepping 
T = 1500;
h = 5e-2;
M = T/h;
every = 1000;
bname = join(['IFRK3',sprintf('%02d',h*1e2)],'');
data = load(join(['../../data/',bname,'.mat'],''));

%% Initial Condition
xIC = randn(data.dp.NxNz,3)/(data.dp.NxNz);


[xf, ES, FS] = IFRK(xIC, M, h, every, data.dp, data.RK);
save('IFRK3h05.mat','xf', 'ES', 'FS');
