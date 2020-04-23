%Time-stepping 
T = 100;
h = 5e-2;
M = 500%T/h;
every = 100;
bname = join(['IFRK3',sprintf('%02d',h*1e2)],'');
data = load(join(['../../data/',bname,'.mat'],''));

%% Initial Condition
xIC = randn(data.dp.NxNz,3)/(9*data.dp.NxNz^2);

xf = IFRK(xIC, M, h, every, data.dp, data.RK);