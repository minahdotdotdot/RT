%% Problem defininig params (including discretization)
N=2^7; 
D=32;
Pr=1;
h=0.05;

%% Define all params 
a_ratio=2;
[dp,pp]=define_params(N,D,Pr,a_ratio);

%% Set up naming strings.
bname=sprintf('I3%03d',h*1e3);
dname=sprintf('D%2dN%02dP%1da%1d',D,log2(dp.Nx),Pr,a_ratio);
pname=sprintf('../../data/%s%sparams.mat',bname,dname);

%% Get or construct linear operator.
worker=4;
L=getL(N,D,Pr,dp,pp,[pname,worker]);

%% Get other operators. We will be using IFRK with RK3. 
bs = 3*2^7; % Size of one matrix dim assigned to each worker.
disp('Now generating exponentiated matrices.')
RK = genIFRKexpL("RK3", h, L, worker, bs, pp, dp);
save(pname,'RK', 'dp','-append');

%% Set IC.
xIC = randn(dp.NxNz,3)/sqrt(dp.NxNz); % Real in physical space.
xIC = dd2Dfft(xIC, dp);

%% Set simulation length.
dT = 2000;
M = dT/h;
every = 100; %Severy = 5; (To save pictures)
Severy=-1;

%% Run simulation.
name='I7'; %
disp('Start running simulation.')
[xf, ES, FS] = IFRK(xIC, M, h, every, Severy, name, dp, RK);
finaldata=sprintf('../../data/%s%s.mat',bname,dname);
save(finaldata,'xf', 'ES', 'FS','dT');

%% Plot fields and energy evolution.
highres=true;
plotfinal(finaldata, dp, bname, dname, highres);
