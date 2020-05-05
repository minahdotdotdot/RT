run dd2D
%Time-stepping 
h = 5e-2;
bname = join(['IFRK3',sprintf('%02d',h*1e2)],'');

%parpool
%linear operator
workers = 8; %3*Nx*Nz = 3(2^{2p+log_2(aspectratio)}) 
%i.e. aspect ratio = 2, p = 7 --> 3NxNz = 3(2^{15})
bs = 3*2^5; % block size: log_2(15 - 3 = 12)
%N=2^7, aspectratio = 2 -> bs = 3*2^7 || 32  sec
%                       -> bs = 3*2^6 || 28  sec
%                       -> bs = 3*2^5 || 5.3 sec
%                       -> bs = 3*2^4 || 6.2 sec

L = genL(pp, dp);
RK = genIFRKexpL("RK3", h, L, workers, bs, pp, dp);
run setupNL
save(join(['../../data/',bname,'.mat'],''),'RK', 'dp');
