run dd2D
%Time-stepping 
h = 5e-2;
bname = join(['IFRK3',sprintf('%02d',h*1e2)],'');

%parpool %parpool only called in fillc
%linear operator
%workers = 8; %3*Nx*Nz = 3(2^{2p+log_2(aspectratio)}) 
%i.e. aspect ratio = 2, p = 7 --> 3NxNz = 3(2^{15})
%bs = 3*2^5; % block size: log_2(15 - 3 = 12)
%N=2^7, aspectratio = 2 -> bs = 3*2^7 || 32  sec
%                       -> bs = 3*2^6 || 28  sec
%                       -> bs = 3*2^5 || 5.3 sec
%                       -> bs = 3*2^4 || 6.2 sec

workers = 16;
%i.e. aspect ratio = 1, p = 10 --> 3NxNz = 3(2^{20})
bs = 3*2^7; % blocksize: log_2(20 - 3 = 17)
%bs = 3*2*7;
%N=2^7 : 3NxNz=3*2^{14} --->  4   sec
%N=2^8 : 3NxNz=3*2^{16} ---> 11.8 sec
%N=2^9 : 3NxNz=3*2^{18} --->



L = genL(pp, dp);
RK = genIFRKexpL("RK3", h, L, workers, bs, pp, dp);
run setupNL
save(join(['../../data/',bname,'.mat'],''),'RK', 'dpNL');
