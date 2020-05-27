run dd2D
%Time-stepping 
h = 5e-3; 
dname=join(['D32N',sprintf('%02d',log2(dp.Nx)),'P1'],'');%dname='D32N07P1'
bname = join(['I3',sprintf('%03d',h*1e3)],'');

%parpool %parpool only called in fillc
%linear operator
%workers = 8; %3*Nx*Nz = 3(2^{2p+log_2(aspectratio)}) 
%i.e. aspect ratio = 2, p = 7 --> 3NxNz = 3(2^{15})
%bs = 3*2^5; % block size: log_2(15 - 3 = 12)
%N=2^7, aspectratio = 2 -> bs = 3*2^7 || 32  sec
%                       -> bs = 3*2^6 || 28  sec
%                       -> bs = 3*2^5 || 5.3 sec
%                       -> bs = 3*2^4 || 6.2 sec

workers = 4;
%i.e. aspect ratio = 1, p = 10 --> 3NxNz = 3(2^{20})
bs = 3*2^7; % Each worker does: TOTAL: 3*2^{20 - 4 = 16} since workers=16=2^4
%                             AT ITER: bs=3*2^7, 2^{16-7} ITERs
%p=log_2(N) :  3NxNz   -> genL  || expL
%  7        : 3*2^{14} ->   25s ||   1.4s
%  8        : 3*2^{16} ->   21s ||   3.5s
%  9        : 3*2^{18} ->   22s ||  20. s
% 10        : 3*2^{20} -> 1637s ||    . s

%par = 0;tic
%L=genL(pp,dp,par);toc
par = workers;tic
delete(gcp('nocreate'))
parpool(workers)
L = genL(pp, dp, par);toc
save(join(['../../data/',bname,dname,'params.mat'],''),'L');
RK = genIFRKexpL("RK3", h, L, workers, bs, pp, dp);
run setupNL
save(join(['../../data/',bname,dname,'params.mat'],''),'RK', 'dpNL','-append');
