%% setup IFRK
run dd2D.m;%_IFsetup.m

%% Time-stepping 
T = 10000;h=5e-3;
dname=join(['D32N',sprintf('%02d',log2(dp.Nx)),'P1'],'');%'D32N07';
bname = join(['I3',sprintf('%03d',h*1e3)],'');
M = T/h;
every = 100; 
data = load(join(['../../data/',bname,dname,'params.mat'],''));
RK = data.RK; dp = data.dpNL;
clear data;
%% Initial Condition
xIC = randn(dp.NxNz,3)/sqrt(dp.NxNz);
%data = load('./I305D32N07.mat');
%xIC = data.xf;
%clear data

[xf, ES, FS] = IFRK(xIC, M, h, every, dp, RK);
save(join(['../../data/',bname,dname,'.mat'],''),'xf', 'ES', 'FS');

%% Plot
x=dp.x/dp.l_o;
z=dp.z/dp.l_o;
t = linspace(0,T,length(data.ES)-1);

run cmap
figure
subplot(221)
imagesc(x, z, Psibox)
title('Psi')
colormap(OrPu/255);
colorbar
subplot(222)
imagesc(x,z,Tbox)
title('T')
colorbar
subplot(223)
imagesc(x,z,Sbox)
colormap(OrPu/255)
title('S')
colorbar
subplot(224)
plot(t,data.ES(1:end-1), 'DisplayName','ES')
hold on;
plot(t, data.FS(1:end-1), 'DisplayName','FS')
legend()
title("Energies")
savefig(join(['../../plots/',bname,dname,'T',sprintf('%05d',T),'.png'],''))
