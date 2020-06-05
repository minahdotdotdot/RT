%% setup IFRK
run dd2D.m;%_IFsetup.m
name = 'A7'
%% Time-stepping 
dT = 100;
%h=5e-3;
dname=join(['D32N',sprintf('%02d',log2(dp.Nx)),'P1'],'');%'D32N07';
bname = join(['I3',sprintf('%03d',h*1e3)],'');
M = dT/h;
every = 1000; Severy = 10; 
data = load(join(['../../data/',bname,dname,'params.mat'],''));
RK = data.RK; dp = data.dpNL;
clear data;
%% Initial Condition
xIC = randn(dp.NxNz,3)/sqrt(dp.NxNz);
xIC = dd2Dfft(xIC, dp);
%data = load('./I305D32N07.mat');
%xIC = data.xf;
%data = load('../../testrealIC.mat');
%xIC = data.xf;
clear data

%% Run Simulation
[xf, ES, FS] = IFRK_saveS(xIC, M, h, every,Severy, name, dp, RK);
%save('../../testreal2IC.mat','xf','ES','FS','dT');
save(join(['../../data/',bname,dname,'.mat'],''),'xf', 'ES', 'FS','dT');

% Plot
%data = load('../../testrealIC.mat');
data = load(join(['../../data/',bname,dname,'.mat'],''));
x=dp.x/dp.l_o;
z=dp.z/dp.l_o;dT=data.dT;
ES = data.ES(1:end-1);
FS = data.FS(1:end-1);
t = linspace(0,dT,length(ES));
%t = linspace(0,dT,length(data.ES)-1);
[Psi, T, S] = boxify3NL(data.xf);
Psibox =real(ifft2(boxify(Psi, dp.Nx, dp.Nz)));
Tbox =real(ifft2(boxify(T, dp.Nx, dp.Nz)));
Sbox =real(ifft2(boxify(S, dp.Nx, dp.Nz)));

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
plot(t, ES, 'DisplayName','ES')
%plot(t,data.ES(1:end-1), 'DisplayName','ES')
hold on;
plot(t, FS,'DisplayName','FS')
%plot(t, data.FS(1:end-1), 'DisplayName','FS')
legend()
title("Energies")
saveas(gcf,join(['../../plots/',bname,dname,'T',sprintf('%05d',dT),'.png'],''))
%savefig(join(['../../plots/',bname,dname,'T',sprintf('%05d',T),'.png'],''))
