clear all
%data = load('../../data/IFRK3h05T1500.mat');
data = load('./I305D32N09.mat');
[Psi, T, S] = boxify3NL(data.xf);
%%
run dd2D
Psibox =real(ifft2(boxify(Psi, dp.Nx, dp.Nz)));
Tbox =real(ifft2(boxify(T, dp.Nx, dp.Nz)));
Sbox =real(ifft2(boxify(S, dp.Nx, dp.Nz)));
%%
run cmap
figure
subplot(221)
imagesc(dp.x, dp.z, Psibox)
title('Psi')
colormap(OrPu/255)
colorbar
subplot(222)
imagesc(dp.x,dp.z,Tbox)
title('T')
colorbar
subplot(223)
imagesc(dp.x,dp.z,Sbox)
colormap(OrPu/255)
title('S')
colorbar
subplot(224)
t = linspace(0,1500,length(data.ES));
plot(t,data.ES, 'DisplayName','ES')
hold on;
plot(t, data.FS, 'DisplayName','FS')
legend()
title("Energies")


%%
t = linspace(0,1500,length(data.ES));
figure
subplot(121)
plot(t, data.ES)
title('ES')
subplot(122)
plot(t, data.FS)
title('FS')