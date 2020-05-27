clear all
%data = load('../../data/IFRK3h05T1500.mat');
data = load('../../data/I350D32N09P1.mat');
[Psi, T, S] = boxify3NL(data.xf);
%%
run dd2D
Psibox =real(ifft2(boxify(Psi, dp.Nx, dp.Nz)));
Tbox =real(ifft2(boxify(T, dp.Nx, dp.Nz)));
Sbox =real(ifft2(boxify(S, dp.Nx, dp.Nz)));
%%
x=dp.x/dp.l_o;
z=dp.z/dp.l_o;
t = linspace(0,4800,length(data.ES)-1);

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


%%
figure
subplot(121)
plot(t, data.ES)
title('ES')
subplot(122)
plot(t, data.FS)
title('FS')