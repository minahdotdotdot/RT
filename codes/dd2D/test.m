%% get data
data = load('../../data/I3050D32N07P1params.mat');
L = data.L; save('../../data/D32N07P1L.mat','L');
dp7=data.dpNL; clear data
data7 = load('../../data/I3050D32N07P1.mat');

data = load('../../data/I3050D32N08P1params.mat');
L = data.L; save('../../data/D32N08P1L.mat','L');
dp8=data.dpNL; clear data
data8 = load('../../data/I3050D32N08P1.mat');

data = load('../../data/I3050D32N09P1params.mat');
L = data.L; save('../../data/D32N09P1L.mat','L');
clear L
clear data
%%
[Psi7, T7, S7] = boxify3NL(data7.xf);
Psibox7 =real(ifft2(boxify(Psi7, dp7.Nx, dp7.Nz)));
Tbox7 =real(ifft2(boxify(T7, dp7.Nx, dp7.Nz)));
Sbox7 =real(ifft2(boxify(S7, dp7.Nx, dp7.Nz)));

[Psi8, T8, S8] = boxify3NL(data8.xf);
Psibox8 =real(ifft2(boxify(Psi8, dp8.Nx, dp8.Nz)));
Tbox8 =real(ifft2(boxify(T8, dp8.Nx, dp8.Nz)));
Sbox8 =real(ifft2(boxify(S8, dp8.Nx, dp8.Nz)));

run cmap
figure
subplot(221)
imagesc(dp7.x/dp7.l_o, dp7.z/dp7.l_o, Sbox7)
title(join([sprintf('%d',dp7.Nx),' grid at time ', sprintf('%d',data7.dT)],''))
colormap(OrPu/255);
colorbar
subplot(222)
imagesc(dp8.x/dp8.l_o, dp8.z/dp8.l_o, Sbox8)
title(join([sprintf('%d',dp8.Nx),' grid at time ', sprintf('%d',data8.dT)],''))
colorbar

[Psi7, T7, S7] = boxify3NL(data7.xf2);
Psibox7 =real(ifft2(boxify(Psi7, dp7.Nx, dp7.Nz)));
Tbox7 =real(ifft2(boxify(T7, dp7.Nx, dp7.Nz)));
Sbox7 =real(ifft2(boxify(S7, dp7.Nx, dp7.Nz)));

[Psi8, T8, S8] = boxify3NL(data8.xf2);
Psibox8 =real(ifft2(boxify(Psi8, dp8.Nx, dp8.Nz)))
Tbox8 =real(ifft2(boxify(T8, dp8.Nx, dp8.Nz)));
Sbox8 =real(ifft2(boxify(S8, dp8.Nx, dp8.Nz)));

subplot(223)
imagesc(dp7.x/dp7.l_o, dp7.z/dp7.l_o, Sbox7)
title(join([sprintf('%d',dp7.Nx),' grid at time ', sprintf('%d',data7.dT+data7.dT2)],''))
colormap(OrPu/255);
colorbar
subplot(224)
imagesc(dp8.x/dp8.l_o, dp8.z/dp8.l_o, Sbox8)
title(join([sprintf('%d',dp8.Nx),' grid at time ', sprintf('%d',data8.dT+data8.dT2)],''))
colorbar
%%
ES7 = [data7.ES(1:end-1);data7.ES2(1:end-1);data7.ES3(1:end-1)];
ES8 = [data8.ES(1:end-1);data8.ES2(1:end-1)];
semilogy(linspace(0,data7.dT+data7.dT2+data7.dT3,length(ES7)), ES7)
hold on;
semilogy(linspace(0,data8.dT+data8.dT2,length(ES8)), ES8)
legend('N=2^7','N=2^8')
title('Spatially Averaged Salt var energies at different discretizations')
%%
display(.5*sum(Sbox7.^2,'all')/(dp7.NxNz))
display(.5*sum(S7.*conj(S7))/(dp7.NxNz^2))
display(.5*sum(S7.^2)/(dp7.NxNz^2))
%display(data7.ES(end-2:end))
%display(.5*sum(Sbox8.^2,'all')/(dp8.NxNz))
%display(.5*sum(abs(S8).^2)/(dp8.NxNz^2))
%display(data8.ES(end-2:end))
