clear all
data = load('../../data/IFRK3h05.mat');
[Psi, T, S] = boxify3NL(data.xf);

run dd2D
Psibox =real(ifft2(boxify(Psi, dp.Nx, dp.Nz)));
Tbox =real(ifft2(boxify(T, dp.Nx, dp.Nz)));
Sbox =real(ifft2(boxify(S, dp.Nx, dp.Nz)));

subplot(131)
contour(dp.x,dp.z,Psibox)
title('Psi')
colorbar
subplot(132)
contour(dp.x,dp.z,Tbox)
title('T')
colorbar
subplot(133)
contour(dp.x,dp.z,Sbox)
title('S')
colorbar

