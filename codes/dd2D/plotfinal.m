function plotfinal(finaldata, dp, bname, dname)
data = load(finaldata);
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
end