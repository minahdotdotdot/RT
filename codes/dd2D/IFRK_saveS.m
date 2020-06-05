function [x, ES, FS]= IFRK_saveS(x, M, h, every, Severy, name, dp, RK)
	ES = zeros(M/every+1,1); FS = zeros(M/every+1,1);kk = 0;
	run cmap
	dx=dp.x/dp.l_o;
	z=dp.z/dp.l_o;
	imk = 0;
	for tt = 1 : M
		x = IFRK_step(x, h, RK, dp);
		if mod(tt, every) == 1
			kk = kk+1; [ES(kk), FS(kk)] = computeE(x, dp);
            if mod(kk, Severy) ==1
			    imk = imk +1;
                display(imk)
                tth =  sprintf('%.2f', tt*h)
                Psi = boxify(x(:,1),dp.Nx,dp.Nz);
                ubox = real(ifft2(-1i*dp.mmm.*Psi));
                wbox = real(ifft2(1i*dp.kkk.*Psi));
                Tbox =real(ifft2(boxify(x(:,2), dp.Nx, dp.Nz)));
			    Sbox =real(ifft2(boxify(x(:,3), dp.Nx, dp.Nz)));
			    f=figure('position',[0,0,750,650])
                subplot(221)
                imagesc(dx,z,ubox);%, clims)
			    colormap(OrPu/255)
			    title(join(['Horizontal velocity at time ', tth],''))
			    colorbar
                subplot(222)
                imagesc(dx,z,wbox);%, clims)
			    colormap(OrPu/255)
			    title(join(['Vertical velocity at time ', tth],''))
			    colorbar
                subplot(223)
                imagesc(dx,z,Tbox);%, clims)
			    colormap(OrPu/255)
			    title(join(['T at time ',tth],''))
			    colorbar
                subplot(224)
                %clims = [-3,3];
			    imagesc(dx,z,Sbox);%, clims)
			    colormap(OrPu/255)
			    title(join(['S at time ',tth],''))
			    colorbar
			    saveas(gcf,join(['../../plots/',name,sprintf('%04d.png',imk)],''));
			    close(f)
			end
			if ismember(1, isnan(x)) || ismember(1, isinf(x))
                ES = ES(1:kk-1); FS=FS(1:kk-1);
				break
			end
		end
	end
end

function update = IFRK_step(x, h, RK, dp)
	s = length(RK.b);
	stages = cell(s,1); PP=0;
	% ks, x, PP are stored as vectors (flattened) within this function.
	stages{1} = flatten(NL(x, dp));
	x = flatten(x);
	for ii = 2 : s
		PP = h*lincomIF(RK.A(ii, 1:ii-1), RK.cx(ii, 1:ii-1), RK.cmat, stages(1:ii-1));
		stages{ii} = flatten(NL(boxify3(RK.cmat(RK.cx(ii,ii))*x + PP, dp.NxNz), dp));
	end
	%new PP for update
	PP = h*lincomIF(RK.b, RK.cx(end,:), RK.cmat, stages);
	update = forcereal(boxify3(RK.cmat(RK.cx(end-1,end))*x+PP, dp.NxNz), dp);
end

%lincom outputs a flat vector
function lincom = lincomIF(A, cx, cmat, k)
	lincom = A(1)*cmat(cx(1))*k{1};
	for ii = 2 : length(k)
		lincom = lincom + A(ii)*cmat(cx(ii))*k{ii};
	end
end
