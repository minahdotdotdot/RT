function final = exRK(IC, A, b, c, M, h, bigL, Nx, Nz, dx, km, every)
	x = IC
	for tt = 1 : M
		x = exRK_step(x, A, b, c, h, bigL, Nx, Nz, dx, km)
		if mod(tt, every) == 1 
			if ismember(1, isnan(x)) || ismember(1, isinf(x));
				break
			end
		end
	end
end


function update = exRK_step(x, A, b, c, h, bigL, Nx, Nz, dx, km)
	stages = cell(length(b),1);
	stages(1) = {LandNL(x, bigL, Nx, Nz, dx, km)};
	for kk = 2 : length(b)
		y = x;
		for ll = 1 : kk-1
			y = y + h*A(kk,ll)*stages{ll};
		end
		stages(kk) = {LandNL(y, bigL, Nx, Nz, dx, km)};
	end
	update = x;
	for kk = 1 : length(b)
		update = update + h*b(kk)*stages{kk};
	end
end

function tend = LandNL(x, bigL, Nx, Nz, dx, km)
	tend = boxify3(bigL*flatten(x)) + NL(x, Nx, Nz, dx, km);
end
