function x = RK4(IC, M, h, bigL, Nx, Nz, dx, km, every)
	x = IC;
	for tt = 1 : M
		x = RK4_step(x, h, bigL, Nx, Nz, dx, km);
		if mod(tt, every) == 1 
			if ismember(1, isnan(x)) || ismember(1, isinf(x))
				display(tt*h)
				break
			end
		end
	end
end

function update = RK4_step(x, h, bigL, Nx, Nz, dx, km)
	update = x;
	k = h * LandNL(x, bigL, Nx, Nz, dx, km);        %k1
	update = update + 1/6*k;
	k = h * LandNL(x + .5*k, bigL, Nx, Nz, dx, km); %k2
    update = update + 1/3*k;
    k = h * LandNL(x + .5*k, bigL, Nx, Nz, dx, km); %k3
    update = update + 1/3*k;  
    k = h * LandNL(x + k, bigL, Nx, Nz, dx, km);    %k4
    update = update + 1/6*k; 
end

function tend = LandNL(x, bigL, Nx, Nz, dx, km)
	tend = boxify3(bigL*flatten(x)) + NL(x, Nx, Nz, dx, km);
end
