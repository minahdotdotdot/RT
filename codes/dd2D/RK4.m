function x = RK4(IC, M, h, bigL, Nx, Nz, km, kk, mm, every)
	x = IC;
        [Psi, T, S] = boxify3NL(x);
	for tt = 1 : M
		x = RK4_step(x, h, bigL, Nx, Nz, km, kk, mm);
		if mod(tt, every) == 1
                        display([norm(x(:,1)) norm(x(:,2)) norm(x(:,3))]) 
			if ismember(1, isnan(x)) || ismember(1, isinf(x))
				display(tt*h)
				break
			end
		end
	end
end

function update = RK4_step(x, h, bigL, Nx, Nz, km, kk, mm)
	update = x;
	k = h * LandNL(x, bigL, Nx, Nz, km, kk, mm);        %k1
	update = update + 1/6*k;
	k = h * LandNL(x + .5*k, bigL, Nx, Nz, km, kk, mm); %k2
    update = update + 1/3*k;
    k = h * LandNL(x + .5*k, bigL, Nx, Nz, km, kk, mm); %k3
    update = update + 1/3*k;  
    k = h * LandNL(x + k, bigL, Nx, Nz, km, kk, mm);    %k4
    update = update + 1/6*k; 
end

function tend = LandNL(x, bigL, Nx, Nz, km, kk, mm)
	tend = boxify3(bigL*flatten(x)) + NL(x, Nx, Nz, km, kk, mm);
end
