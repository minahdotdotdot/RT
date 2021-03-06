function [x, ES, FS] = RK4(x, M, h, every, L, dp);
	ES = []; FS = [];
	for tt = 1 : M
		x = RK4_step(x, h, L, dp);
		if tt == 1
            tic
        end
		if mod(tt, every) == 1
			toc
            %display(tt*h) 
            [ES(end+1), FS(end+1)] = computeE(x, dp);
			if ismember(1, isnan(x)) || ismember(1, isinf(x))
				%display(tt*h)
				break
			end
			tic
		end
	end
end

function update = RK4_step(x, h, L, dp)
	update = x;
	k = h * LandNL(x, L, dp);        %k1
	update = update + 1/6*k;
	k = h * LandNL(x + .5*k, L, dp); %k2
    update = update + 1/3*k;
    k = h * LandNL(x + .5*k, L, dp); %k3
    update = update + 1/3*k;  
    k = h * LandNL(x + k, L, dp);    %k4
    update = update + 1/6*k; 
end

function tend = LandNL(x, L, dp); 
	tend = boxify3(L*flatten(x), dp.NxNz) + NL(x, dp);
end
