function [x, ES, FS]= ARK(x, M, h, every, dp, arks)
	ES = zeros(M/every+1,1); FS = zeros(M/every+1,1);kk = 0;
	% make this better
	invd = inv(1 - h*arks.d*L)
	%Initialize stages cell array.
	stages=repmat(flatten(x), 1, length(arks.b));
	for tt = 1 : M
		x = ARK_step(x, h, arks, dp, stages);
        if tt == 1
        end
		if mod(tt, every) == 1
			kk = kk+1; [ES(kk), FS(kk)] = computeE(x, dp);
			if ismember(1, isnan(x)) || ismember(1, isinf(x))
				break
			end
		end
	end
end

function update = ARK_step(x, h, arks, dp, stages)
	s = length(RK.b);
	PP=0;
	% ks, x, PP are stored as vectors (flattened) within this function.
	x = flatten(x);
	stages(:,1) = x;
	for ii = 2 : s
		ks[i,:] = invd .* (x + h*(
            L .* lincom(RKT.Ae(i-1,1:i-1), stages(:,1:i-1)
            + lincomN(RKT.Ai(i-1,1:i-1), stages(:,1:i-1), dp)
            )
        )
	end
	%new PP for update
	PP = h*lincomIF(RK.b, RK.cx(end,:), RK.cmat, stages);
	update = boxify3(...
		zhat + h*( L .* lincom(RKT.b, ks) + lincomN(RKT.b, ks, fP, k, NLfunc)),...
		dp.Nx, dp.Nz);
end


function tmp = lincomN(A, stages, dp)
    tmp = A(1)*NL(boxify3(stages(:,1)), dp)
    for i = 2 : length(A)
        tmp = tmp + A(i) * flatten(NL(boxify3(stages(:,i)), dp));
    end
    return tmp
end