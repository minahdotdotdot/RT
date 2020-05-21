function [x, ES, FS]= ARK(x, M, h, every, dp, arks, L, invd)
	ES = zeros(M/every+1,1); FS = zeros(M/every+1,1);kk = 0;
	%Initialize stages cell array.
	stages=repmat(flatten(x), 1, length(arks.b));
	for tt = 1 : M
		x = ARK_step(x, h, arks, dp, stages, L, invd);
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

function update = ARK_step(x, h, arks, dp, stages, L, invd)
	s = length(arks.b);
	PP=0;
	% ks, x, PP are stored as vectors (flattened) within this function.
	x = flatten(x);
	stages(:,1) = x;
	for ii = 2 : s
		stages(:,ii) = invd * (x + h*(...
            L * sum(arks.Ae(ii-1,1:ii-1).*stages(:,1:ii-1),2)... %matlab does this fast.
            + lincomN(arks.Ai(ii-1,1:ii-1), stages(:,1:ii-1), dp)...
            )...
        );
    end
	update = boxify3(...
		x + h*(L*sum(arks.b.*stages,2) + lincomN(arks.b, stages, dp)),...
		dp.NxNz);
end


function tmp = lincomN(A, stages, dp)
    tmp = A(1)*flatten(NL(boxify3(stages(:,1),dp.NxNz),dp));
    for ii = 2 : length(A)
        tmp = tmp + A(ii) * flatten(NL(boxify3(stages(:,ii), dp.NxNz), dp));
    end
end
