include("readwrite.jl")
α=1/2;
β=0;
zhat = readCfile("IFRK3-000100");
include("setupMMT.jl")
E = zeros(size(zhat)[1])
HNL = zeros(size(zhat)[1])
for i = 1 : size(zhat)[1]
	E[i] = norm(ifft((abs.(k).^(α/2)) .*zhat[i,:],2))^2;
	HNL[i] = norm(ifft((abs.(k).^(β/4)) .*zhat[i,:],2))^4;
end
ϵ = E ./ HNL;