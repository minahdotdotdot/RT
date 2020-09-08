include("readwrite.jl")
#=zhat = readCfile("IFRK3-000100");
include("setupMMT.jl")
E = zeros(size(zhat)[1])
HNL = zeros(size(zhat)[1])
for i = 1 : size(zhat)[1]
	E[i] = sum(abs.(ifft((abs.(k).^(fP.α/2)) .*zhat[i,:])).^2);
	HNL[i] = .5*fP.λ*sum(abs.(ifft((abs.(k).^(fP.β/4)) .*zhat[i,:])).^4);
end
ϵ = HNL ./ E;
=#
ϵ=readdlm("../txtfiles/eps.txt");
#=using PyPlot, LaTeXStrings
plot(1:length(ϵ),ϵ)
axvline(43, c=:red, label="end of transient")
legend()
title(L"\epsilon=L/H_{NL}"*", ground truth solution")
=#