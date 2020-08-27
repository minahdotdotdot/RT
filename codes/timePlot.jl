using DelimitedFiles
include("readwrite.jl")
na="IFRK3-000100"
sol = readCfile(na);
N = 2^13;
totenergy = Float64.(sum(sol .* conj.(sol), dims=2)/N^2);
#semilogy(0:length(totenergy)-1,totenergy, label="h=0.0001")
T=200;
semilogy(0:T-1,totenergy[1:T], label="h=0.0001")
axvline(43, c=:red, label="end of transient")
ylabel("Total Energy")
xlabel("Time")
title("Total energy for MMT solved with IFRK3")
legend()
