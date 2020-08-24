using DelimitedFiles
include("readwrite.jl")
na="IFRK3-000100"
sol = readCfile(na);
N = 2^13;
totenergy = sum(sol .* conj.(sol), dims=2)/N^2;
semilogy(0:length(totenergy),totenergy)

