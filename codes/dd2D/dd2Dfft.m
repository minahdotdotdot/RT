function xICfft = dd2Dfft(xIC, dp);
    [Psi, T, S] = boxify3NL(xIC);
    Psi =vectorize(fft2(boxify(Psi, dp.Nx, dp.Nz)),dp.NxNz);
    T =vectorize(fft2(boxify(T, dp.Nx, dp.Nz)),dp.NxNz);
    S =vectorize(fft2(boxify(S, dp.Nx, dp.Nz)),dp.NxNz);
    xICfft=[Psi T S];
end
