function genRationalApprox(name)
    nf = load(join(['../data/Lhc_',name,'.mat'],''));
    scheme=nf.scheme;
    h=nf.h; deg=nf.deg;
    file = load(join(['../data/',name,'.mat'],''));
    crat = file.crat;
    s = size(file.x);
    for i = 1 : s(1)
        for j = 1: s(1)-1
            Z = exp(h*file.x(i,j)*file.L);
            [r, pol, res, zer, zj, fj, wj] = aaa(Z,file.L,'degree', deg);
            crat(i,j) = {r(h*file.x(i,j)*file.L)};
        end
    end
    save(join(['../data/',name,'R.mat'],''), 'crat');
end
