function genRationalApprox(name)
    nf = load(join(['../data/Lhc_',name,'.mat'],''));
    scheme=nf.scheme;
    h=nf.h; deg=nf.deg;
    file = load(join(['../data/',name,'.mat'],''));
    crat = file.crat;
    crat2 = file.crat;
    s = size(file.x);
    uniquex = unique(file.x);
    L = [];
    for i = 1 : length(uniquex)
        L = [L h*uniquex(i)*file.L];
    end
    Z = exp(L);
    [r, pol, res, zer, zj, fj, wj] = aaa(Z,L,'degree', deg);
    c = containers.Map('KeyType', 'double', 'ValueType', 'any');
    c2 = containers.Map('KeyType', 'double', 'ValueType', 'any');
    for i = 1 : length(uniquex)
        L = h*uniquex(i)*file.L;
        c(uniquex(i)) = r(L) - r(0) +1;
        c2(uniquex(i)) = r(L);
    end

    for i = 1 : s(1)
        for j = 1: s(1)-1
            crat(i,j) = {c(file.x(i,j))};
            crat2(i,j) = {c2(file.x(i,j))};
        end
    end
    save(join(['../data/',name,'R.mat'],''), 'crat');
    save(join(['../data/',name,'R2.mat'],''), 'crat2');

    %{
    c = containers.Map('KeyType', 'double', 'ValueType', 'any');
    c2 = containers.Map('KeyType', 'double', 'ValueType', 'any');
    for i = 1 : length(uniquex)
        L = h*uniquex(i)*file.L;
        Z = exp(L);
        [r, pol, res, zer, zj, fj, wj] = aaa(Z,L,'degree', deg);
        c(uniquex(i)) = r(L) - r(0) +1;
        c2(uniquex(i)) = r(L);
    end

    for i = 1 : s(1)
        for j = 1: s(1)-1
            crat(i,j) = {c(file.x(i,j))};
            crat2(i,j) = {c2(file.x(i,j))};
        end
    end
    save(join(['../data/',name,'R.mat'],''), 'crat');
    save(join(['../data/',name,'R2.mat'],''), 'crat2');
%}
end
