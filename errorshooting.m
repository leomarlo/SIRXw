
wi =2; ki =2;
for bi = length(brange):-1:1
    disp(bi)
    bla = simRinf{wi,ki,bi};
    if std(bla) <0.000001
        idx = ones(length(bla),1);
        C=[bla(1),bla(1)];
%             [idx,C] = kmeans(bla',2);
    else
        [idx,C] = kmeans(bla',2);
    end
    [xmean,xsdev,weight] = corrected(bla,0.05,0.05)
end