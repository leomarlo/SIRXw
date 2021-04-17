function [x_mean,x_std,weight,upperones,condition] = corrected(X,thr_ku, thr_rho)
    if std(X) <0.000001
        idx = 2*ones(length(X),1);
        C=[X(1),X(1)];
        Ku =0;
%             [idx,C] = kmeans(bla',2);
    else
        [idx,C] = kmeans(X',2);
        Ku = kurtosis(X)*(std(X)^2);
    end
        
    if C(1)>C(2)
        idx = 3 - idx;
    end
    upperones = idx==2;
    rho = sum(upperones)/length(idx);
    condition = Ku>thr_ku && rho > thr_rho && sum(upperones)>2;
    if Ku>thr_ku && rho > thr_rho && sum(upperones)>2
        x_mean = mean(X(upperones));
        x_std = std(X(upperones));
        weight = rho;
    else
        x_mean = mean(X);
        x_std = std(X);
        weight = 1;
    end
end