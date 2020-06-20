function [av,bs] = movingaverage(estimates,weights,brange,betaN,sigmas)
    bs = brange(1):(brange(end)-brange(1))/(betaN-1):brange(end);
    av = ones(1,length(bs));
    dbeta = brange(2)-brange(1);
    i = 1;
    for x=bs
        modulation = zeros(1,length(estimates));
        for j = 1:length(estimates)
            modulation(j) = gaussianf(brange(j),x,sigmas * dbeta) * weights(j);
        end
        modulation = modulation/sum(modulation);
        av(i)= modulation * estimates';  
        i=i+1;
    end
end