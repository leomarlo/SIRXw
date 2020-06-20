
function JB = JarqueBera(X)
    n = prod(size(X));
    S = skewness(X(:));
    K = kurtosis(X(:));
    JB = (n/6)*(S^2 + (1/4)*((K-3)^2));
end