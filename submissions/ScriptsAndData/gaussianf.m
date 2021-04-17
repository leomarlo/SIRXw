function y = gaussianf(x,mu,sd)
    y = 1/(2*pi*sd)*exp(-(x-mu).^2/(2*sd^2));
end