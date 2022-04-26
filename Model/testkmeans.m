close all; clear all;

function [x2, n, b] = compute_xpdf(x)
  x2 = reshape(x, 1, prod(size(x)));
  [n, b] = hist(x2, 40);
  % This is definitely not probability density function
  x2 = sort(x2);
  % downsampling to speed up computations
  x2 = interp1 (1:length(x2), x2, 1:1000:length(x2));
end

nboot = 500;
sample_size = [256 256];

% Unimodal
sample2d = normrnd(0.0, 10.0, sample_size);

[xpdf, n, b] = compute_xpdf(sample2d);
% [dip, p_value, xlow, xup] = HartigansDipTest(xpdf); 
[dip, p_value, xlow,xup]=HartigansDipSignifTest(xpdf,nboot)

figure;
subplot(1,2,1);
bar(n, b)
title(sprintf('Probability of unimodal %.2f', p_value))

% Bimodal
sample2d = sign(sample2d) .* (abs(sample2d) .^ 0.5);

[xpdf, n, b] = compute_xpdf(sample2d);
[dip, p_value, xlow, xup] = HartigansDipSignifTest(xpdf, nboot); 

subplot(1,2,2);
bar(n, b)
title(sprintf('Probability of unimodal %.2f', p_value))

print -dpng modality.png