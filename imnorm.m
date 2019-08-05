function A = imnorm(A,pct)
if nargin<2
    pct = 0.05;
end
num = numel(A);

idx = round(num.*pct);
idx = max(5,idx);

tmp = sort(A(:),'descend');
tmp = tmp((~isnan(tmp) & ~isinf(tmp)));
scale = mean(tmp(2:idx));
A = A./scale;