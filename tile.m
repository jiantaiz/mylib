function T = tile(M, dim, dir)
% T = tile(M, dim, dir)
% tile the dim dimension of image stack M along row or col.
% M: image stack e.g a 3D volume with 10 slices size(M) = [32 32 10]
% dim: tile the dim dimention of M.
% dir: col (1) or row (2)
if nargin < 3
    dir = 2;
end
sz = size(M);

% [M,perm,nshifts] = shiftdata(M,dim);


for k = 1:sz(dim)
    idx = repmat({':'}, 1, numel(sz));
    idx{dim} = k;
    C{k} = M(idx{:}); 
end

T = cat(dir, C{:});



