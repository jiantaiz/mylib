function T = tile(M, dim, dir)
% MOSAIC tiles image stack M in a row/col.
%    T = tile(M, dim, dir) 
%       M: image stack e.g a 3D volume with 10 slices size(M) = [32 32 10]
%       dim: the dimention of M to tile.
%       dir: 1 - colum, 2 - row (default)
%    Example:
%    ... 
%
%    Subfunctions: 
%    See also: 

% AUTHOR    : Yi Sui
% DATE      : 05/16/2017

if nargin < 3
    dir = 2;
end
sz = size(M);

for k = 1:sz(dim)
    idx = repmat({':'}, 1, numel(sz));
    idx{dim} = k;
    C{k} = M(idx{:}); 
end

T = cat(dir, C{:});
