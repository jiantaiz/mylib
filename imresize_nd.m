function x= imresize_nd(A,num,dim,varargin)
% IMRESIZE_ND ...
%    x= IMRESIZE_ND(A,num,dim,varargin) ...
%
%    Example:
%    % resize the 3rd dimention to 128 pixels
%    x= IMRESIZE_ND(A,128,3);
%    % resize dimentions 1 2 3  to [128 128 96] pixels
%    x= IMRESIZE_ND(A,[128 128 96]);
%    % double the size on dim [1 2 3]
%    x= IMRESIZE_ND(A,2,[1 2 3]);
%    Subfunctions: 
%    See also: imresize

% AUTHOR    : Yi Sui
% DATE      : 10/31/2018
%%

if numel(num)>1
    x=A;
    for k=1:numel(num);
        x = imresize_nd(x,num(k),k,varargin{:});
    end
    return;
end

if numel(dim)>1
    sz=size(A);
    newsz = sz(dim) .* num; 
    x=A;
    for k = 1:numel(dim)
        x = imresize_nd(x,newsz(k),dim(k),varargin{:});
    end
    return;
end

if num == size(A,dim) % same image size, do nothing.
    x=A;
    return;
end

[x,perm,nshifts]  = shiftdata(A,dim);


[Nrow,Ncol,Nz]=size(x);

x = imresize (x,[num, Ncol],varargin{:});
x = unshiftdata(x,perm,nshifts);
