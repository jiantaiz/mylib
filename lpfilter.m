function cimgs = lpfilter(cimgs,dim,alpha)


if nargin <2
    dim = 3;
end

if nargin <3
    alpha = 3;
end

if numel(dim)>1
    for k=1:numel(dim)
        cimgs = lpfilter(cimgs,dim(k),alpha);
    end
    return;
end
cimgs(isnan(  cimgs))=0;
CIMGS = fftshift(fft(fftshift(cimgs,dim),[],dim),dim);

sz = size(cimgs);
% gw = gausswin(sz(dim),alpha);


N = sz(dim)-1;
n = (0:N)'-N/2;
gw = exp(-(1/2)*(alpha*n/(N/2)).^2);



if dim>1
    gw2 = reshape(gw,[ones(1,dim-1), sz(dim)]);
    CIMGS = bsxfun(@times, CIMGS, gw2);
else
    CIMGS = bsxfun(@times, CIMGS, gw); % to make coder happy
end

cimgs = ifftshift(ifft(ifftshift(CIMGS,dim),[],dim),dim);