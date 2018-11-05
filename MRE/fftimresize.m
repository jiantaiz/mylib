function im =  fftimresize(im,newsize)
% FFTIMRESIZE ...
%    im =  FFTIMRESIZE(im,newsize) ...
%
%    Example:
%    ... 
%
%    Subfunctions: 
%    See also: 

% AUTHOR    : Yi Sui
% DATE      : 05/16/2017
%%

sz = size(im);

newsize=[newsize, sz(3:end)];

zp = zeros(newsize);

IM = fftshift(fft2(im));

x0 = floor((newsize(1) - sz(1))/2);
y0 = floor((newsize(1) - sz(1))/2);
if x0>=0
    zp( x0+1:x0+sz(1), y0+1:y0+sz(2),: ) = IM(:,:,:);
else
    x0 = -x0;
    zp(:,:,:) = IM(x0+1:end-x0, y0+1:end-y0,: );
end
% if y0>=0
    
y0 = -y0;
zp(:,:,:) = IM(x0+1:end-x0, y0+1:end-y0,: );
im = ifft(ifftshift(zp));


% 
% 
% 
% IM = fft(IM,[],2);
