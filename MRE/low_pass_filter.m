function I = low_pass_filter(cimgs,M)
% LOW_PASS_FILTER ...
%    I = LOW_PASS_FILTER(cimgs,M) ...
%
%    Example:
%    ... 
%
%    Subfunctions: 
%    See also: 

% AUTHOR    : Yi Sui
% DATE      : 05/16/2017
%%
Kx = ifftshift(fft(fftshift(cimgs,1),[],1),1);
Kxy = ifftshift(fft(fftshift(Kx,2),[],2),2);
Kxyz = ifftshift(fft(fftshift(Kxy,3),[],3),3);

sz = size(Kxy);

% M=128;
w=hanning(sz(1)-M)';
w=[zeros(1,M/2) w zeros(1,M/2)];

[Wx, Wy] = meshgrid(w,w);
W = Wx.*Wy;

Kxy_fil = bsxfun(@times, Kxy,W);


Ix = ifftshift(ifft(fftshift(Kxy_fil,1),[],1),1);
I = ifftshift(ifft(fftshift(Ix,2),[],2),2);


% imdisp(angle(cat(2,I(:,:,1:10),cimgs(:,:,1:10))))
% imdisp(Kxy_fil(:,:,1))
