function w = first_harmonic(w, n, delay,mask)
%FIRST_HARMONIC     reutrns first harmonic component of a wave.
%   w_1 = FIRST_HARMONIC(w, n) is the first harmonic component of 
%   w(with a size of [nx ny nz nt]) interplated to n time points, 
%   thus w_1 has a size of [nx ny nz n];
% by Yi Sui
w=squeeze(w);
if nargin <4
    mask=[];
end
if nargin <3
    delay=0;
end
if nargin <2
    n=16;
end
if size(w,5)>1
    w=[w(:,:,:,:,1),w(:,:,:,:,2),w(:,:,:,:,3)];
%     w=[w(:,:,:,:,1);w(:,:,:,:,2);w(:,:,:,:,3)];

end
freq = 2; % first harmonic
% freq = 4; % second harmonic
if size(w,4) > 2
    w_f=ifft(w,[],4); 
    w_f =(w_f(:,:,:,freq))*2;
elseif size(w,4) == 2
    w_f = w(:,:,:,1) + 1i*w(:,:,:,2);
else
    if isreal(w);
        warning('expecting complex number for w');
    else
    w_f = w*2;
    end
end

for j = 1:abs(n)
    w(:,:,:,j) = real(exp(-1i*(j-1)*2*pi/n - 1i*delay).*w_f(:,:,:,1));
end


if ~isempty(mask)
    w = bsxmul(w,mask);
end