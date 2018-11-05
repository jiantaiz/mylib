function wout = slim_harmonics(w, n, delay,mask)
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
if size(w,4) > 2
    w_f=ifft(w,[],4).*2;
end

for j = 1:abs(n)
    wout(:,:,:,j,1) = real(exp(-1i*(j-1)*2*pi/n - 1i*delay).*w_f(:,:,:,2));
    wout(:,:,:,j,2) = real(exp(-1i*(j-1)*2*pi/n - 1i*delay).*w_f(:,:,:,3));
    wout(:,:,:,j,3) = real(exp(-1i*(j-1)*2*pi/n - 1i*delay).*w_f(:,:,:,4));
end


if ~isempty(mask)
    wout = bsxmul(wout,mask);
end