function img_out = rm_linear_ph(img,a,b,phi)
%remove li,near phase on a image by shiftting k-space
if nargin <3; phi=0;end
img_out = rm_linear_ph_helper(img,a,b);

img_out = img_out .* exp(1i.*phi);



function out = rm_linear_ph_helper(m,a,b)
    M = fftshift(fft2(fftshift(m)));
    
    M = circshift(M,a,1);
    M = circshift(M,b,2);
    
    out = ifftshift(ifft2(ifftshift(M)));
    
    


