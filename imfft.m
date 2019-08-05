function X = imfft(x,dim,mode)
if nargin<3
    mode=1;
end
if mode == 1
    X = fftshift(fft(fftshift(x,dim),[],dim),dim);
elseif mode == -1
    X = ifftshift(ifft(ifftshift(x,dim),[],dim),dim);
else
    error('fft: mode=1, ifft mode=-1');
end

