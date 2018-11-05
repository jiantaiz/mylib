function ph_bmf = get_shear_wave(ph_uw_fft,mask,time_steps)
% GET_SHEAR_WAVE ...
%    ph_bmf = GET_SHEAR_WAVE(ph_uw_fft,mask,time_steps) ...
%
%    Example:
%    ... 
%
%    Subfunctions: 
%    See also: 

% AUTHOR    : Yi Sui
% DATE      : 05/16/2017
%%[nx, ny, nz,m] = size(ph_uw_fft);N=3;for n = 0:N-1    phoff = 2*pi*n/N;    img(:,:,:,:,n+1) = ph_uw_fft.*exp(1i.*phoff);endimg = real(img);ph_bmf = zeros(nx,ny,nz,m,time_steps);for sl = 1:nz    bmf = remove_bulk_motion(img(:,:,sl,:,:),repmat(mask(:,:,sl,:),[1 1 1 1,size(img,5)]),'poly',3);    ph_bmf_fft = ifft(bmf,[],5);    if time_steps > 3        N=time_steps;        for n = 0:N-1            phoff = 2*pi*n/N;            ph_bmf(:,:,sl,:,n+1) = real(2.*ph_bmf_fft(:,:,:,:,2).*exp(1i.*phoff));        end    endend
