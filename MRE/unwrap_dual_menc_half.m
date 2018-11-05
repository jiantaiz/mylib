function [phs_uw, phs] = unwrap_dual_menc_half(cimgs,T,N_medfil,high_vib,z_filter)
% UNWRAP_DUAL_MENC unwrap phase maps from dual menc aquisition
%    [phs_uw, phs, phs_est] = UNWRAP_DUAL_MENC(cimgs,T,N_medfil,high_vib,z_filter)
%    INPUT:
%       cimgs(x,y,z,t,dir): contains complex images.
%       T:  3x3 the transform matrix from low menc to high menc
%       i.e.  phs_hi = T*phs_low
%       N_medfil: number of times to apply median filter.
%    OUTPUT:
%       phs_uw: unwrapped phase maps
%       phs: wrapped phase maps
%       phs_est: estimated phase maps from low menc
%    Limitation: only work with 4 time offsets acq. calc_ph_first_harmonic
%    halves dynamic range.
%    See also: calc_menc

% AUTHOR    : Yi Sui
% DATE      : 06/09/2017
%%
if nargin<3
    N_medfil=0;
end
if nargin<4
    high_vib = [0 0 0];
end

if numel(high_vib) == 1
    high_vib = [high_vib high_vib high_vib];
end
if nargin<5
    z_filter= 0;
end

phs = angle(cimgs(:,:,:,:,1:2:end).*conj(cimgs(:,:,:,:,2:2:end)));
sz  = size(cimgs);
big_ph_est = complex( zeros( [sz(1:3) 3]));

%x dir
dir = 1;
s1 = cimgs(:,:,:,:,(dir-1)*2+1);
s2 = cimgs(:,:,:,:,(dir-1)*2+2);
ph = (s1.*s2.*conj(cimgs(:,:,:,:,5)).*conj(cimgs(:,:,:,:,6))); % x-z
ph  = calc_ph_first_harmonic(ph,[],4);
phs_low(:,:,:,1,dir) = real(ph);
phs_low(:,:,:,2,dir) = imag(ph);

ph = calc_ph_first_harmonic(s1,s2,4);
phs_high(:,:,:,1,dir) = real(ph);
phs_high(:,:,:,2,dir) = imag(ph);

%y dir
dir = 2;
s1 = cimgs(:,:,:,:,(dir-1)*2+1);
s2 = cimgs(:,:,:,:,(dir-1)*2+2);
ph = (s1.*s2.*conj(cimgs(:,:,:,:,5)).*conj(cimgs(:,:,:,:,6))); % x-z
ph  = calc_ph_first_harmonic(ph,[],4);
phs_low(:,:,:,1,dir) = real(ph);
phs_low(:,:,:,2,dir) = imag(ph);

ph = calc_ph_first_harmonic(s1,s2,4);
phs_high(:,:,:,1,dir) = real(ph);
phs_high(:,:,:,2,dir) = imag(ph);


r = real(T(1,1));%menc ratio , only works when half integer e.g 7.5,8.5, etc.

N1 = (phs_low * r - phs_high)/2/pi;
N1rd = round(N1);
%if small menc already wrapped. 
N2 = N1 + (N1>0)*(-r) + (N1<=0)*r;
N2rd = round(N2);

N = (abs(N1 - N1rd) < abs(N2 - N2rd) ).* N1rd + (abs(N1 - N1rd) >= abs(N2 - N2rd) ).* N2rd;
phs_uw = phs_high + N.*2*pi;
% 
% phs_z_est = cimgs(:,:,:,:,5).*cimgs(:,:,:,:,6);
% 
% if high_vib(3)
%     
%     phs_z_est = angle(phs_z_est);
%     phs_z_est = phs_z_est - phs_z_est(:,:,:,[1 1:end-1],:);
%     phs_z_est = unwrap(phs_z_est,[],4);
%     big_ph_est(:,:,:,3) = phs_z_est(:,:,:,1)-phs_z_est(:,:,:,3) + 1i*(phs_z_est(:,:,:,2)-phs_z_est(:,:,:,4));
%     big_ph_est(:,:,:,3) = big_ph_est(:,:,:,3).* conj(T(3,3));
% else
%     big_ph_est(:,:,:,3) = calc_ph_first_harmonic(phs_z_est,[],4).* conj(T(3,3));
%     
% end
% 
% 
% if z_filter == 1
% %     big_ph_est = lpfilter(big_ph_est,3,2.5);
% elseif z_filter == 2
% %     big_ph_est = ordfilt3D_2(big_ph_est,14);
% elseif z_filter == 3
% %     big_ph_est = lpfilter(big_ph_est,3,4);
% %     big_ph_est = ordfilt3D_2(big_ph_est,14);
%     big_ph_est = medfilt3_nd(big_ph_est);
% end
%     phs_hi_est = big_ph_est;
% 

% phs_hi_est = big_ph_est;
% phs_est = zeros([sz(1:3) 4, 3]);
% for k = 1:4
%     phs_est(:,:,:,k,:) = real(0.5 * phs_hi_est.*exp(1i.*-(k-1)*pi/2));
% end
% phs_est = real(phs_est);
% tol = 1;
% phs_uw = myunwrap(phs,phs_est, 1,[],tol);
for n=1:N_medfil
    if verLessThan('matlab','9.1.0')
        phs_uw_mf = ordfilt3D_2(phs_uw,14);
    else
        phs_uw_mf = medfilt3_nd(phs_uw);
    end
    phs_uw = myunwrap(phs,phs_uw_mf, 1,[],tol);
end


