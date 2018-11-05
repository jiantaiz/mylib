function [phs_uw, phs, phs_est, mag, dimg]= unwrap_dsat(matfile,cp_file,g,fat_phase_correction)
% UNWRAP_DSAT ...
%    [phs_uw, phs, mag, dimg]= UNWRAP_DSAT(matfile,cp_file,g,fat_phase_correction) 
%     unwrap dual saturation aquisition, the first 4 time points are water,
%     the next 4 time points are fat. matfile contains cimgs and dinfos.
%
%    Example:
%    ... 
%
%    Subfunctions: 
%    See also: 

% AUTHOR    : Yi Sui
% DATE      : 04/10/2018

if nargin<3
    g=0.777;
end
if nargin<4
    fat_phase_correction=1;
end

if isstruct(matfile)
    cimgs = matfile.cimgs;
    dinfos = matfile.dinfos;
else
    load(matfile,'cimgs','dinfos');
end
te = dinfos(1).EchoTime;% get echo time
freq_motion = dinfos(1).Private_0019_10bc;
b0 = dinfos(1).MagneticFieldStrength;
clear dinfos;
%calculate menc
meg_num =1;
[T,menc,menc_mtx]=calc_dual_menc_tmat(cp_file,g,freq_motion,te,meg_num);

N_medfilt =1; % number of times to apply median filter.
high_vib = [1 1 0]; % x-,y-,z- high vibration flag.
z_filter = 0
brain_mask2 = [];

h2o = cimgs(:,:,:,1:end/2,:);
fat = cimgs(:,:,:,end/2+1:end,:);

mag_fat = mean(abs(fat(:,:,:,:)),4);
mask_fat = create_mask(mag_fat,1);

tmp = bsxfun(@times,fat.*conj(h2o),mask_fat);
% ov(angle(tmp))
tmp2=tmp(abs(tmp)>0);
ph_diff = angle(mean(tmp2)); % constant phase difference between water and fat due to acquision.

clear tmp tmp2 mask_fat mag_fat


if b0 == 3 %3T
    f_fat = 440;
elseif b0 == 1.5
    f_fat = 220;
else
    f_fat = b0*3.5*42.58;
end
if fat_phase_correction == 1
%     cimgs_cmb = abs(h2o).*h2o + abs(fat).*fat.*exp(1i.*(te*1e-3*f_fat + pi)); %fat phase correction 
    cimgs_cmb = abs(h2o).*h2o + abs(fat).*fat.*exp(1i.*-ph_diff); %fat phase correction 
    ph_diff
    angle(exp(1i.*(te*1e-3*f_fat)))
    angle(exp(1i.*(te*1e-3*f_fat+pi)))
else
    cimgs_cmb = abs(h2o).*h2o + abs(fat).*fat; %fat phase correction
end
[phs_uw0, phs0, phs_est]  = unwrap_dual_menc(cimgs_cmb,T,N_medfilt,high_vib,z_filter);

dimgs_cmb = h2o(:,:,:,:,1:2:end).*conj(h2o(:,:,:,:,2:2:end)) + fat(:,:,:,:,1:2:end).*conj(fat(:,:,:,:,2:2:end)); % combine phase difference 
phs = angle(dimgs_cmb);
phs_uw = myunwrap(phs,phs_uw0,1);
if nargout >3
    mag = sqrt(abs(dimgs_cmb));
end
if nargout >4
    dimg = mag.*exp(1i*phs);
end
