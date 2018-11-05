function [phs_uw,cimgs,phs,phs_est,magn_w,magn_f,frac_f] = unwrap_fat_water_dual_menc(cimgs_w,cimgs_f,T,varargin)
% UNWRAP_FAT_WATER_DUAL_MENC unwrap skull+brain images acquired using two scans at either water
% or fat peak. the final cimgs are combined with least squares 
%    phs_uw = UNWRAP_FAT_WATER_DUAL_MENC(cimgs_w,cimgs_f,T)
%    [phs_uw,cimgs,phs,phs_est,magn_w,magn_f,frac_f] = UNWRAP_FAT_WATER_DUAL_MENC(cimgs_w, cimgs_f,T,N_medfil,high_vib,z_filter) ...
%
%    See also: unwrap_dual_menc

% AUTHOR    : Yi Sui
% DATE      : 07/05/2017
%%
if isempty(cimgs_f) % no fat images
    cimgs = cimgs_w;
else %combine fat and water
    magn_w = mean(abs(cimgs_w(:,:,:,:)),4);
    magn_f = mean(abs(cimgs_f(:,:,:,:)),4);
    frac_f = magn_f./(magn_f+magn_w); %fat fraction
    
    %least squares combination
    % cimgs = (bsxfun(@times,cimgs_f,frac_f) ...
    %       +  bsxfun(@times,cimgs_w,1-frac_f) );
    % cimgs = bsxfun(@rdivide,cimgs,frac_f.^2+(1-frac_f).^2);
    
    % cimgs = cimgs_w + cimgs_f;
    %unwrap
    
    cimgs = abs(cimgs_w).*cimgs_w + abs(cimgs_f).*cimgs_f;
    cimgs = sqrt(abs(cimgs)).*exp(1i*angle(cimgs));
    
end
[phs_uw, phs, phs_est]  = unwrap_dual_menc(cimgs,T,varargin{:});
