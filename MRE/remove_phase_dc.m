function ph=remove_phase_dc(pdimg,dim)
% REMOVE_PHASE_DC remove dc component from MRE phase different images along the time dimension
%    ph=REMOVE_PHASE_DC(pdimg,dim) pdimg is complex-valued. dim specifies
%    the time dimension, May fail when dynamic range > 2pi
%
%    See also: 

% AUTHOR    : Yi Sui
% DATE      : 05/16/2017

img_mean = mean(pdimg,dim);
ph=angle(bsxfun(@times,pdimg,conj(img_mean)));
