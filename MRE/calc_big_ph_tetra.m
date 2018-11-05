function [ph_big, pdimg] = calc_big_ph_tetra(cimgs,mode);
% CALC_BIG_PH_TETRA calculate phase images of a tetra MRE scan
%    [ph_big, pdimg] = CALC_BIG_PH_TETRA(cimgs,mode); 
%    See also: 

% AUTHOR    : Yi Sui
% DATE      : 05/16/2017
%%
if strcmpi(mode ,'balanced') || mode == 1
    % X Y Z;X -Y Z; -X Y -Z; -X -Y Z
    
    dir1 = cimgs(:,:,:,:,1);
    dir2 = cimgs(:,:,:,:,2);
    dir3 = cimgs(:,:,:,:,3);
    dir4 = cimgs(:,:,:,:,4);
    ph_big(:,:,:,:,1) = dir1.*dir2.*conj(dir3).*conj(dir4);
    ph_big(:,:,:,:,2) = dir1.*conj(dir2).*dir3.*conj(dir4);
    ph_big(:,:,:,:,3) = dir1.*conj(dir2).*conj(dir3).*dir4;
 
    if (nargout > 1) pdimg = ph_big; end
        
    ph_big = angle(ph_big);
    
end
