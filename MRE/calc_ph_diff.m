function [ph, ph_uw, dimg_cmb] = calc_ph_diff(cimgs,dualsat)
%calculate phase difference images
%cimgs(nx,ny,nz,nt,ndir)

if nargin<2
    dualsat=0;
end

if ischar(cimgs)
    load(cimgs,'cimgs')
end

dimg = cimgs(:,:,:,:,1:2:end).* conj(cimgs(:,:,:,:,2:2:end));
sz = size(dimg);
% ph = angle(dimg);
% dimg = sqrt(abs(dimg)).*exp(1i.*ph);
if dualsat == 1
    dimg_cmb = dimg(:,:,:,1:sz(4)/2,:) + dimg(:,:,:,sz(4)/2+1:end,:);
elseif dualsat == 2
    dimg_cmb = dimg(:,:,:,1:2:end,:) + dimg(:,:,:,2:2:end,:);
else
    dimg_cmb = dimg;
end
ph = angle(dimg_cmb);
if nargout>1
    ph_uw = unwrap3D_ssh(ph,1);
% ph_uw=[];
end

if nargout>2
    dimg_cmb = sqrt(abs(dimg_cmb)).*exp(1i.*ph);
end


