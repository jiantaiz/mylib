function ortho_view(img,voxel_size, clim,opt)
% ortho_view(img,voxel_size,clim,opt)
% orthographic view of 3D or 4D images.


if iscell(img)
    img  = cat(2,img{:});
else
    img = img(:,:,:,:);
end
if nargin <2 || isempty(voxel_size)
    voxel_size = [1 1 1];
end
if nargin < 3
    clim = [];
end
if nargin <4
    opt = [];
end

[nx,ny,nz,nt] = size(img);
origin = floor([nx ny nz]/2);
datatype = 64;%float64
description= '';
img = rot90(img,-1); %rotate 90 deg clockwise to address the coordinate difference b/n matlab and nii.
nii = make_nii(img, voxel_size, origin, datatype, description);

opt.glblocminmax = clim;
view_nii (nii,opt)
addRoiToolbar;