function ov(img,voxel_size, clim,opt)
% ortho_view(img,voxel_size,clim,opt)
% orthographic view of 3D or 4D images.


if iscell(img) % interpolate to the same size 
    for k=1:numel(img)
        sz(k,1) = size( img{k},1);
        sz(k,2) = size( img{k},2);
        sz(k,3) = size( img{k},3);
    end
    newsz = max(sz);
    for k=1:numel(img)
        img{k} = single(img{k});
        if any(newsz ~= sz(k,:))
            img{k} = imresize_nd(img{k},newsz(1:3));
        end
    end
    tile_size = [size(img,2) size(img,1) size(img,3)];
%     img  = cat(2,img{:});
    img = cell2mat(img);
    
else
    img = img(:,:,:,:);
    tile_size = [];
end
if strcmpi(class(img),'double') || strcmpi(class(img),'logical');
    img = single(img);
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

opt.tile_size = tile_size;
opt.glblocminmax = clim;
view_nii (nii,opt)
addRoiToolbar;