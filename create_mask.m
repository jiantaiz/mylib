function mask = create_mask (magn,k)
% create mask for 3D volume
if nargin<2
    k=1.2;
end
sz = size(magn);
[nx, ny, nz] = size(magn);
magn = reshape(magn,[nx, ny, 1, nz])./max(magn(:));
mask = im2bw(magn, k*graythresh(magn));
mask = reshape(mask,sz);
mask = imclose (mask,strel('disk',4));
end