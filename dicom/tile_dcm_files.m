function IM = tile_dcm_files(dicom_dir, clim, n)
% read in a dicom directory and tile the images in a row.
% IM: one tiled image
% clim: [amin amax] display window level.
% n: size of colormap default 256
if nargin <2
    clim =[];
end
if nargin <3
    n=256;
end
IM = read_mre_dicom(dicom_dir);
IM = IM(:,:,:,2);
IM  = tile(IM,3,2);
if ~isempty(clim)
    IM = mat2gray(IM, clim)*n;
end