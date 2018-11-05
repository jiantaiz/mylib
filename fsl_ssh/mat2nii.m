function fname  = mat2nii(mat,voxsize,fname)
if nargin<2
    voxsize = [1 1 1];
end
if nargin<3
    fname= [tempname,'.nii']; %create a tmp file name
end
dim = ndims(mat);
mat = flipud(mat);
mat = permute(mat, [2, 1, 3:dim]); %adjust image (ij) coord to physical (x y) coord
nii=make_nii(mat,voxsize);
save_nii(nii,fname);