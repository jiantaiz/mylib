function mat  = nii2mat(fname)
nii = load_nii(fname);
mat = nii.img;
dim = ndims(mat);
mat = permute(mat, [2, 1, 3:dim]); %adjust image (ij) coord to physical (x y) coord
mat = flipud(mat);
