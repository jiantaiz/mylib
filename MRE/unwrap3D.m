function M = unwrap3D(mat, DO_MEAN_SHIFT,varargin)
% UNWRAP3D is a wrapper of Josh's unix program GC_UNWRAP (graph cut unwrap)
%    M = UNWRAP3D(mat, DO_MEAN_SHIFT)
%    See also: 

% AUTHOR    : Yi Sui
% DATE      : 05/16/2017

if ispc
    warning('Require linux os, using unwrap3D_ssh instead.');
    M=unwrap3D_ssh(mat,DO_MEAN_SHIFT,varargin{:});    
    return;
end

if nargin<2
    DO_MEAN_SHIFT = 0;
end
tmp_file_in = tempname();
tmp_file_out = tempname();
mat2bin('',tmp_file_in,mat);

% GC_UNWRAP IN_FILENAME OUT_FILENAME Nx Ny Nz Nt DO_MEAN_SHIFT

[Nx, Ny, Nz, Nt]=size(mat);

cmd = sprintf('GC_UNWRAP %s %s %d %d %d %d %d',tmp_file_in,tmp_file_out,Nx, Ny, Nz, Nt, DO_MEAN_SHIFT );
system(cmd);

dims=size(mat);
fid = fopen(tmp_file_out,'r');
M = fread(fid,'single');
M = single(reshape(M,dims));
fclose(fid);
