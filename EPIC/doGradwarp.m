function gradwarpImage = doGradwarp(rawImage,pfileFullPath,gwcoefs,gradwarpMethod)
% DOGRADWARP ...
%    gradwarpImage = DOGRADWARP(rawImage,pfileFullPath,gwcoefs) ...

%    Example:
%    ...
%
%    Subfunctions:
%    See also:

% AUTHOR    : Yi Sui
% DATE      : 03/18/2019
%%

if nargin<3 || isempty(gwcoefs)
    gwcoefs = get_gw_coils_mr55();
end

if nargin <4 
    gradwarpMethod = '2DGradWarp';
end
[pathstr,name,ext] = fileparts(pfileFullPath);
switch lower(ext)
    case '.7'
        pfileHandle = GERecon('Pfile.Load', pfileFullPath);
        header = GERecon('Pfile.Header', pfileHandle);
        nslices = pfileHandle.slices;
    case '.h5'
            archive = GERecon('Archive.Load', pfileFullPath); 
            nslices = archive.Slices;

end







%    NOTE: GERecon('Gradwarp') takes first dim (column) as X, second dim(row) as Y. So
%    the rawImage needs to be rotated to comply with this convention.
rawImage = rot90(rawImage);

[nx, ny, nz, nv] = size(rawImage);
gradwarpImage = zeros(size(rawImage));

switch lower(gradwarpMethod)
    case {'2dgradwarp', '2d'}
        disp('2D Gradwarp')
        for slice = 1:nslices
            if (strcmpi(ext,'.7'))
                corners = GERecon('Pfile.Corners', slice);
            elseif (strcmpi(ext,'.h5'))
                sliceInfo = GERecon('Archive.Info', archive,  slice);
                corners = sliceInfo.Corners;
            end
            
            parfor vol = 1:nv
                
                oneSlice = rawImage(:,:,slice,vol);
                re = GERecon('Gradwarp', real(oneSlice), corners, 'SphericalHarmonicCoefficients', gwcoefs);
                im = GERecon('Gradwarp', imag(oneSlice), corners, 'SphericalHarmonicCoefficients', gwcoefs);
                gradwarpImage(:,:,slice,vol) = re + 1i*im;
            end
        end
    case {'3dgradwarp', '3d'}
        disp('3D Gradwarp')
        if (strcmpi(ext,'.7'))
            corners(1) = GERecon('Pfile.Corners', 1);
            corners(2) = GERecon('Pfile.Corners', nz);
        elseif (strcmpi(ext,'.h5'))
            sliceInfo = GERecon('Archive.Info', archive,  1);
            corners(1) = sliceInfo.Corners;
            sliceInfo = GERecon('Archive.Info', archive,  nz);
            corners(2) = sliceInfo.Corners;
        end

        for vol = 1:nv
            oneVol = rawImage(:,:,:,vol);
            re = GERecon('Gradwarp', real(oneVol), corners, 'SphericalHarmonicCoefficients', gwcoefs);
            im = GERecon('Gradwarp', imag(oneVol), corners, 'SphericalHarmonicCoefficients', gwcoefs);
            gradwarpImage(:,:,:,vol) = re + 1i*im;
        end
    otherwise
        error('gradwarpMethod = 2DGradwarp or 3DGradwarp');
end

gradwarpImage = rot90(gradwarpImage,-1);