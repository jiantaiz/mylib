function [SensitivityRatios,VirtualCoils,info]=process_asset_cal(pfile)
% PROCESS_ASSET_CAL generate AssetCalibration.h5 from GE CAL scan pfile
%    [SensitivityRatios,VirtualCoils,info]=PROCESS_ASSET_CAL(pfile) ...
%
%    Example:
%    ... 
%
%    Subfunctions: 
%    See also: 

% AUTHOR    : Yi Sui
% DATE      : 08/23/2018
%%
[pathstr,name,ext] = fileparts(pfile); 
if (strcmpi(ext,'.h5'))
    h5file = pfile;
else % if not h5 file, then assume it is pfile, generate AssetCalibration.h5 using GERecon
    GERecon('Calibration.Process', pfile);
    h5file = fullfile(pathstr,'AssetCalibration.h5');
end
info = h5info(h5file);
SensitivityRatios = h5read(h5file,'/Regular/Data/SensitivityRatios');
SensitivityRatios = SensitivityRatios.real + 1i*SensitivityRatios.imag;
SensitivityRatios = permute(SensitivityRatios,[1 2 4 3]);
% SensitivityRatios = rot90(SensitivityRatios);
VirtualCoils = h5read(h5file,'/Regular/Data/VirtualCoils');
% VirtualCoils = rot90(VirtualCoils);

names={...
    'CalibrationMax'
    'CalibrationMean'
    'CalibrationType'
    'Channels'
    'CoilID'
    'Exam'
    'FirstSliceCornerPoints'
    'Landmark'
    'LastSliceCornerPoints'
    'PatientEntry'
    'PatientPosition'
    'RatioMean'
    'ScanCenter'
    'Slices'
    'XRes'
    'YRes'
    };
for k =1:numel(names);
info.(names{k}) = h5read(h5file,['/Regular/Info/',names{k}]);
end
