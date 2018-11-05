function PS=print_dicom_tags(dcm,fname)
if exist('gems-dicom-dict.txt','file')
    dinfo = dicominfo(dcm,'dictionary','gems-dicom-dict.txt');
else
    dinfo = dicominfo(dcm);
end

% if isfield(dinfo,'User25_User48')
%     dinfo.User25_User48 = sprintf('%g ',dinfo.User25_User48');
% end

addpath('C:\Users\m165355\OneDrive\Mayo\Matlab\online_scripts\bspm-master\dependencies')
PS=printstruct(dinfo,'maxarray',24,'write',fname);
