function [imgs, TI, dinfos]  = read_t1_dicom(folder,dcm_pattern,savemat)
% READ_MRE_DICOM reads in MRE data from GE dicom folder and 
% save in cimgs(x,y,z,IQ_channel,dir,ph_offsets).
%
%    [cimgs dinfo]  = READ_MRE_DICOM(folder,dcm_pattern) 
%         dcm_name_pattern = '*.dcm'(default) or 'IM*' etc.
%
%    Example:
%    [cimgs dinfo]  = READ_MRE_DICOM('.','*.dcm') 
%
%    See also: read_mre_dicom2(parfor version)

% AUTHOR    : Yi Sui
% DATE      : 01/06/2017


if nargin <2
    dcm_pattern ='*.dcm';
end

if nargin <3
    savemat = 0;
end

dcm = dir(fullfile(folder,dcm_pattern));

for k=numel(dcm):-1:1;
    dinfo=dicominfo(fullfile(folder,dcm(k).name));
    dinfos(k)=dinfo;
    ti(k) = dinfo.InversionTime;
    z(k) = dinfo.InStackPositionNumber;% slice
    im(:,:,k) = dicomread(fullfile(folder,dcm(k).name));   
end
TI = unique(ti);

for k=1:numel(dcm)
    t = find (TI == ti(k));
    imgs(:,:,z(k),t)=im(:,:,k); 
end
    
imgs = squeeze(imgs);
if savemat == 1
    seno = dinfo.SeriesNumber;
    sedesc = strrep(dinfo.SeriesDescription, ' ', '_');
    sedesc = sedesc;
    fname = sprintf('s%04d_%s',seno,sedesc);
    fname = matlab.lang.makeValidName(fname);
    imgs = single(imgs);
    save(fname,'imgs','dinfos','TI');
end
