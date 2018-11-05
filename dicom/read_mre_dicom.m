function [cimgs, dinfos]  = read_mre_dicom(folder,dcm_pattern,savemat)
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

% 4/4/2018 : cimgs(x,y,z,c,dir,t) -> cimgs(:,:,:,z,c,t,dir)


if nargin <2
    dcm_pattern ='*.dcm';
end

if nargin <3
    savemat = 0;
end

dcm = dir(fullfile(folder,dcm_pattern));

dinfo = dicominfo(fullfile(folder,dcm(1).name));

z0 = double(dinfo.SliceLocation);
thick = double(dinfo.SliceThickness);
numofslices = double(dinfo.Private_0021_104f);
zlocs = z0+[0:numofslices-1].*thick;

for k=numel(dcm):-1:1;
    
%     dinfo(k) = dicominfo(fullfile(folder,dcm(k).name));
%     z = dinfo(k).InStackPositionNumber;% slice
%     t = dinfo(k).TemporalPositionIdentifier; %time offset step
%     d = dinfo(k).EchoNumber;% direction
%     c = dinfo(k).Private_0043_102f; %channel I(2),Q(3)
    dinfo = dicominfo(fullfile(folder,dcm(k).name));
    
    if isfield(dinfo,'InStackPositionNumber')
       z = dinfo.InStackPositionNumber;% slice
    else
        zloc = double(dinfo.SliceLocation);
        z= find(abs(zlocs-zloc)<0.0001);
    end
    
    if isfield(dinfo,'TemporalPositionIdentifier')
        t = dinfo.TemporalPositionIdentifier; %time offset step
    else
        t=1;
    end
    d = dinfo.EchoNumber;% direction
    nd(d) = 1; % for removing all zero dimension later
    c = dinfo.Private_0043_102f; %channel I(2),Q(3)
    if c(1)>1
        c=c(1)-1;
    else
        c=1;
    end
    
    im=dicomread(fullfile(folder,dcm(k).name));    
    cimgs(:,:,:,z,c,t,d) = im ;
    dinfos(k)=dinfo;
end
cimgs=cimgs(:,:,:,:,:,:,nd==1);
if size(cimgs,5)>1
    cimgs = double(cimgs(:,:,:,:,1,:,:)) + 1i .* double(cimgs(:,:,:,:,2,:,:));
end
cimgs = squeeze(cimgs);
if savemat == 1
    seno = dinfo.SeriesNumber;
    sedesc = strrep(dinfo.SeriesDescription, ' ', '_');
    sedesc = sedesc;
    fname = sprintf('s%04d_%s',seno,sedesc);
    fname = matlab.lang.makeValidName(fname);
    cimgs = single(cimgs);
    save(fname,'cimgs','dinfos');
end
