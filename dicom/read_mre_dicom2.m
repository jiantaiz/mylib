function [cimgs dinfo]  = read_mre_dicom2(folder,dcm_name_pattern)
% [cimgs dinfo]  = READ_MRE_DICOM2(folder,dcm_name_pattern) 
%
% READ_MRE_DICOM2 reads in MRE data from GE dicom folder and 
% save in % cimgs(x,y,z,IQ_channel,dir,ph_offsets).
%
% dcm_name_pattern = '*.dcm'(default) or 'IM*' etc. 
%
% This is the parfor version of read_mre_dicom
% by Yi Sui
% 1/1/2017

if nargin <2
    dcm_name_pattern ='*.dcm';
end
dcm = dir(fullfile(folder,dcm_name_pattern));
% dinfo = dicominfo(fullfile(folder,dcm(1).name));

N = numel(dcm);
z=zeros(N,1);
t=zeros(N,1);
d=zeros(N,1);
c=zeros(N,1);
sz = size(dicomread(fullfile(folder,dcm(1).name)));
IM = zeros([sz N]);
parfor k=1:N;
    
    
    dinfo(k) = dicominfo(fullfile(folder,dcm(k).name));
    z(k) = dinfo(k).InStackPositionNumber;% slice
    t(k) = dinfo(k).TemporalPositionIdentifier; %time offset step
    d(k) = dinfo(k).EchoNumber;% direction
    tmp_c = dinfo(k).Private_0043_102f; %channel I(2),Q(3)
    c(k) = tmp_c(1);
    
    IM(:,:,k) = dicomread(fullfile(folder,dcm(k).name));
end
for k=1:N;
    cimgs(:,:,z(k),c(k)-1,d(k),t(k)) = IM(:,:,k);
end


if size(cimgs,4)>1
    cimgs = double(cimgs(:,:,:,1,:,:)) + 1i .* double(cimgs(:,:,:,2,:,:));
end
cimgs = squeeze(cimgs);