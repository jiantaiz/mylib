function write_mre_dicom(folder,cimgs,dinfos,se_desc,se_no)
% READ_MRE_DICOM reads in MRE data from GE dicom folder and 
% save in % cimgs(x,y,z,IQ_channel,dir,ph_offsets).
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


% if nargin <2
%     dcm_pattern ='*.dcm';
% end
% dcm = dir(fullfile(folder,dcm_pattern));
% dinfo = dicominfo(fullfile(folder,dcm(1).name));
    
if nargin<5
   se_no = [];
end
mkdir(folder);

parfor k=1:numel(dinfos);
if mod(k,10)==0
    disp(k)
end
    
%     dinfo(k) = dicominfo(fullfile(folder,dcm(k).name));
%     z = dinfo(k).InStackPositionNumber;% slice
%     t = dinfo(k).TemporalPositionIdentifier; %time offset step
%     d = dinfo(k).EchoNumber;% direction
%     c = dinfo(k).Private_0043_102f; %channel I(2),Q(3)
    dinfo = dinfos(k);  
   [~, fname,ext]=fileparts( dinfo.Filename);
   if isempty(ext)
       ext = '.dcm';
   end
   dcm_name = fullfile(folder,[fname ext]);
   
   
    z = dinfo.InStackPositionNumber;% slice
    if isfield(dinfo,'TemporalPositionIdentifier')
    t = dinfo.TemporalPositionIdentifier; %time offset step
    else
        t=1;
    end
    d = dinfo.EchoNumber;% direction
%     nd(d) = 1; % for removing all zero dimension later
    c = dinfo.Private_0043_102f; %channel I(2),Q(3)
    if c(1)>1
        c=c(1)-1;
    else
        c=1;
    end
    if c == 1
    im= real(cimgs(:,:,z,d,t));
    elseif c==2
    im= imag(cimgs(:,:,z,d,t));
    end
    im = int16(im);
    dinfo.SeriesDescription = se_desc;
    if ~isempty(se_no)
        dinfo.SeriesNumber = se_no;
    end
    dinfo.SeriesInstanceUID=[dinfo.SeriesInstanceUID,'0'];
    dicomwrite(im,dcm_name,dinfo,'CreateMode','copy','WritePrivate',true,'VR','explicit');
end
