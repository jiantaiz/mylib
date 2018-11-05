function prot = getGEProtocolFromDicom(dcm,outputdir,outputname)
%    prot = GETGEPROTOCOLFROMDICOM(dcm,outputdir,outputname)
%         extract GE protocol info from DICOM and save to a protocol text file (LxProtocol)
%    INPUT:
%         dcm: DICOM file
%         outputdir: default './'
%         outputname: default 'LxProtocol'
%    OUTPUT:
%         prot: protocol content in text
%
%    Example:
%    GETGEPROTOCOLFROMDICOM('i001.dcm') 
%           creates LxProtocol in current directory
%    GETGEPROTOCOLFROMDICOM('i001.dcm','the/outputdir/dir') 
%           creates LxProtocol in the output dir
%    GETGEPROTOCOLFROMDICOM('i001.dcm','the/outputdir/dir','outputname') 
%           specifies output dir and output protocol name
%    prot = GETGEPROTOCOLFROMDICOM('i001.dcm')
%           returns the protocol content in text, but not create a file.
%
%    See also: 

% AUTHOR    : Yi Sui
% DATE      : 04/20/2018
%%

if nargin<2
    outputdir = './';
end
if nargin<3
    outputname = 'LxProtocol';
end
dinfo = dicominfo(dcm);
data = dinfo.Private_0025_101b; %Protocol data in gzip formate
tmp = [tempname(),'.gz'];
fid = fopen(tmp,'w');
fwrite(fid,data(5:end)); %write to a temporary .gz file
fclose(fid);
fname =gunzip(tmp);% gunzip file

if nargout>0
    fid = fopen(fname{1},'r');
    prot = fread(fid,'*char')';
    fclose(fid);
end
if nargout<1 || nargin>1
    movefile(fname{1},fullfile(outputdir,outputname)); %move to destination
    fprintf('Protocol has been saved in %s\n',fullfile(outputdir,outputname));
end
