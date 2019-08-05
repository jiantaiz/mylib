function reorgnize_folder(d)
% reorgnize the archive folder by Kevin's son-of-recon
fprintf('Orgnizing %s\n',d);
folder = dir(d);
folder = folder([folder.isdir]);

for k=3:numel(folder)
    reorgnize_folder(fullfile(d,folder(k).name));
end

dcm = dir(fullfile(d,'*.MR'));
if isempty(dcm)
    dcm = dir(fullfile(d,'*.dcm'));
end
if ~isempty(dcm)
    mkdir(fullfile(d,'dcm'));
    
    if exist('gems-dicom-dict.txt','file')
        dinfo = dicominfo(fullfile(d,dcm(1).name),'dictionary','gems-dicom-dict.txt');
    else
        dinfo = dicominfo(fullfile(d,dcm(1).name));
    end
    seno = dinfo.SeriesNumber;
    sedesc = strrep(dinfo.SeriesDescription, ' ', '_');
    se_name = sprintf('s%04d_%s',seno,sedesc);
    se_name = matlab.lang.makeValidName(se_name);
    try
      getGEProtocolFromDicom(fullfile(d,dcm(1).name),d);
    catch
      disp('Error when get protocol from dicom');  
    end
    PS=printstruct(dinfo,'maxarray',24,'write',fullfile(d,'dicom_tags.txt'));
    
    if (isunix) %similar performance as java
        cmd = sprintf('mv %s %s',fullfile(d,'*.MR'), fullfile(d,'dcm',filesep));
        system(cmd);
        cmd = sprintf('rename .MR .MR.dcm %s',fullfile(d,'dcm','*'));
        system(cmd);        
    else
        for k=1:numel(dcm)
            java.io.File(fullfile(d,dcm(k).name)).renameTo(java.io.File(fullfile(d,'dcm',[dcm(k).name,'.dcm'])));%much faster
        end
        
        
    end
    
    
else
    se_name = [];
end

dcm = dir(fullfile(d,'epidcm*.proc'));
if ~isempty(dcm)
    mkdir(fullfile(d,'proc'));
    
    if isunix
        cmd = sprintf('mv %s %s',fullfile(d,'*.proc'), fullfile(d,'proc',filesep));
        system(cmd);
        cmd = sprintf('rename .proc .proc.dcm %s',fullfile(d,'proc','*'));
        system(cmd);  
    else
        for k=1:numel(dcm)
            java.io.File(fullfile(d,dcm(k).name)).renameTo(java.io.File(fullfile(d,'proc',[dcm(k).name,'.dcm'])));%much faster
        end
    end
end

if ~isempty(se_name)
    if isunix
        cmd = sprintf('mv %s %s',d, fullfile(d,'..',se_name));
        system(cmd);
    else
        java.io.File(d).renameTo(java.io.File(fullfile(d,'..',se_name))); %much faster
    end
end

