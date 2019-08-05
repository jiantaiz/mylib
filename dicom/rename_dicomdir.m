function  ex = rename_dicomdir(d)
f = dir(d);
f = f(3:end); %remove . and ..
ex='';

dcm = f(~[f.isdir]);
fname='';

for k = 1:numel(dcm)
    if  isdicom(fullfile(d,dcm(k).name))
        fname = dcm(k).name;
        break;
    end
end


if ~isempty(dcm) && ~isempty(fname) &&  isdicom(fullfile(d,fname))
    dinfo = dicominfo(fullfile(d,fname));
    seno = dinfo.SeriesNumber;
    sedesc = strrep(dinfo.SeriesDescription, ' ', '_');
%     sedesc = matlab.lang.makeValidName(sedesc);
    dirname = sprintf('s%04d_%s',seno,sedesc);
    dirname = matlab.lang.makeValidName(dirname);
    
    exno = dinfo.StudyID;
    try
    exdesc = strrep(dinfo.StudyDescription,' ','_');
    catch
    exdesc = 'N_A';
    end
    exdate = dinfo.StudyDate;
    station = dinfo.StationName;
    ex = sprintf('e%s_%s_%s_%s',exno,exdesc,station,exdate);
    ex = matlab.lang.makeValidName(ex);

%     movefile(d, fullfile(d,'..',dirname));
    java.io.File(d).renameTo(java.io.File(fullfile(d,'..',dirname)));%much faster


end

dirs = f([f.isdir]);
if ~isempty(dirs)
    for k = 1:numel(dirs)
        try
        ex_dirname = rename_dicomdir(fullfile(d,dirs(k).name));
        catch
            continue;
        end
    end
    fullfile(d,'..',ex_dirname)
    ex_dirname
    java.io.File(d).renameTo(java.io.File(fullfile(d,'..',ex_dirname)));%much faster using java

end