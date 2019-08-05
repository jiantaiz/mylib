function  ex = rename_pfiledir(d)
f = dir(d);
f = f(3:end); %remove . and ..
ex='';

pfile = dir(fullfile(d,'*.7'));

if ~isempty(pfile)
    pfileHandle = GERecon('Pfile.Load', fullfile(d,pfile(1).name),'No-Anonymize');
    header = GERecon('Pfile.Header',pfileHandle);
    seno =  header.SeriesData.se_no;
    sedesc = strrep(nonzeros(header.SeriesData.se_desc)', ' ', '_');

%     sedesc = matlab.lang.makeValidName(sedesc);
    dirname = sprintf('s%04d_%s',seno,sedesc);
    
    dirname = matlab.lang.makeValidName(dirname);
    
    
    exno = header.ExamData.ex_no;
%     try
%     exdesc = strrep(header.ExamData.ex_desc,' ','_');
%     catch
%     exdesc = 'N_A';
%     end
    patname = header.ExamData.patnameff;
    patname = patname(patname>0);
    exdate = datetime(1970,1,1,0,0,0,'Format','MMddyyyy') +header.ExamData.ex_datetime/60/60/24;
    station =header.ExamData.ex_suid;
    ex = sprintf('e%d_%s_%s_%s',exno,patname,station,char(exdate));
    ex = matlab.lang.makeValidName(ex);

%     movefile(d, fullfile(d,'..',dirname));
    java.io.File(d).renameTo(java.io.File(fullfile(d,'..',dirname)));%much faster


end

dirs = f([f.isdir]);
if ~isempty(dirs)
    for k = 1:numel(dirs)
        ex_dirname = rename_pfiledir(fullfile(d,dirs(k).name));
    end
    fullfile(d,'..',ex_dirname)
%     ex_dirname
    java.io.File(d).renameTo(java.io.File(fullfile(d,'..',ex_dirname)));%much faster using java

end