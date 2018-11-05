function mat2bin (filepath,filename,M)
%function mat2bin (filepath,filename,M)

fid = fopen(fullfile(filepath,filename),'w');
fwrite(fid,M,'float32');
fclose(fid);
