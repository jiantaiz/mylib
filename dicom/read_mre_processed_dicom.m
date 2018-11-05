function s=read_mre_processed_dicom(folder,savemat)
if nargin<1
    folder = 'proc';
end
if nargin<2
    savemat ='';
end

curlx = read_mre_dicom(folder,'*s4*.dcm',0);
curly = read_mre_dicom(folder,'*s5*.dcm',0);
curlz = read_mre_dicom(folder,'*s6*.dcm',0);
di = read_mre_dicom(folder,'*s8*.dcm',0);
osssnr =  read_mre_dicom(folder,'*s17*.dcm',0);
mag = read_mre_dicom(folder,'*s14*.dcm',0);
mask = create_mask(double(mag),0.8);
s.curlx = single(curlx);
s.curly = single(curly);
s.curlz = single(curlz);
s.mag = single(mag);
s.mask = single(mask);
s.di = single(di);
s.osssnr = single(osssnr);
if ~isempty(savemat)
    save(savemat, 's');
end