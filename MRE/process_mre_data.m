function [ph, ph_uw, mag, mask, dimg] = process_mre_data(matfile,save_to_file)
if nargin<2
    save_to_file=0;
end
load(matfile,'cimgs');
dimg = cimgs(:,:,:,:,1:2:end).*conj(cimgs(:,:,:,:,2:2:end));
ph = angle(dimg);
ph_uw = unwrap3D_ssh(ph,1);
mag = mean(abs(cimgs(:,:,:,:)),4);
mask = create_mask(mag,0.5);
if save_to_file
    fprintf('saving results to %s\n',matfile);
    save (matfile, 'dimg','ph','ph_uw','mag','mask','-append');
end