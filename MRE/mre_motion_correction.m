
% clear
% load('cimgs.mat')
pdimg = cimgs(:,:,:,:,1:2:end) .* conj(cimgs(:,:,:,:,2:2:end));
pdimg = sqrt(abs(pdimg)).*exp(1i*angle(pdimg));

[pdimg_corr, omat]=ssh_motion_correction(pdimg(:,:,:,:),'spline');
pdimg_corr = reshape(pdimg_corr,size(pdimg));
pdimg_corr = flipud(pdimg_corr);
%%

% ov(abs(pdimg))
% ov(angle(pdimg))
TH=2;
FOV=240;
coord=0;
f=80;
FOVz_int = TH*size(pdimg,3);

[magn, divergence, curlx, curly, curlz] = AA_Curl(pdimg(:,:,:,:,1), pdimg(:,:,:,:,2), pdimg(:,:,:,:,3),TH,FOV , coord);
[magn, divergence, curlxc, curlyc, curlzc] = AA_Curl(pdimg_corr(:,:,:,:,1), pdimg_corr(:,:,:,:,2), pdimg_corr(:,:,:,:,3),TH,FOV , coord);

[di_recon cdi_recon] = di_JDT_AM([FOV/1000 FOV/1000 FOVz_int/1000], f,[7 7 7], 3, ones(size(curlx)),2, curlx, curly, curlz);

[di_recon_c cdi_recon_c] = di_JDT_AM([FOV/1000 FOV/1000 FOVz_int/1000], f,[7 7 7], 3, ones(size(curlxc)),2, curlxc, curlyc, curlzc);

lfe_recon = lfe([FOV/1000 FOV/1000 FOVz_int/1000], f, 3, curlx, curly, curlz);
lfe_recon_c = lfe([FOV/1000 FOV/1000 FOVz_int/1000], f, 3, curlxc, curlyc, curlzc);

save pdimg_mot_corr magn pdimg_corr omat di_recon cdi_recon di_recon_c cdi_recon_c curlx curly curlz curlxc curlyc curlzc