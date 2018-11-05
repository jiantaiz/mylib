%% calculate menc, and transform matrix
freq_motion = 60; %vibration frequency
g=0.777;
pathstr='.';
cp_file = fullfile(pathstr,'..','cp','cornerPoints.all');
dcm =fullfile(pathstr,'i001.dcm');
info = dicominfo(dcm);
te = info.EchoTime;% get echo time

%calculate menc
[T,menc,menc_mtx]=calc_dual_menc_tmat(cp_file,g,freq_motion,te);

%% unwrap
load('cimgs.mat')
load('mask_sb.mat')
N_medfilt =1; % number of times to apply median filter.
high_vib = [0 1 0]; % x-,y-,z- high vibration flag.
z_filter = 0; 
brain_mask2 = [];
[phs_uw, phs, phs_est]  = unwrap_dual_menc(cimgs,T,N_medfilt,high_vib,z_filter);

%%
mask = mask >0;
mask2 = mask(:,:,:,ones(1,4),ones(1,3));
ov([phs_uw.*mask2, phs_est.*mask2, 1*(phs_uw.*mask2-phs_est.*mask2)],[],[-5 5])
%%
mse = mean((phs_uw(mask2) - phs_est(mask2)).^2); %mean squre error
msd = sqrt(mse)%mean squre deviation
save unwrapdata.mat