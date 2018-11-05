 f = rdir('*\*\*\cimgs.mat');
%   f = f(1:4);
 clear OMAT max_motion
 try
     parfor k=1:numel(f)
         k
         c = load(f(k).name);
         [vol_corr, omat]=ssh_motion_correction(abs(c.cimgs(:,:,:,[1 end])));
         OMAT(:,:,:,k) = omat;
         max_motion(k) = max(abs(omat(1:3,end,2)));
     end
     save motion_detect OMAT f max_motion
 catch err
     save motion_detect_inter OMAT f max_motion err
 end
