%2D unwrap test


% pdimg = cimgs(:,:,:,:,1:2:end).*conj(cimgs(:,:,:,:,2:2:end));

sl =5;
t=6
mask2 = M(:,:,sl);
mask2(64,64) = 23;
pdimg = exp(1i.*double(phs_z(:,:,sl,t))).*mask2;

sl = 3;
t=1;
% mask = sqrt(abs(pdimg)) >100;

% mask = M(:,:,sl);
mask = M;
mask2 = M;

mask2(64,64,:) = 3;
pdimg = bsxfun(@times, exp(1i.*double(phs_z)), mask2);
[nx ny nz nt nd] = size(pdimg)

pickstart=0;

for sl = 1:nz
    sl
    parfor t = 1:nt
        t
        [ph_uw(:,:,sl,t), ph_qual(:,:,sl,t)]=unwrap2D(pdimg(:,:,sl,t),mask(:,:,sl),pickstart);
    end
end
%%
imdisp(cat(2,angle(pdimg(:,:,6,:)),ph_uw(:,:,6,:)))

%%
for se = [5 6 7 8]
% se = 4;
if ~exist(sprintf('s%d\\cimgs_org.mat',se))
    copyfile(sprintf('s%d\\cimgs.mat',se),sprintf('s%d\\cimgs_org.mat',se))
end
 load (sprintf('s%d\\cimgs.mat',se));
 
 cimgs = cimgs(1:2:end,1:2:end,:,:,:);
 
 save (sprintf('s%d\\cimgs.mat',se),'cimgs');
end