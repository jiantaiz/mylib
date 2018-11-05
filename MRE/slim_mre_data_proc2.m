for se = [4 5 6 7 8]
% se = 4;
data(se)  = load (sprintf('s%d\\results2.mat',se));
end
%% check phase unwrap results
sl=6;
se = 8;
pdimg = data(se).pdimg;
ph_uw = data(se).ph_uw;
imdisp(cat(2,angle(pdimg(:,:,sl,:)),ph_uw(:,:,sl,:)),[-pi pi]);
set(gcf,'name',sprintf('se%d',se))


%% compare slim results

se1 = 4;
se2 = 5;
% sl=10
ph_uw_fft1 = data(se1).ph_uw_fft;
ph_uw_fft2 = data(se2).ph_uw_fft;
f= 2:4
sl = 10
im = cat(2,ph_uw_fft1(:,:,sl,f),ph_uw_fft2(:,:,sl,f).*exp(1i*0));
imdisp(real(im),[-pi pi])
%% check magnitude images
sl= 10:20
pdimg1 = data(se1).pdimg;
pdimg2 = data(se2).pdimg;
im  = cat(2,pdimg1(:,:,sl,1),pdimg2(:,:,sl,1));
imdisp(sqrt(im))

%% interplate time series
N=3;
clear img2 img1
for n = 0:N-1
    phoff = 2*pi*n/N;
    img1(:,:,:,:,n+1) = ph_uw_fft1.*exp(1i.*phoff);
    img2(:,:,:,:,n+1) = ph_uw_fft2.*exp(1i.*(phoff));

end
sl =10
% bulk motion filter
img1 = real(img1);
img2 = real(img2);
mask = abs(pdimg1).^(1/2) ;
mask = mask > 700;

% ph_bmf1 = remove_bulk_motion(img1(:,:,sl,:,:),repmat(mask(:,:,sl,:),[1 1 1 1,32]),'average');
% ph_bmf2 = remove_bulk_motion(img2(:,:,sl,:,:),repmat(mask(:,:,sl,:),[1 1 1 1,32]),'average');


% ph_bmf1 = remove_bulk_motion(img1(:,:,sl,:,:),repmat(mask(:,:,sl,:),[1 1 1 1,N]),'bhp',4,5);
% ph_bmf1 = remove_bulk_motion(ph_bmf1,repmat(mask(:,:,sl,:),[1 1 1 1,N]),'blpf',40,5);
% ph_bmf2 = remove_bulk_motion(img2(:,:,sl,:,:),repmat(mask(:,:,sl,:),[1 1 1 1,N]),'bhp',4,5);
% ph_bmf2 = remove_bulk_motion(ph_bmf2,repmat(mask(:,:,sl,:),[1 1 1 1,N]),'blpf',40,5);

ph_bmf1 = remove_bulk_motion(img1(:,:,sl,:,:),repmat(mask(:,:,sl,:),[1 1 1 1,N]),'poly',2);

ph_bmf2 = remove_bulk_motion(img2(:,:,sl,:,:),repmat(mask(:,:,sl,:),[1 1 1 1,N]),'poly',2);

ph_bmf1_fft = fft(ph_bmf1,[],5);
ph_bmf2_fft = fft(ph_bmf2,[],5);
%

N=32;
 clear img2 img1 ph_bmf1 ph_bmf2;
 
for n = 0:N-1
    phoff = 2*pi*n/N;
    ph_bmf1(:,:,:,:,n+1) = ph_bmf1_fft(:,:,:,:,2).*exp(1i.*phoff);
    ph_bmf2(:,:,:,:,n+1) = ph_bmf2_fft(:,:,:,:,2).*exp(1i.*phoff);

end

%%

% imdisp(real(cat(2,ph_bmf1(:,:,1,2,:),ph_bmf1(:,:,1,3,:),ph_bmf1(:,:,1,4,:))))
t = 1:N;
sl = 11

im1 =real(cat(2,ph_bmf1(:,:,sl,2,t),ph_bmf1(:,:,sl,3,t),ph_bmf1(:,:,sl,4,t)));
im2 =real(cat(2,ph_bmf2(:,:,sl,2,t),ph_bmf2(:,:,sl,3,t),ph_bmf2(:,:,sl,4,t)));
imdisp(cat(1,im1,im2),[-pi/8 pi/8],'size',[1 1]);
%%
f=2;
im1 =real(cat(2,ph_bmf1_fft(:,:,1,2,f),ph_bmf1_fft(:,:,1,3,f),ph_bmf1_fft(:,:,1,4,f)));
im2 =real(cat(2,ph_bmf2_fft(:,:,1,2,f),exp(1i.*-pi/4).*ph_bmf2_fft(:,:,1,3,f),exp(1i.*0)*ph_bmf2_fft(:,:,1,4,f)));
imdisp(cat(5,im1,im2),[-pi/8 pi/8])

%%
% ph_bmf1 = remove_bulk_motion(ph_uw_fft1(:,:,sl,:,:),repmat(mask(:,:,sl,:),[1 1 1 1,32]),'average');


ph_bmf1 = get_shere_wave(ph_uw_fft1,mask,8);
ph_bmf2 = get_shere_wave(ph_uw_fft2,mask,8);



%%
tic
for se = [4 5 6 7 8]
    % se = 4;
    ph_uw_fft = data(se).ph_uw_fft;
    pdimg = data(se).pdimg;
    mask = sqrt(abs(pdimg)) > 200 & ~isnan(ph_uw_fft);
    ph_bmf= get_shear_wave(ph_uw_fft(:,:,:,2:4),mask,8);
    data(se).ph_bmf = ph_bmf;
end

toc

%%
sl = 4:3:20;
f=1;
se = 8
im = data(se).ph_bmf(:,:,sl,f,:);
im = mosaic(im,3,2);
%  im = ph_bmf_tmp(:,:,sl,f,:);
% im = mask(:,:,sl,f,:);
%  im = ph_uw_fft(:,:,sl,:);
imdisp((im),'size',[1 1])

%%
sl = 4:3:20;
f=2;
im = data(se).ph_uw_fft(:,:,sl,f);
im = mosaic(im,3,2);
%  im = ph_bmf_tmp(:,:,sl,f,:);
% im = mask(:,:,sl,f,:);
%  im = ph_uw_fft(:,:,sl,:);
imdisp(imag(im))
%%
tic
for se = [4 5 6 7 8]
    ph_uw = data(se).ph_uw;
    pdimg = data(se).pdimg;
    mask = sqrt(abs(pdimg)) > 400 & ~isnan(ph_uw);
    mask = mask(:,:,:,1);
    ph_uw (isnan(ph_uw)) = 0;
    sz = size(ph_uw);
    fd = fopen(sprintf('se%d_ph_uw_%dx_%dy_%dz_%dt.phs',se,sz(1:4)),'w');
    fwrite(fd,ph_uw,'float32');
    fclose(fd);
    
    fd = fopen(sprintf('se%d_mask_%dx_%dy_%dz.r4',se,sz(1:3)),'w');
    fwrite(fd,mask,'float32');
    fclose(fd);
end
%%
tic
for se = [4 5 6 7 8]
    ph_uw = data(se).ph_uw;
    sz = size(ph_uw);
    fd = fopen(sprintf('se%d_ph_uw_%dx_%dy_%dz_%dt_BMF1.0.phs',se,sz(1:4)),'r');
    data(se).ph_bmf_DL = reshape(fread(fd,'float32'), sz );
    
    fclose(fd);
    data(se).ph_bmf_fft = ifft(data(se).ph_bmf_DL , [] , 4);
end
toc
%%
se=4
ph_bmf1_fft = ifft(data(se).ph_bmf,[],5);
ph_bmf2_fft = ifft(data(se).ph_bmf_DL,[],4);
ph_bmf2_fft = reshape(ph_bmf2_fft(:,:,:,2:4),[256,256,20,3,1]);
%

N=16;
 clear img2 img1 ph_bmf1 ph_bmf2;
 
for n = 0:N-1
    phoff = 2*pi*n/N;
    ph_bmf1(:,:,:,:,n+1) = ph_bmf1_fft(:,:,:,:,2).*exp(1i.*phoff);
    ph_bmf2(:,:,:,:,n+1) = ph_bmf2_fft(:,:,:,:,1).*exp(1i.*phoff);

end

%
%%
% imdisp(real(cat(2,ph_bmf1(:,:,1,2,:),ph_bmf1(:,:,1,3,:),ph_bmf1(:,:,1,4,:))))
t = 1:N;
sl = 15

im1 =real(cat(2,ph_bmf1(:,:,sl,1,t),ph_bmf1(:,:,sl,2,t),ph_bmf1(:,:,sl,3,t)));
im2 =real(cat(2,ph_bmf2(:,:,sl,1,t),ph_bmf2(:,:,sl,2,t),ph_bmf2(:,:,sl,3,t)));
imdisp(cat(1,im1,im2),[-pi/8 pi/8],'size',[1 1]);
%%
for se = [4 5 6 7 8]
    
    ph_bmf2_fft = ifft(data(se).ph_bmf_DL,[],4);
    ph_bmf2_fft = ph_bmf2_fft(:,:,:,2:4);
    %
    
    N=8;
    clear img2 img1 ph_bmf1 ph_bmf2 ph_bmf_t;
    
    for n = 0:N-1
        phoff = 2*pi*n/N;
        ph_bmf_t(:,:,:,:,n+1) = real(2*ph_bmf2_fft(:,:,:,:).*exp(-1i.*phoff));
        
    end
    data(se).ph_bmf_t = ph_bmf_t;
end


%%
sl = 3
im1 = mosaic(data(4).ph_bmf_t(:,:,sl,:,end:-1:1),4);
im2 = mosaic(data(8).ph_bmf_t(:,:,sl,:,end:-1:1),4);
imdisp(cat(1,im1,im2),[-pi/2 pi/2],'size',[1 1],'map','awave');






