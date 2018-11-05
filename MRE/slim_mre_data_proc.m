cimgs = s16newslimsag.cimgs;
sl = 3;
phase_shift = 0;
dir=1;
timeoffset=1:8;
a=0;
b=0;
c=0;
ph = calc_small_ph_tetra(cimgs(:,:,sl,:,:),phase_shift,a,b,c);

imdisp(ph(:,:,:,:))

ph = angle(cimgs(:,:,sl,:,1).*conj(cimgs(:,:,sl,:,2).*exp(1i.*-0)));
ph = unwrap(ph,[],4);
imdisp(ph(:,:,:,:,:))

ph = angle(cimgs(:,:,:,:,1).*(cimgs(:,:,:,:,2).*exp(1i.*-0)));

imdisp(ph(:,:,sl,:,:))
ph = angle(cat(2, cimgs(:,:,:,:,1),(cimgs(:,:,:,:,2))));
imdisp(ph(:,:,sl,:,:))
%%
% cimgs = s13newslim.cimgs;
% cimgs = s16newslimsag.cimgs;
% cimgs = s6oldslim.cimgs;
for s=[7]
    clear cimgs
    load (sprintf('s%d\\cimgs.mat',s));
    
    % cimgs = fftInterpolate(cimgs,2).*imresize(abs(cimgs)>0.1,2,'nearest');
    clear ph_uw ph_qual
    dir = 1;
    rotang = pi;
    time_offset = 1:size(cimgs,4);
    pdimg = cimgs(:,:,:,time_offset,dir).* conj(cimgs(:,:,:,time_offset,dir+1)) .*exp(1i.*rotang);
    mask = sqrt(abs(pdimg)) > 400;
    % mask=[];
    %
    [nx, ny, nz,nt] = size(pdimg);
    for sl = 6
        sl
        parfor t=1:nt;
            [ph_uw(:,:,sl,t), ph_qual(:,:,sl,t)]=unwrap2D(pdimg(:,:,sl,t),mask(:,:,sl,t));
        end
    end
    ph_fft = ifft(angle(pdimg),[],4);

    ph_uw_fft = ifft(ph_uw,[],4);
    
%     save( sprintf('s%d\\results.mat',s), 'ph_uw_fft' , 'ph_uw', 'pdimg','mask','ph_qual');
end
%%
s6oldslim.ph_uw_fft = ph_uw_fft;
s6oldslim.ph_uw = ph_uw;
save s6oldslim

%%

s16newslimsag.ph_uw = ph_uw;
s16newslimsag.ph_uw_fft = ph_uw_fft;
%%
im = cat(2,ph_uw_fft(:,:,sl,:),ph_fft(:,:,sl,:)); 
imdisp(abs(im))
%%

imdisp(cat(2,angle(pdimg(:,:,sl,:)),ph_uw(:,:,sl,:)));
imdisp(cat(2,angle(pdimg(:,:,sl,:)),ph_uw(:,:,sl,:), unwrap(ph_uw(:,:,sl,:),[pi*1.2],4) ) );



imdisp(real(cat(2,s14oldslim.ph_uw_fft(:,:,sl,:),s13newslim.ph_uw_fft(:,:,sl,:)   )))



%%
N=32;
% clear img2
for n  = 0:N-1
    phoff = 2*pi*n/N;
    img1(:,:,:,:,n+1) = ph_uw_fft.*exp(1i.*phoff);
    
end
%%
imdisp(real(cat(2,img2(:,:,sl,2,:),img2(:,:,sl,3,:),img2(:,:,sl,4,:))))
imdisp(real(cat(2,img1(:,:,sl,2,:),img2(:,:,sl,2,:))))



%%
mask = abs(pdimg).^(1/2) ;
mask = mask > 700;
imdisp(mask(:,:,1:20))
%%
% ph_bmf = remove_bulk_motion(ph_uw,mask,'gaussian',14);
ph_bmf = remove_bulk_motion(ph_uw,mask,'bhp',3,7);
%%
im = cat(2,angle(pdimg(:,:,sl,:).*exp(1i.*pi*1.1)),ph_uw(:,:,sl,:),ph_bmf(:,:,sl,:));
imdisp(real(im));

%%
ph_bmf_fft = ifft(real(ph_bmf),[],4);
imdisp(real(ph_bmf_fft(:,:,5,:)))

%%
% ph_bmf = remove_bulk_motion(img2(:,:,sl,:,:),repmat(mask(:,:,sl,:),[1 1 1 1,32]),'average');
ph_bmf = remove_bulk_motion(img2(:,:,sl,:,:),repmat(mask(:,:,sl,:),[1 1 1 1,32]),'bhp',3,5);
ph_bmf = remove_bulk_motion(ph_bmf,repmat(mask(:,:,sl,:),[1 1 1 1,32]),'blpf',25,5);

imdisp(real(cat(2,ph_bmf(:,:,1,2,:),ph_bmf(:,:,1,3,:),ph_bmf(:,:,1,4,:))))



%%
im= cat(2,ph_bmf(:,:,1,2,:),ph_bmf(:,:,1,3,:),ph_bmf(:,:,1,4,:));
im = max(im,[],5);
imdisp(im,im);



%% remove 2pi phase jump 
for se=[4 5 6 7 8 ]
    load (sprintf('s%d\\results.mat',se));
    [nx, ny, nz,nt] = size(pdimg);
    mask = sqrt(abs(pdimg)) > 600;
    clear med1 med2
    threshold = pi;
    for sl = 1:nz
        for t = 1:nt
            s1 = ph_uw (:,:,sl,t);
            s1 = s1(mask(:,:,sl,t));
            med1(t) = median(s1,'omitnan');
            s2 = angle(pdimg (:,:,sl,t));
            s2 = s2(mask(:,:,sl,t));
            med2(t) = median(s2,'omitnan');
        end
        idx = find((med1-med2) > threshold);
        ph_uw(:,:,sl,idx) = ph_uw(:,:,sl,idx) - 2*pi;
        
        idx = find((med1-med2) < -threshold);
        ph_uw(:,:,sl,idx) = ph_uw(:,:,sl,idx) + 2*pi;
    end
    ph_uw_fft = ifft(ph_uw,[],4);
    
    save (sprintf('s%d\\results2.mat',se), 'pdimg', 'mask', 'ph_uw', 'ph_qual', 'ph_uw_fft');
end

%%




imdisp(cat(2,angle(pdimg(:,:,sl,:)),ph_uw(:,:,sl,:)));


ph_uw_fft2 = ifft(ph_uw,[],4);


imdisp(real(cat(2,ph_uw_fft2(:,:,sl,:),ph_uw_fft(:,:,sl,:))))