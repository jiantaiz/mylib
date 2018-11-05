function ph = calc_small_ph(cimgs,sl,phase_shift)
a=1;b=3;c=pi/2;
a=0;b=0;c=0;
dir = 1;
ydata1 =squeeze( cimgs(:,:,sl,:,dir));% positive meg
ydata2 =squeeze( cimgs(:,:,sl,:,dir+1)); % negtive meg
ydata12(:,:,:,1) = rm_linear_ph(ydata1.*ydata2,a,b,c);

dir = 3;
ydata1 =squeeze( cimgs(:,:,sl,:,dir));% positive meg
ydata2 =squeeze( cimgs(:,:,sl,:,dir+1)); % negtive meg
ydata12(:,:,:,2) = rm_linear_ph(ydata1.*ydata2,a,b,c);

dir = 5;
ydata1 =squeeze( cimgs(:,:,sl,:,dir));% positive meg
ydata2 =squeeze( cimgs(:,:,sl,:,dir+1)); % negtive meg
ydata12(:,:,:,3) = rm_linear_ph(ydata1.*ydata2,a,b,c);
%%
ph = angle(ydata12(:,:,:,[1 2 3]).*conj(ydata12(:,:,[5:8,1:4],[2 3 1])));
% ph = angle(ydata12(:,:,:,:)); 
ph = unwrap(ph,[],3);
ph_fft = ifft(ph,[],3);
ph_fft(:,:,1,:)=0;
disp 'small phase dc=0';
sz = size(ph);
% phase_shift = 0.1999;

for k=1:sz(3)/2;
    ph_fft(:,:,k+1,:) = ph_fft(:,:,k+1,:).*exp(1i.*k*phase_shift);
    ph_fft(:,:,end-k+1,:) = ph_fft(:,:,end-k+1,:).*exp(-1i.*k*phase_shift);
end



ph = real(fft(ph_fft,[],3));
ph =unwrap(ph,[],3);
ph = mean(ph,4);
% imdisp(ph (:,:,:),'size',[nan 8])
%%
% imdisp(mean(ph,4),'size',[nan 8]);