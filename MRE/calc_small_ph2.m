function ph = calc_small_ph2(cimgs,sl,phase_shift,dir,timeoffset)
if nargin <4
 dir = 5;
end
if nargin <5
    timeoffset = 1:size(cimgs,4);
end
ydata1 =squeeze( cimgs(:,:,sl,timeoffset,dir));% positive meg
ydata2 =squeeze( cimgs(:,:,sl,timeoffset,dir+1)); % negtive meg
ydata12 = rm_linear_ph(ydata1.*ydata2,1,3,pi/2);
%%
ph = angle(ydata12);
% ph = angle(ydata12(:,:,:,:)); 
ph = unwrap(ph,[],3);
ph_fft = ifft(ph,[],3);
ph_fft(:,:,1,:)=0;
disp 'set small phase dc=0';
sz = size(ph);
% phase_shift = 0.1999;

for k=1:sz(3)/2;
    ph_fft(:,:,k+1,:) = ph_fft(:,:,k+1,:).*exp(1i.*k*phase_shift);
    ph_fft(:,:,end-k+1,:) = ph_fft(:,:,end-k+1,:).*exp(-1i.*k*phase_shift);
end



ph = real(fft(ph_fft,[],3));
ph =unwrap(ph,[],3);
% ph = mean(ph,4);
% imdisp(ph (:,:,:),'size',[nan 8])
%%
% imdisp(mean(ph,4),'size',[nan 8]);