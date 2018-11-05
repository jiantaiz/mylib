function ph = calc_small_ph_tetra(cimgs,phase_shift,a,b,c)
% cimgs [x,y,z,timeoffset,dir]
% dir = 5;
% ydata1 =squeeze( cimgs(:,:,sl,:,dir));% positive meg
% ydata2 =squeeze( cimgs(:,:,sl,:,dir+1)); % negtive meg
% ydata = squeeze( cimgs(:,:,sl,:,:)); % [x,y,z,timeoffset,dir]negtive meg

%a = 4 ;b=-1; c=0;
ydata = rm_linear_ph(prod(cimgs,5),a,b,c);
%%
ph = angle(ydata);
% ph = angle(ydata12(:,:,:,:)); 
% ph = unwrap(ph,[],4);
ph_fft = ifft(ph,[],4);
ph_fft(:,:,:,1,:)=0;
disp 'set small phase dc=0';
sz = size(ph);
% phase_shift = 0.1999;

for k=1:sz(4)/2;
    ph_fft(:,:,:,k+1,:) = ph_fft(:,:,:,k+1,:).*exp(1i.*k*phase_shift);
    ph_fft(:,:,:,end-k+1,:) = ph_fft(:,:,:,end-k+1,:).*exp(-1i.*k*phase_shift);
end

ph = real(fft(ph_fft,[],4));
% ph =unwrap(ph,[],4);

% imdisp(mean(ph,4),'size',[nan 8]);