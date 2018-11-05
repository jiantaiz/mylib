t = 1:4;
phase_shift =0;
a= 2 ;
b= 0 ;
c = pi/2;
sl = 3;
clear small_ph big_ph
mask = abs(cimgs) > 100;
cimgs= cimgs.*mask;
for dir =[1 3 5]
%     small_ph (:,:,:,:,(dir+1)/2) = calc_small_ph_tetra(cimgs(:,:,:,t,dir:dir+1),phase_shift,a,b,c);
    
%     img = cimgs(:,:,:,t,dir).*(cimgs(:,:,:,t,dir+1));
     img = cimgs(:,:,:,t,dir);
    img = sqrt(abs(img)).*exp(1i.*angle(img) ) .* exp(1i.*0);
    small_ph(:,:,:,:,(dir+1)/2) = angle(img);
%     small_ph(:,:,:,:,(dir+1)/2) = remove_phase_dc(img,4);
%     small_ph2(:,:,:,:,(dir+1)/2)=angle(rm_linear_ph(), a,b,c)); 

    img = cimgs(:,:,:,t,dir).* conj(cimgs(:,:,:,t,dir+1));
    img = sqrt(abs(img)).*exp(1i.*angle(img) ).* exp(1i.*pi);
    big_ph(:,:,:,:,(dir+1)/2)=remove_phase_dc(img,4);
    
    
    
%     big_ph(:,:,:,:,(dir+1)/2) = angle(cimgs(:,:,:,t,dir).* conj(cimgs(:,:,:,t,dir+1)));
    
    
    
end
%%


small_ph_fft = 2*ifft(small_ph,[],4);
big_ph_fft  = 2*ifft(big_ph,[],4);

% small_ph_fft = small_ph_fft(:,:,:,2,:);


phx = small_ph_fft(:,:,:,2,1);
phy = small_ph_fft(:,:,:,2,2);
phz = small_ph_fft(:,:,:,2,3); % in Forier domain

big_ph_est(:,:,:,1,1) = (phx - pol2cplx(0.648, 6.16).* phz) .* pol2cplx(18.083, 6.76) ;
big_ph_est(:,:,:,1,2) = (phy - pol2cplx(0.648, 6.16).* phz) .* 18.88 ;
big_ph_est(:,:,:,1,3) = phz .* pol2cplx(6.8707, 11.092);

%%
sl=5
imdisp(real(cat(2,small_ph_fft(:,:,sl,2,:),big_ph_fft(:,:,sl,2,:),0*big_ph_est(:,:,sl,1,:) )),[-pi pi])
imdisp(imag(cat(2,small_ph_fft(:,:,sl,2,:),big_ph_fft(:,:,sl,2,:),0*big_ph_est(:,:,sl,1,:) )),[-pi pi])
imdisp(abs(cat(2,small_ph_fft(:,:,sl,2,:),big_ph_fft(:,:,sl,2,:),0*big_ph_est(:,:,sl,1,:) )),[-pi pi])


%%
N=8
dir =2;
clear im im2 im3 im4;
sl = 1

for k = 1:N
%     im(:,:,k) = big_ph_est(:,:,sl,1,dir).*exp(-1i*(k-1)*2*pi/N);
     im2(:,:,k) = big_ph_fft(:,:,sl,2,dir).*exp(-1i*(k-1)*2*pi/N);
    im3(:,:,k) = small_ph_fft(:,:,sl,2,dir).*exp(-1i*(k-1)*2*pi/N);
    im4(:,:,k) = -img_fft(:,:,sl)/2.*exp(1i*(k-1)*2*pi/N);
end
imdisp(cat(2,real(im2),real(im3),real(im4)));

imdisp(imag(cat(2,im2,im3,im4)));



