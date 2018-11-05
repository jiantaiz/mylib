
dir = 1;

t = 1:32;

menc_ratio = 10;
phase_shift =0;
a= 0 ;
b= 0 ;
c = 0;
sl = 8;
for dir =[1 3 5]
    small_ph (:,:,:,:,(dir+1)/2) = calc_small_ph_tetra(cimgs(:,:,:,t,dir:dir+1),phase_shift,a,b,c);
end
%%
imdisp((small_ph(:,:,sl,t,2)  ));
%%
clear pdimg
for dir =[1 3 5]

    pdimg(:,:,:,:,(dir+1)/2) = cimgs(:,:,:,t,dir).* conj(cimgs(:,:,:,t,dir+1));
    

end
big_ph = angle(pdimg.*exp(1i.*pi));
big_ph_fft = ifft(big_ph,[],4);
big_ph_fft(:,:,:,1,:) = 0;
big_ph_2 = fft(big_ph_fft,[],4);

%%
imdisp(big_ph(:,:,sl,:,2));

%%

%%
menc_ratio = 7;
small_ph2 = small_ph.*menc_ratio;

uw = myunwrap(big_ph(:,:,:,t,3),small_ph2(:,:,:,t,3),1);
imdisp(cat(2,uw(:,:,sl,t),big_ph(:,:,sl,t,3),small_ph2(:,:,sl,t,3)))
  
%%
phi=[0,nan,pi,nan,0];
t = 1:32;
clear big_ph_2;
for dir =[1 3 5]

pdimg = cimgs(:,:,:,t,dir).* conj(cimgs(:,:,:,t,dir+1));
    


big_ph = angle(pdimg.*exp(1i.*phi(dir)));
big_ph_fft = ifft(big_ph,[],4);
big_ph_fft(:,:,:,1,:) = 0;
big_ph_2(:,:,:,:,(dir+1)/2) = fft(big_ph_fft,[],4);
end


%%
clear big_ph_2
[big_ph,pdimg] = calc_big_ph_tetra(cimgs,1);

phi=[0,pi,0];

for dir = [1 2 3]
    big_ph = angle(pdimg(:,:,:,:,dir).*exp(1i.*phi(dir)));
    big_ph_fft = ifft(big_ph,[],4);
    big_ph_fft(:,:,:,1,:) = 0;
    big_ph_2(:,:,:,:,dir) = fft(big_ph_fft,[],4);
end
%%

ph = reshape(ph_uw_5hz.ph_uw,[128 128 8 4 8 3]);
ph_fft_5hz = ifft(ph,[],4);

ph = reshape(ph_uw_0hz.ph_uw,[128 128 8 4 8 3]);
ph_fft_0hz = ifft(ph,[],4);
%%
sl = 6;
f = 2;
dir =3;


angdf = (ph_fft_5hz(:,:,sl,f,:,dir)) - (ph_fft_0hz(:,:,sl,f,:,dir))  ;
angdf = cat(2,ph_fft_0hz(:,:,sl,f,:,dir),ph_fft_5hz(:,:,sl,f,:,dir),angdf) ;

% imdisp(std(abs(angdf).^2,[],6))
imdisp(real(angdf))

%%
imdisp(ph(:,:,sl,:,1,dir))



