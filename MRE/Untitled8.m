
load('cimgs.mat');
load('phs.mat')
clear small_ph big_ph big_ph_est
mask = abs(cimgs) > 90;
cimgs= cimgs.*mask;
phs_x =phs_x.*mask(:,:,:,:,1);
phs_y =phs_y.*mask(:,:,:,:,1);
phs_z =phs_z.*mask(:,:,:,:,1);
% cimgs(:,:,:,t,4) =  cimgs(:,:,:,t,4).*(-1);  % y-dir rotate pi
t=1:4
 
 sl=4
big_ph = angle(cimgs(:,:,:,:,1:2:end).*conj(cimgs(:,:,:,:,2:2:end)));
small_ph = angle(cimgs(:,:,:,:,1:2:end).*cimgs(:,:,:,:,2:2:end));

imdisp(big_ph(:,:,sl,:,:),[-2*pi, 2*pi],'size',[3,4]);set(gcf,'name', 'big_ph'); colormap (awave)
imdisp(small_ph(:,:,sl,:,:),[-2*pi, 2*pi],'size',[3,4]);set(gcf,'name', 'small_ph'); colormap (awave)


%%


% for dir =[1 3 5]
%         img1 = cimgs(:,:,:,t,dir);
%         img2 = cimgs(:,:,:,t,dir+1);
%         big_ph1f(:,:,:,(dir+1)/2) = calc_ph_first_harmonic(img1,img2,4);
%         small_ph1f(:,:,:,(dir+1)/2) = calc_ph_first_harmonic(img1,conj(img2),4);
% 
% end

% 
% small_ph1f(:,:,:,2) = shift_phase(small_ph1f(:,:,:,2), pi);% Y remove a pi constant DC component.
% phx = small_ph1f(:,:,:,1);
% phy = small_ph1f(:,:,:,2);
% phz = small_ph1f(:,:,:,3); % in Forier domain

% big_ph_est(:,:,:,1) = (phx - pol2cplx(0.648, 6.16).* phz);% .* pol2cplx(18.083, 6.76) ;
% big_ph_est(:,:,:,2) = (phy - pol2cplx(0.648, 6.16).* phz);% .* 18.88 ;
% big_ph_est(:,:,:,3) = phz; %.* pol2cplx(6.8707, 11.092);

% big_ph_est(:,:,:,1) = (phx - pol2cplx(0.648, 6.16).* phz).* pol2cplx(18.083, 6.76) ;
% big_ph_est(:,:,:,2) = (phy - pol2cplx(0.648, 6.16).* phz) .* 18.88 ;
% big_ph_est(:,:,:,3) = phz .* pol2cplx(6.8707, 11.092);

%  big_ph_est(:,:,:,1) = (phx - phz).*conj(menc_sub(1)./(menc_add(1))) ;
%  big_ph_est(:,:,:,2) = (phy - phz).*conj(menc_sub(2)./(menc_add(2)));
%   big_ph_est(:,:,:,3) = phz.* conj(menc_sub(3,3)./(menc_add(3,3))) ;
 
 
 
%  big_ph_est(:,:,:,1) = (phx - phz).*T(1,1) ;
%  big_ph_est(:,:,:,2) = (phy - phz).*T(2,2);
%  big_ph_est(:,:,:,3) = phx.*T(3,1) + phy.*T(3,2) + phz.*T(3,3);
%  big_ph_est(:,:,:,3) = big_ph_est(:,:,:,3).*exp(1i.* -0.45*pi);
 
 
 
%  big_ph_est(:,:,:,3) = phz.*(-T(3,3))  ;
 
 
%  imdisp(cat(2,real(small_ph1f(:,:,sl,:)),real(big_ph1f(:,:,sl,:)),real(big_ph_est(:,:,sl,:))))
%  set(gcf, 'name', 'first harmonic (real): small | big | big_est');
 %%
%  imdisp(cat(2,abs(small_ph1f(:,:,sl,:)),abs(big_ph1f(:,:,sl,:)),abs(big_ph_est(:,:,sl,:))  ))
%  set(gcf, 'name', 'first harmonic (abs): small | big | big_est');
 %%
%  imdisp(cat(2,imag(small_ph1f(:,:,sl,:)),imag(big_ph1f(:,:,sl,:)),imag(big_ph_est(:,:,sl,:))  ))
%  set(gcf, 'name', 'first harmonic (imag): small | big | big_est');
%%



phs_x_est = angle(cimgs(:,:,:,:,1).*cimgs(:,:,:,:,2).*conj(cimgs(:,:,:,:,5)).*conj(cimgs(:,:,:,:,6)));
% phs_x_est = phs_x_est.*T(1,1);
% big_ph_est(:,:,:,1) = calc_ph_first_harmonic(phs_x_est,[],4).*T(1,1);
big_ph_est(:,:,:,1) = phs_x_est(:,:,:,1)-phs_x_est(:,:,:,3) + 1i*(phs_x_est(:,:,:,2)-phs_x_est(:,:,:,4));
big_ph_est(:,:,:,1) = big_ph_est(:,:,:,1).*T(1,1);

phs_y_est = angle( cimgs(:,:,:,:,3).*cimgs(:,:,:,:,4).*conj(cimgs(:,:,:,:,5)).*conj(cimgs(:,:,:,:,6)));
% phs_y_est = phs_y_est *T(2,2);
% big_ph_est(:,:,:,2) = calc_ph_first_harmonic(exp(1i*angle(phs_y_est)),[],4) *T(2,2);
big_ph_est(:,:,:,2) = phs_y_est(:,:,:,1)-phs_y_est(:,:,:,3) + 1i*(phs_y_est(:,:,:,2)-phs_y_est(:,:,:,4));
big_ph_est(:,:,:,2) = big_ph_est(:,:,:,2).*T(2,2);


phs_z_est = angle(cimgs(:,:,:,:,5).*cimgs(:,:,:,:,6));
big_ph_est(:,:,:,3) = phs_z_est(:,:,:,1)-phs_z_est(:,:,:,3) + 1i*(phs_z_est(:,:,:,2)-phs_z_est(:,:,:,4));

big_ph_est(:,:,:,3) = calc_ph_first_harmonic(cimgs(:,:,:,:,5),conj(cimgs(:,:,:,:,6)),4);

big_ph_est(:,:,:,3) = big_ph_est(:,:,:,3).*conj(T(3,3)) + big_ph_est(:,:,:,2)./T(2,2).*conj(T(3,2)) + big_ph_est(:,:,:,1)./T(1,1).*conj(T(3,1)) ;


for k = 1:4
    phs_x_est(:,:,:,k) = 0.5* big_ph_est(:,:,:,1).*exp(1i.*-(k-1)*pi/2);
    phs_y_est(:,:,:,k) = 0.5* big_ph_est(:,:,:,2).*exp(1i.*-(k-1)*pi/2);
    phs_z_est(:,:,:,k) = 0.5* big_ph_est(:,:,:,3).*exp(1i.*-(k-1)*pi/2);
end



tol = 1.1

    phs_x_uw = myunwrap(phs_x,real(phs_x_est), 1,[],tol);
    phs_y_uw = myunwrap(phs_y,real(phs_y_est), 1,[],tol);
    phs_z_uw = myunwrap(phs_z,real(phs_z_est), 1,[],tol);
    
    s_small=load('\\mr-cimstore\mre-cim\Ziying\For_Yi\2017_0301_e52212_PVA_Menc\s20_60pct_0.2meg\phs.mat');
    
    s_small.phs_x =s_small.phs_x.*mask(:,:,:,:,1);
    s_small.phs_y =s_small.phs_y.*mask(:,:,:,:,1);
    s_small.phs_z =s_small.phs_z.*mask(:,:,:,:,1);
    
    menc_big = 1.6;
    menc_small = 0.4;
    scale=[1 1 1];
    phs_x_uw2 = myunwrap(phs_x,scale(1).*menc_big/menc_small*s_small.phs_x, 1);
    phs_y_uw2 = myunwrap(phs_y,scale(2).*menc_big/menc_small*s_small.phs_y, 1);
    phs_z_uw2 = myunwrap(phs_z,scale(3).*menc_big/menc_small*s_small.phs_z, 1);
%     imdisp(cat(2,phs_x(:,:,4,:),real(phs_x_est(:,:,4,:)),real(phs_x_uw(:,:,4,:))))
%     imdisp(cat(2,phs_y(:,:,4,:),real(phs_y_est(:,:,4,:)),real(phs_y_uw(:,:,4,:))))
%     imdisp(cat(2,phs_z(:,:,4,:),real(phs_z_est(:,:,4,:)),real(phs_z_uw(:,:,4,:))))
%     
%%
 load('phs_3duw.mat')
  imdisp(cat(2,phs_x(:,:,sl,:),scale(1).*menc_big/menc_small*s_small.phs_x(:,:,sl,:),real(phs_x_uw2(:,:,sl,:)),real(phs_x_est(:,:,sl,:)),real(phs_x_uw(:,:,sl,:)), x_uw(:,:,sl,:)),[-9 9]);set(gcf,'name','X: big | small | unwrap | small_dual | unwrap_dual')
  colormap(awave)
  imdisp(cat(2,phs_y(:,:,sl,:),scale(2).*menc_big/menc_small*s_small.phs_y(:,:,sl,:),real(phs_y_uw2(:,:,sl,:)),real(phs_y_est(:,:,sl,:)),real(phs_y_uw(:,:,sl,:)), y_uw(:,:,sl,:)),[-9 9]);set(gcf,'name','Y: big | small | unwrap | small_dual | unwrap_dual')
  colormap(awave)
  imdisp(cat(2,phs_z(:,:,sl,:),scale(3).*menc_big/menc_small*s_small.phs_z(:,:,sl,:),real(phs_z_uw2(:,:,sl,:)),real(phs_z_est(:,:,sl,:)),real(phs_z_uw(:,:,sl,:)), z_uw(:,:,sl,:)),[-9 9]);set(gcf,'name','Z: big | small | unwrap | small_dual | unwrap_dual')
  colormap(awave)
  
    
    %%
    imdisp(cat(2,phs_x(:,:,sl,:),1.7/0.3*s22.phs_x(:,:,sl,:),real(phs_x_est(:,:,sl,:)),real(phs_x_uw(:,:,sl,:))))
%%
    imdisp(cat(2,phs_y(:,:,sl,:),1.7/0.3*s22.phs_y(:,:,sl,:),real(phs_y_est(:,:,sl,:)),real(phs_y_uw(:,:,sl,:))))
    
    
    
