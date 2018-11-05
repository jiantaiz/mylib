function [phs_uw, phs, phs_est] = unwrap_dual_menc2(cimgs,T,N_medfil)
% UNWRAP_DUAL_MENC unwrap phase maps from dual menc aquisition
%    [phs_uw, phs, phs_est] = UNWRAP_DUAL_MENC(cimgs,T,N_medfil)
%    INPUT:
%       cimgs(x,y,z,t,dir): contains complex images.
%       T:  3x3 the transform matrix from low menc to high menc 
%       i.e.  phs_hi = T*phs_low
%       N_medfil: number of times to apply median filter.
%    OUTPUT:
%       phs_uw: unwrapped phase maps
%       phs: wrapped phase maps
%       phs_est: estimated phase maps from low menc
%    Limitation: only work with 4 time offsets acq.
%    See also: calc_menc

% AUTHOR    : Yi Sui
% DATE      : 06/09/2017
%%
if nargin<3
    N_medfil=0;
end
cimgs= double(cimgs);
phs = angle(cimgs(:,:,:,:,1:2:end).*conj(cimgs(:,:,:,:,2:2:end)));
% 
%  phs_low_menc = calc_ph_first_harmonic(cimgs(:,:,:,:,1:2:end).*cimgs(:,:,:,:,2:2:end),[],4);
% phs_low_menc = ordfilt3D_2(phs_low_menc,14);
% [phs_low_menc,perm,nshifts] = shiftdata(phs_low_menc,5);
% 
% phs_hi_est = conj(T)*phs_low_menc(:,:); % matrix time
% phs_hi_est = reshape(phs_hi_est,size(phs_low_menc));
% phs_hi_est = unshiftdata(phs_hi_est,perm,nshifts);
% phs_hi_est = ordfilt3D_2(phs_hi_est,14);

phs_x_est = (cimgs(:,:,:,:,1).*cimgs(:,:,:,:,2).*conj(cimgs(:,:,:,:,5)).*conj(cimgs(:,:,:,:,6))); % x-z
big_ph_est(:,:,:,1) = calc_ph_first_harmonic(phs_x_est,[],4).*T(1,1);

dphs_est(:,:,:,:,1) =angle( phs_x_est.*conj(phs_x_est(:,:,:,[1 1:end-1]))) .*abs(T(1,1));


phs_y_est = cimgs(:,:,:,:,3).*cimgs(:,:,:,:,4).*conj(cimgs(:,:,:,:,5)).*conj(cimgs(:,:,:,:,6));%y-z
big_ph_est(:,:,:,2) = calc_ph_first_harmonic(phs_y_est,[],4).*T(2,2);

dphs_est(:,:,:,:,2) =angle( phs_y_est.*conj(phs_y_est(:,:,:,[1 1:end-1]))).*abs(T(2,2));


big_ph_est(:,:,:,3) = calc_ph_first_harmonic(cimgs(:,:,:,:,5),conj(cimgs(:,:,:,:,6)),4);
% II = big_ph_est(:,:,:,3).*(T(3,3)) + big_ph_est(:,:,:,2)./T(2,2).*conj(T(3,2)) +  big_ph_est(:,:,:,1)./T(1,1).*conj(T(3,1)) ;
big_ph_est(:,:,:,3)= big_ph_est(:,:,:,3).*conj(T(3,3));% + big_ph_est(:,:,:,2)./T(2,2).*conj(T(3,2)) +  big_ph_est(:,:,:,1)./T(1,1).*conj(T(3,1)) ;
% big_ph_est(:,:,:,3) = real(RR)+1i*imag(RR);
% big_ph_est(:,:,:,3) = big_ph_est(:,:,:,3).*(T(3,3)) + big_ph_est(:,:,:,2).*conj(T(3,2)) + big_ph_est(:,:,:,1).*conj(T(3,1)) ;

% big_ph_est(:,:,:,3) = big_ph_est(:,:,:,3).*(T(3,3)) + big_ph_est(:,:,:,2)./T(2,2).*(T(3,2)) + big_ph_est(:,:,:,1)./T(1,1).*(T(3,1)) ;
% big_ph_est(:,:,:,3) = big_ph_est(:,:,:,3).* (T(3,3)) + phs_low_menc(:,:,:,2).*(T(3,2)) + phs_low_menc(:,:,:,1).*(T(3,1)) ;
% big_ph_est(:,:,:,3) = big_ph_est(:,:,:,3).* (T(3,3)) ;

% phs_hi_est = ordfilt3D_2(big_ph_est,14);
phs_hi_est  = big_ph_est;
for k = 1:4
    phs_est(:,:,:,k,:) = 0.5 * phs_hi_est.*exp(1i.*-(k-1)*pi/2);
end
phs_est = real(phs_est);
% dphs_est = phs_est - phs_est (:,:,:,[1 1:end-1],:);
dphs_est(:,:,:,:,3)= phs_est(:,:,:,:,3) - phs_est (:,:,:,[1 1:end-1],3);
dphs_est = ordfilt3D_2(dphs_est,14);
tol = 1;
% phs_uw = myunwrap(phs,phs_est, 1,[],tol);
phs_uw = phs;
for k=2:4;
    phs_uw(:,:,:,k,:) = myunwrap(phs_uw(:,:,:,k,:),phs_uw(:,:,:,k-1,:)+dphs_est(:,:,:,k,:), 1,[],tol);
end

jmp=round(mean(phs_uw,4)./2./pi );
phs_uw = phs_uw-2*pi*jmp(:,:,:,ones(1,4),:);
for n=1:N_medfil
phs_uw_mf = ordfilt3D_2(phs_uw,14);
phs_uw = myunwrap(phs,phs_uw_mf, 1,[],tol);
end


