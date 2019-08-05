function [phs_uw, phs, phs_est] = unwrap_dual_menc(cimgs,T,N_medfil,high_vib,z_filter)
% UNWRAP_DUAL_MENC unwrap phase maps from dual menc aquisition
%    [phs_uw, phs, phs_est] = UNWRAP_DUAL_MENC(cimgs,T,N_medfil,high_vib,z_filter)
%    INPUT:
%       cimgs(x,y,z,t,dir): contains complex images.
%       T:  3x3 the transform matrix from low menc to high menc
%       i.e.  phs_hi = T*phs_low
%       N_medfil: number of times to apply median filter.
%       high_vib: 1x3 indicates which direction (matrix dimention) has high vibration
%       z_filter: apply z_filter or not
%    OUTPUT:
%       phs_uw: unwrapped phase maps
%       phs: wrapped phase maps
%       phs_est: estimated phase maps from low menc
%    Limitation: only work with 4 time offsets acq. calc_ph_first_harmonic
%    halves dynamic range.
%    See also: calc_menc

% AUTHOR    : Yi Sui
% DATE      : 06/09/2017
%%
if nargin<3
    N_medfil=0;
end
if nargin<4
    high_vib = [0 0 0];
end

if numel(high_vib) == 1
    high_vib = [high_vib high_vib high_vib];
end
if nargin<5
    z_filter= 0;
end

cimgs = single(cimgs);
sz  = size(cimgs); %sz=[nx,ny,nz,nt,ndir]

if sz(4)~=4 && sz(4)~=8
    error('only support time offsets of 4 or 8');
    return;
end
if sz(4) == 8
    [phs_uw1, phs1, phs_est1] = unwrap_dual_menc(cimgs(:,:,:,1:2:end,:),T,N_medfil,high_vib,z_filter);
    [phs_uw2, phs2, phs_est2] = unwrap_dual_menc(cimgs(:,:,:,2:2:end,:),T,N_medfil,high_vib,z_filter);
    phs_uw = cat(4,phs_uw1,phs_uw2);
    phs_uw = phs_uw(:,:,:,[1 5 2 6 3 7 4 8],:);
    if nargout>1
        phs = cat(4,phs1,phs2);
        phs = phs(:,:,:,[1 5 2 6 3 7 4 8],:);
    end
    if nargout>2
        phs_est = cat(4,phs_est1,phs_est2);
        phs_est = phs_est(:,:,:,[1 5 2 6 3 7 4 8],:);
    end
    return;
end

phs = angle(cimgs(:,:,:,:,1:2:end).*conj(cimgs(:,:,:,:,2:2:end)));


%  cimgs = lpfilter(cimgs,3,2.5);
%  cimgs = lpfilter(cimgs,2,2);
% cimgs = lpfilter(cimgs,3,2);

%
%  phs_low_menc = calc_ph_first_harmonic(cimgs(:,:,:,:,1:2:end).*cimgs(:,:,:,:,2:2:end),[],4);
% phs_low_menc = ordfilt3D_2(phs_low_menc,14);
% [phs_low_menc,perm,nshifts] = shiftdata(phs_low_menc,5);
%
% phs_hi_est = conj(T)*phs_low_menc(:,:); % matrix time
% phs_hi_est = reshape(phs_hi_est,size(phs_low_menc));
% phs_hi_est = unshiftdata(phs_hi_est,perm,nshifts);
% phs_hi_est = ordfilt3D_2(phs_hi_est,14);


big_ph_est = complex( zeros( [sz(1:3) 3]));
phs_x_est = (cimgs(:,:,:,:,1).*cimgs(:,:,:,:,2).*conj(cimgs(:,:,:,:,5)).*conj(cimgs(:,:,:,:,6))); % x-z

phs_x_est = lpfilter(exp(1i*angle(phs_x_est)),3,2);

if high_vib(1) %if low menc already has wrapping along time dimention when using calc_ph_first_harmonic which doubles the phase;
%     phs_x_est = unwrap(angle(phs_x_est) ,[],4);
    phs_x_est = angle(phs_x_est);
    big_ph_est(:,:,:,1) = phs_x_est(:,:,:,1)-phs_x_est(:,:,:,3) + 1i*(phs_x_est(:,:,:,2)-phs_x_est(:,:,:,4));
    big_ph_est(:,:,:,1) = big_ph_est(:,:,:,1).*T(1,1);
else
    big_ph_est(:,:,:,1) = calc_ph_first_harmonic(phs_x_est,[],4).*T(1,1);
end

phs_y_est = cimgs(:,:,:,:,3).*cimgs(:,:,:,:,4).*conj(cimgs(:,:,:,:,5)).*conj(cimgs(:,:,:,:,6));%y-z
phs_y_est = lpfilter(exp(1i*angle(phs_y_est)),3,2);

if high_vib(2)
%     phs_y_est = unwrap(angle(phs_y_est) ,[],4);
    phs_y_est = angle(phs_y_est);
%     phs_y_est = unwrap3D_ssh(phs_y_est,1);
    big_ph_est(:,:,:,2) = phs_y_est(:,:,:,1)-phs_y_est(:,:,:,3) + 1i*(phs_y_est(:,:,:,2)-phs_y_est(:,:,:,4));
    big_ph_est(:,:,:,2) = big_ph_est(:,:,:,2).*T(2,2);
else
    big_ph_est(:,:,:,2) = calc_ph_first_harmonic(phs_y_est,[],4).*T(2,2);
end


phs_z_est = cimgs(:,:,:,:,5).*cimgs(:,:,:,:,6);
phs_z_est = lpfilter(exp(1i*angle(phs_z_est)),3,2);

if high_vib(3)
    
    phs_z_est = angle(phs_z_est);
    phs_z_est = phs_z_est - phs_z_est(:,:,:,[1 1:end-1],:);
    phs_z_est = unwrap(phs_z_est,[],4);
    big_ph_est(:,:,:,3) = phs_z_est(:,:,:,1)-phs_z_est(:,:,:,3) + 1i*(phs_z_est(:,:,:,2)-phs_z_est(:,:,:,4));
    big_ph_est(:,:,:,3) = big_ph_est(:,:,:,3).* conj(T(3,3));
else
    big_ph_est(:,:,:,3) = calc_ph_first_harmonic(phs_z_est,[],4).* conj(T(3,3));
%     big_ph_est(:,:,:,3) = calc_ph_first_harmonic(phs_z_est,[],4).* T(3,3);
    
end



%  big_ph_est(:,:,:,3) = big_ph_est(:,:,:,3) + big_ph_est(:,:,:,2)./T(2,2).*(T(3,2)) + big_ph_est(:,:,:,1)./T(1,1).*(T(3,1)) ;
 big_ph_est(:,:,:,3) = big_ph_est(:,:,:,3) + big_ph_est(:,:,:,2)./T(2,2).*(T(3,2)) + big_ph_est(:,:,:,1)./T(1,1).*(T(3,1)) ;


% II = big_ph_est(:,:,:,3).*(T(3,3)) + big_ph_est(:,:,:,2)./T(2,2).*conj(T(3,2)) +  big_ph_est(:,:,:,1)./T(1,1).*conj(T(3,1)) ;
% big_ph_est(:,:,:,3)= big_ph_est(:,:,:,3).*conj(T(3,3));% + big_ph_est(:,:,:,2)./T(2,2).*conj(T(3,2)) +  big_ph_est(:,:,:,1)./T(1,1).*conj(T(3,1)) ;
% big_ph_est(:,:,:,3) = real(RR)+1i*imag(RR);
% big_ph_est(:,:,:,3) = big_ph_est(:,:,:,3).*(T(3,3)) + big_ph_est(:,:,:,2).*conj(T(3,2)) + big_ph_est(:,:,:,1).*conj(T(3,1)) ;

% big_ph_est(:,:,:,3) = big_ph_est(:,:,:,3).*(T(3,3)) + big_ph_est(:,:,:,2)./T(2,2).*(T(3,2)) + big_ph_est(:,:,:,1)./T(1,1).*(T(3,1)) ;
% big_ph_est(:,:,:,3) = big_ph_est(:,:,:,3).* (T(3,3)) + phs_low_menc(:,:,:,2).*(T(3,2)) + phs_low_menc(:,:,:,1).*(T(3,1)) ;
% big_ph_est(:,:,:,3) = big_ph_est(:,:,:,3).* (T(3,3)) ;

if z_filter == 1
%     big_ph_est = lpfilter(big_ph_est,3,2.5);
elseif z_filter == 2
    big_ph_est = ordfilt3D_2(big_ph_est,14);
elseif z_filter == 3
%     big_ph_est = lpfilter(big_ph_est,3,4);
%     big_ph_est = ordfilt3D_2(big_ph_est,14);
%     big_ph_est = medfilt3_nd(big_ph_est);
end
    phs_hi_est = big_ph_est;


% phs_hi_est = big_ph_est;
phs_est = zeros([sz(1:3) 4, 3]);
for k = 1:4
    phs_est(:,:,:,k,:) = real(0.5 * phs_hi_est.*exp(1i.*-(k-1)*pi/2));
end
% phs_est = real(phs_est);
tol = 1;
phs_uw = myunwrap(phs,phs_est, 1,[],tol);



% estimated phase may have a global background shift, which can cause
% errors on some local points. here I remove the shift from phs_est and
% unwrap again.
if(1)
    magn = mean(abs(cimgs(:,:,:,:)),4);
    sz = size(magn);
    sz = [sz(1:2) 1 sz(3)];
    magn = reshape(magn,sz)./max(magn(:));
    mask = im2bw(magn, 0.9*graythresh(magn));
    mask = mask(:,:,:);
%     ov(single(mask));
    dph = bsxfun(@times,(phs_uw-phs_est),mask);
    [~, ~, nv] = size(phs_uw);
    for k =1:nv
        tmp = dph(:,:,k); 
        tmp = tmp(tmp~=0);
        tmp = tmp(~isnan(tmp));
        dph_m = mean(tmp(:),'omitnan');
        phs_est(:,:,k) = phs_est(:,:,k)+dph_m;
    end
    phs_uw = myunwrap(phs,phs_est, 1,[],tol);
    clear mask magn dph
end

%median filter to remove isolated error points

for n=1:N_medfil
    if verLessThan('matlab','9.1.0')
        phs_uw_mf = ordfilt3D_2(phs_uw,14);
    else
        phs_uw_mf = medfilt3_nd(phs_uw);
    end
    phs_uw = myunwrap(phs,phs_uw_mf, 1,[],tol);
end