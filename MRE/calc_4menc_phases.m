function [PH, menc] = calc_4menc_phases(cimgs,menc_mtx)
% CALC_4MENC_PHASES calculates phase maps at 4 menc for multi_menc unwrap
%    [PH, menc] = CALC_4MENC_PHASES(cimgs,menc_mtx) 
%          INPUT:
%                 cimgs(x,y,z,t,dir+-), currently only work with 4 time
%                 offsets
%                 menc_mtx: 3x3x4 menc matrix at 4 different gradient settings
%    See also: myunwrap_multi_menc

% AUTHOR    : Yi Sui
% DATE      : 05/25/2017


% 1st(smallest) menc
% inv(menc_mtx(:,:,1))
% dx = (ph_addx - ph_addz)*k, so is dy
cimgs_org = cimgs;
% cimgs = low_pass_filter(cimgs,64);
m = menc_mtx(:,:,4)*inv(menc_mtx(:,:,1));% transform matrix form 1 to 4
menc(:,1) =1./ diag(abs(m));
ph_add(:,:,:,:,1) = calc_ph_first_harmonic(cimgs(:,:,:,:,1).*cimgs(:,:,:,:,2).*conj(cimgs(:,:,:,:,5)).*conj(cimgs(:,:,:,:,6)),[],4);%x
ph_add(:,:,:,:,2) = calc_ph_first_harmonic(cimgs(:,:,:,:,3).*cimgs(:,:,:,:,4).*conj(cimgs(:,:,:,:,5)).*conj(cimgs(:,:,:,:,6)),[],4);%y
%z
mz = m(3,:)./abs(m(3,3));
p1 = calc_ph_first_harmonic(cimgs(:,:,:,:,1).*cimgs(:,:,:,:,2),[],4);
p2 = calc_ph_first_harmonic(cimgs(:,:,:,:,3).*cimgs(:,:,:,:,4),[],4);
p3 = calc_ph_first_harmonic(cimgs(:,:,:,:,5).*cimgs(:,:,:,:,6),[],4);

ph_add(:,:,:,:,3)   = p3.*mz(3) ; %phase differenc correction, ignore contribution from x, y
% ph_add(:,:,:,:,3)   = p1.*mz(1) + p2.*mz(2) + p3.*mz(3) ; %phase differenc correction


% 2nd menc
m = menc_mtx(:,:,4)*inv(menc_mtx(:,:,2));% transform matrix form 2 to 4
menc(:,2) = 1./diag(abs(m));

ph_neg = calc_ph_first_harmonic(cimgs(:,:,:,:,2:2:end),[],4);% no correction, ignore small (~10%) contribution from z
for k = 1:3
    ph_neg (:,:,:,:,k) = ph_neg (:,:,:,:,k) .* (m(k,k)./ abs(m(k,k)));%phase correction
end


% 3rd menc
m = menc_mtx(:,:,4)*inv(menc_mtx(:,:,3));% transform matrix form 3 to 4
menc(:,3) = 1./diag(abs(m));

ph_pos = calc_ph_first_harmonic(cimgs(:,:,:,:,1:2:end),[],4);
for k = 1:3
    ph_pos (:,:,:,:,k) = ph_pos (:,:,:,:,k) .* (m(k,k)./ abs(m(k,k)));
end
% 4th(largest) menc
cimgs = cimgs_org;
m = menc_mtx(:,:,4)*inv(menc_mtx(:,:,4));% transform matrix form 4 to 4
menc(:,4) = 1./diag(abs(m));
ph_sub = calc_ph_first_harmonic(cimgs(:,:,:,:,1:2:end),cimgs(:,:,:,:,2:2:end),4);


Ph_cplx = cat(6,ph_add,ph_neg,ph_pos,ph_sub);
PH_r = real(Ph_cplx);%real
PH_i = imag(Ph_cplx);%imag
PH = cat(4,PH_r ,PH_i);
