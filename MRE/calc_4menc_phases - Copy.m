function PHs = calc_4menc_phases(cimgs)

PHx = calc_4menc_phases_helper (cimgs,1, 2);
PHy = calc_4menc_phases_helper (cimgs,3, 4);
PHz = 

phz_pos = calc_ph_first_harmonic(cimgs(:,:,:,:,5),[],4);
phz_neg = calc_ph_first_harmonic(cimgs(:,:,:,:,6),[],4);
phz_add = calc_ph_first_harmonic(cimgs(:,:,:,:,5).*cimgs(:,:,:,:,6),[],4);
phz_sub = calc_ph_first_harmonic(cimgs(:,:,:,:,5),cimgs(:,:,:,:,6),4);
Phz_cplx = cat(6,ph_add,-ph_neg,ph_pos,ph_sub);
PHz_r = real(Phz_cplx);%real
PHz_i = imag(Phz_cplx);%imag
PHz = cat(4,PHz_r ,PHz_i);

function PH = calc_4menc_phases_helper (cimgs,id_pos, id_neg)
% id_pos = 1;
% id_neg = 2;

phz = calc_ph_first_harmonic(cimgs(:,:,:,:,5).*cimgs(:,:,:,:,6),[],4 );

ph_pos = calc_ph_first_harmonic(cimgs(:,:,:,:,id_pos),[],4);
ph_pos = ph_pos- phz./2;
ph_pos = ph_pos ./ (menc_pos_cplx(1,1)./abs(menc_pos_cplx(1,1)));
 

ph_neg = calc_ph_first_harmonic(cimgs(:,:,:,:,id_neg),[],4);
ph_neg = ph_neg - phz./2;
ph_neg = ph_neg ./ (menc_neg_cplx(1,1)./abs(menc_neg_cplx(1,1)));
ph_neg = -ph_neg;

ph_add = calc_ph_first_harmonic(cimgs(:,:,:,:,id_pos).*cimgs(:,:,:,:,id_neg).*conj(cimgs(:,:,:,:,5)).*conj(cimgs(:,:,:,:,6)),[],4);
ph_add = ph_add .* (menc_add(1,1)./abs(menc_add(1,1)));

ph_sub = calc_ph_first_harmonic(cimgs(:,:,:,:,id_pos),cimgs(:,:,:,:,id_neg),4);

Ph_cplx = cat(6,ph_add,-ph_neg,ph_pos,ph_sub);
PH_r = real(Ph_cplx);%real
PH_i = imag(Ph_cplx);%imag
PH = cat(4,PH_r ,PH_i);