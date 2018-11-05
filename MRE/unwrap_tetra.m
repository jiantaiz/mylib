function [my_ph_uw, ph_big, ph_small] = unwrap_tetra(cimgs,small_menc)
% e13143s4 = load('\\mr-cimstore\mre-cim\Ziying\For_Yi\2016_0927_e13143_ZY_5+60Hz\s4_5hz=9T_60hz=4T_pos1\cimgs.mat');
phase_shift=-1.1;
phase_shift=angle(small_menc);
small_menc_abs = abs(small_menc);
% timeoffset = 17:20;
ph_small = calc_small_ph_tetra(cimgs,phase_shift);
ph_small = squeeze(ph_small)./small_menc_abs;

ph_big = calc_big_ph_tetra(cimgs,'balanced');
ph_big = squeeze(ph_big(:,:,:,:,3));%currently, only work on z-dirction.
% 
% N=2;
% fft_points = 8;
% v = [-N:N];
% for k = 1:fft_points
%     in{k} =v;% {v,v,v,v};
% end
% uw_kernal = combvec(in{:}).*2*pi;


my_ph_uw12 = myunwrap(ph_big,ph_small,1);
%%
my_ph_uw = remove_2pi_dc(my_ph_uw12,4);

%
my_ph_uw_tmp = my_ph_uw;
sz = size(my_ph_uw);
% threshold=0.1;
min_val=ones(sz(1:3))*10000;
for dc =-2*pi:pi/64:2*pi

    my_ph_uw_fft = ifft(my_ph_uw_tmp,[],4);
    m3=abs(my_ph_uw_fft(:,:,:,3:end-1)).^2;
    mask  = m3 < min_val*0.9;
    
    min_val(mask) = m3(mask);

    mask2 = repmat(mask,[1,1,1,size(my_ph_uw_fft,4)]);

    my_ph_uw(mask2) = my_ph_uw_tmp(mask2);
    
    my_ph_uw_tmp = myunwrap(ph_big,ph_small+dc,1);
end
%
my_ph_uw=remove_2pi_dc(my_ph_uw,4);
%%
% I1 = remove_2pi_dc(my_ph_uw12,4);
% I2 = remove_2pi_dc(my_ph_uw,4);
% imdisp(cat(2,I1(:,:,:),I2(:,:,:),mask2(:,:,:)*10))


