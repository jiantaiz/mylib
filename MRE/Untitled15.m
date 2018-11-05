




% e13143s7 = load('\\mr-cimstore\mre-cim\Ziying\For_Yi\2016_0927_e13143_ZY_5+60Hz\s7_5hz=9T_60hz=4T_pos2\cimgs.mat');

e13143s4 = load('\\mr-cimstore\mre-cim\Ziying\For_Yi\2016_0927_e13143_ZY_5+60Hz\s4_5hz=9T_60hz=4T_pos1\cimgs.mat');

cimgs = e13143s4.cimgs;

mask = abs(cimgs) > 50;

mask = all(mask,5);

% cimgs = cimgs.*mask;

small_menc = 1/6*exp(1i*-1.1);
sz = size(cimgs);
sz = sz(1:end-1);
ph_uw_z = zeros(sz);
ph_small_z=zeros(sz);
ph_big_z=zeros(sz);
for k=1:9
    time_offset = (k-1)*4+1 : (k-1)*4+4 ;
    [my_ph_uw, ph_big, ph_small] = unwrap_tetra(cimgs(:,:,:,time_offset,:),small_menc);
    
    
    ph_uw_z(:,:,:,time_offset) = my_ph_uw;
    ph_small_z(:,:,:,time_offset)=ph_small;
    ph_big_z(:,:,:,time_offset)=ph_big;
end
%%
my_ph_uw = my_ph_uw.*mask(:,:,:,time_offset);

imdisp(cat(2,my_ph_uw(:,:,:),ph_big(:,:,:),ph_small(:,:,:)))
%%
ph_uw_z2=ph_uw_z.*mask;
ph_small_z2=ph_small_z.*mask;

sl = 1;
time_offset = 1:36;
imdisp(cat(2,ph_uw_z2(:,:,sl,time_offset),ph_big_z(:,:,sl,time_offset),ph_small_z2(:,:,sl,time_offset)))


sz = size(ph_uw_z2);
%%
ph_uw_z3 = reshape(ph_uw_z2,[128,128,4,4,9]);
ph_small_z3 = reshape(ph_small_z2,[128,128,4,4,9]);

%%

ph_uw_z3_fft = ifft(ph_uw_z3,[],4);
ph_small_z3_fft = ifft(ph_small_z3,[],4);
%%
sl=2
imdisp( abs(cat(2,ph_uw_z3_fft(:,:,sl,2,:),ph_small_z3_fft(:,:,sl,2,:))))
    
%%  
imdisp( angle(cat(2,ph_uw_z3_fft(:,:,sl,2,:),ph_small_z3_fft(:,:,sl,2,:))))


