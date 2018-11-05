

sl =4;
phase_shift = 0;
dir =2;
timeoffset = 1:36;

big_ph = angle(cimgs(:,:,:,timeoffset,dir) .* conj(cimgs(:,:,:,timeoffset,dir+1)));

imdisp(big_ph(:,:,sl,:))

%%
dir =1;
sl =4;
timeoffset = 1:36;
%%
a = 1 ;b=0;phi=+1;

cimgs = rm_linear_ph(e13181s5.cimgs,a,b,phi);


small_ph_x = angle(cimgs(:,:,:,timeoffset,1) .* (cimgs(:,:,:,timeoffset,2)));
small_ph_y = angle(cimgs(:,:,:,timeoffset,3) .* (cimgs(:,:,:,timeoffset,4)));
small_ph_z = angle(cimgs(:,:,:,timeoffset,5) .* (cimgs(:,:,:,timeoffset,6)));
%
imdisp(( cat(2, small_ph_x(:,:,2,:),small_ph_y(:,:,2,:),small_ph_z(:,:,2,:),abs(cimgs(:,:,2,:,1))/1000*13)))
%%
imdisp(small_ph(:,:,sl,:))

%%
ph = calc_small_ph2(cimgs,sl,phase_shift,dir,timeoffset);


imdisp(ph)

%%
for offset_8hz = 1:9
    t = (offset_8hz-1)*4+1 :(offset_8hz-1)*4+4;
    for dir = [1 3 5]
    small_ph(:,:,:,t,(dir+1)/2) = calc_small_ph_tetra(cimgs(:,:,:,t,dir:dir+1),0,1,0,2);
    end
end
