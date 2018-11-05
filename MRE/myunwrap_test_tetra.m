
% e13143s4 = load('\\mr-cimstore\mre-cim\Ziying\For_Yi\2016_0927_e13143_ZY_5+60Hz\s4_5hz=9T_60hz=4T_pos1\cimgs.mat');
e13143s7 = load('\\mr-cimstore\mre-cim\Ziying\For_Yi\2016_0927_e13143_ZY_5+60Hz\s7_5hz=9T_60hz=4T_pos2\cimgs.mat');

cimgs = e13143s7.cimgs;
sl=1:4;
phase_shift=-1.1;
timeoffset = 17:20;
ph_small = calc_small_ph_tetra(cimgs(:,:,sl,1:4,:),phase_shift);
ph_small = squeeze(ph_small).*6;

ph_big = calc_big_ph_tetra(cimgs,'balanced');
ph_big = squeeze(ph_big(:,:,sl,timeoffset,3));
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
I1 = remove_2pi_dc(my_ph_uw12,4);
I2 = remove_2pi_dc(my_ph_uw,4);
imdisp(cat(2,I1(:,:,:),I2(:,:,:),mask2(:,:,:)*10))



%%
try close(1000); end;
figure(1000);
plot([1 4], [pi -pi; pi -pi],'k--');
grid on;
hold on
plt = plot([1,2,3,4,5], ones(5));
legend({'pi','-pi','big menc','my unwrap','small menc'});
xlabel('time offsets');
% imdisp([ph_big,my_ph_uw, unwrap(my_ph_uw,[],3),ph_small,ph_small_uw,ph_small_2],[-2*pi 2*pi],'size',[nan,1])
% imdisp([ph_big,my_ph_uw1,my_ph_uw2,my_ph_uw12,my_ph_uw1+my_ph_uw2,ph_small_2,z_3d_uw2,unwrap(ph_big,[],3)],[-5*pi 5*pi],'size',[nan,1]);
k=ifft(my_ph_uw,[],3);k2 = ifft(ph_small,[],3);imdisp(abs([k k2]))
imdisp([ph_big,my_ph_uw12,my_ph_uw,ph_small],[-5*pi 5*pi],'size',[nan,1]);

colormap(awave)
linestyle={'-','-.','--','--'};
clear ydata
while 1
    [xi, yi, but] = ginput(1);      % get a point
    if ~isequal(but, 1)             % stop if not button 1
        break
    end
    xi = round(mod(xi,size(ph_big,2))); 
    yi = round(yi);
    
    ydata(:,1) =squeeze( ph_big(yi,xi,:));% 
    ydata(:,2) =squeeze( my_ph_uw12(yi,xi,:)); % 
    ydata(:,3) =squeeze( ph_small(yi,xi,:)); % 
    ydata(:,4) =squeeze( my_ph_uw(yi,xi,:)); % 
    xdata = [1:size(ydata,1)]';

    for k=1:4
    set(plt(k),'xdata',xdata,'ydata',ydata(:,k),'LineStyle',linestyle{k});
    end
    set(plt(1).Parent,'ylim',[-14,14]);
%       ydata_uw=squeeze(s3.x_ux(yi,xi,sl,:));
%     fft(ydata_uw)
%     set(plt2,'xdata',xdata,'ydata',ydata_uw)
%     
    
     ang = angle(ifft(ydata))
     mag = abs(ifft(ydata))
end












