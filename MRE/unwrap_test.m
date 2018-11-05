% cimgs =  fliplr( rot90(s0712.cimgs));
% cimgs =  fliplr( rot90(s8_cw.cimgs));
% e13143s4 = load('\\mr-cimstore\mre-cim\Ziying\For_Yi\2016_0927_e13143_ZY_5+60Hz\s4_5hz=9T_60hz=4T_pos1\cimgs.mat')
cimgs = s8_cw.cimgs;
% cimgs = e13143s4.cimgs;
sl=42;
dir = 5;
ydata1 =squeeze( cimgs(:,:,sl,:,dir));% positive meg
ydata2 =squeeze( cimgs(:,:,sl,:,dir+1)); % negtive meg
imdisp(ydata1(:,:,1:2))
% z_3d_uw=imresize(squeeze(z_uw(:,:,floor((sl-1)/48*72)+1,:)),[128,128]);
% z_3d_uw = fliplr( rot90(z_3d_uw));
% ydata1 = rm_linear_ph(ydata1,1,-1,pi/2);

% ydata2 = rm_linear_ph(ydata2,1,-1,pi/2);
%     ydata = [angle(ydata1.*conj(ydata2)),angle(ydata1),angle(ydata2)] ;

ydata12 = rm_linear_ph(ydata1.*ydata2,1,3,-pi-1);
ph_small = angle(ydata12);

% imdisp(angle([ydata1,ydata2,ydata12,ydata1.*conj(ydata2)]))

%
% ph_mean =  repmat( mean(ph_small,3), [1,1,size(ph_small,3)]);
sz = size(ph_small);

%ph_small_2 = 11*(ph_small - ph_mean);
ph_small_uw = unwrap(ph_small,[],3);
ph_small_fft = ifft(ph_small_uw,[],3);
% ph_small_fft(:,:,[1,3:end-1])=0;
ph_small_fft(:,:,1)=0;
phase_shift = 0.1999;

% phase_shift = -2.8403;

% phase_shift = 0.3013;
% phase_shift = 0.12;

% ph_small_2 = calc_small_ph(cimgs,sl,phase_shift)*10.1/2;
 ph_small_2 = calc_small_ph2(cimgs,sl,phase_shift)*10.1;

%%


 ph_big= angle(ydata1.*conj(ydata2));
 ph_big1= angle(rm_linear_ph(ydata1,1,1,0));
 ph_big2= angle(rm_linear_ph(ydata2,1,1,0));


% imdisp([ph_big ph_big1 ph_big2 ph_small_2],[-pi,pi])
%
% ph_big= angle(ydata1.*conj(ydata2).*exp(1i.*pi));

% ph_big = ph_big(90:110,40:end-40,:);

% ph_small_2 = ph_small_2(90:110,40:end-40,:);


N=2;
fft_points = 8;
v = [-N:N];
for k = 1:fft_points
    in{k} =v;% {v,v,v,v};
end
uw_kernal = combvec(in{:}).*2*pi;

% my_ph_uw = myunwrap(ph_big,ph_small_2,1,uw_kernal);

my_ph_uw12 = myunwrap(ph_big,ph_small_2,1);
% my_ph_uw1 = myunwrap(ph_big1,ph_small_2,0.5);
% my_ph_uw2 = myunwrap(-ph_big2,ph_small_2,0.5);
%  z_3d_uw2 = myunwrap(z_3d_uw,ph_small_2,1);
% z_3d_uw2 = unwrap(z_3d_uw,[],3);
% z_3d_uw2 = z_3d_uw;
%

%%
my_ph_uw = remove_2pi_dc(my_ph_uw12);

%
my_ph_uw_tmp = my_ph_uw;
sz = size(my_ph_uw);
threshold=0.1;
min_val=ones(sz(1:2))*10000;
for dc =-2*pi:pi/64:2.*pi

    my_ph_uw_fft = ifft(my_ph_uw_tmp,[],3);
    m3=sum(abs(my_ph_uw_fft(:,:,3:end-1)).^2,3);
    mask  = m3 < min_val*0.96;
    
    min_val(mask) = m3(mask);

    mask2 = repmat(mask,[1,1,size(my_ph_uw_fft,3)]);

    my_ph_uw(mask2) = my_ph_uw_tmp(mask2);
    
    my_ph_uw_tmp = myunwrap(ph_big,ph_small_2+dc,1);
end
%
my_ph_uw=remove_2pi_dc(my_ph_uw);





%%
% 
% my_ph_uw_fft = ifft(my_ph_uw2,[],3);
% 
% sc = sign(my_ph_uw_fft(:,:,[1])) .* round(abs(my_ph_uw_fft(:,:,1))/2/pi);
%  my_ph_uw_fft(:,:,[1]) =  my_ph_uw_fft(:,:,[1]) - sc*2*pi;
% %   my_ph_uw_fft(:,:,1) = 0;
% my_ph_uw = fft(my_ph_uw_fft,[],3);
% % my_ph_uw2=my_ph_uw;

% my_ph_uw = remove_2pi_dc(my_ph_uw12);
my_ph_uw1 = remove_2pi_dc(my_ph_uw1);
my_ph_uw2 = remove_2pi_dc(my_ph_uw2);
try close(1000); end;
figure(1000);
plt = plot([1,2,3,4,5], ones(5));
legend({'big menc','my unwrap','pipeline unwrap','small menc'});
xlabel('time offsets');
% imdisp([ph_big,my_ph_uw, unwrap(my_ph_uw,[],3),ph_small,ph_small_uw,ph_small_2],[-2*pi 2*pi],'size',[nan,1])
% imdisp([ph_big,my_ph_uw1,my_ph_uw2,my_ph_uw12,my_ph_uw1+my_ph_uw2,ph_small_2,z_3d_uw2,unwrap(ph_big,[],3)],[-5*pi 5*pi],'size',[nan,1]);
imdisp([ph_big,my_ph_uw12,ph_small_2,my_ph_uw],[-5*pi 5*pi],'size',[nan,1]);

colormap(awave)
linestyle={'-','-.','--','--'};
%
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
    ydata(:,3) =squeeze( my_ph_uw(yi,xi,:)); % 
    ydata(:,4) =squeeze( ph_small_2(yi,xi,:)); % 
    xdata = [1:size(ydata,1)]';

    for k=1:4
    set(plt(k),'xdata',xdata,'ydata',ydata(:,k),'LineStyle',linestyle{k});
    end
    set(plt(1).Parent,'ylim',[-14,14]);
    set(plt(5),'xdata',[1 8],'ydata',[pi pi])
%       ydata_uw=squeeze(s3.x_ux(yi,xi,sl,:));
%     fft(ydata_uw)
%     set(plt2,'xdata',xdata,'ydata',ydata_uw)
%     
    
     ang = angle(ifft(ydata))
     amp = abs(ifft(ydata))
    
end












