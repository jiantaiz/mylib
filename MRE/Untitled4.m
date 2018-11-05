sl = 10;
dir = 1;

figure(200);
axes;
plt = plot([1]);hold on;
plt2= plot([2 3 5]);
hold off;
grid on;

 ydata1 =squeeze( s3.cimgs(:,:,sl,:,dir));% positive meg
 ydata2 =squeeze( s3.cimgs(:,:,sl,:,dir+3)); % negtive meg
 ydata = angle(ydata1.*conj(ydata2));
 ydata = abs(ydata1);
 imdisp(ydata);
% imdisp(s3.cimgs(:,:,sl,1,1) );
while 1
    
    [xi, yi, but] = ginput(1);      % get a point
    if ~isequal(but, 1)             % stop if not button 1
        break
    end
    xi = round(xi); 
    yi = round(yi);
    
    ydata1 =squeeze( s3.cimgs(yi,xi,sl,:,dir));% positive meg
    ydata2 =squeeze( s3.cimgs(yi,xi,sl,:,dir+3)); % negtive meg
    
    ydata = angle(ydata1.*conj(ydata2));
    fft(ydata)
    
    xdata = 1:numel(ydata);
    set(plt,'xdata',xdata,'ydata',ydata);
    ydata_uw=squeeze(s3.x_ux(yi,xi,sl,:))
    fft(ydata_uw)
    set(plt2,'xdata',xdata,'ydata',ydata_uw)
    
    
     
    
end