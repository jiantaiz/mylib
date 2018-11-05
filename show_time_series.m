function show_time_series(fig,ylim,xdata,Ntile)
% show_time_series(fig,ylim)
if nargin<4
    Ntile = 1;
end
    try close(1000); end;
    state = get(fig,'UserData');
    I = state.I;
    
    
%     figure(1000);
%     plot([1 4], [pi -pi; pi -pi],'k--');
%     grid on;
%     hold on
%     plt = plot([1,2,3,4,5], ones(5));
%     legend({'pi','-pi','big menc','my unwrap','small menc'});
%     xlabel('time offsets');
%     colormap(awave)
    figure(1000);
%     plt = plot([1 2],[1 2]);
    plt = plot(1:Ntile,ones(Ntile));
    grid on;
    linestyle={'-','-.','--','--'};
%     fig = gcf;
%     assignin('base','fig_9465748',fig)
sz=size(I);
    while 1
%         figure(fig);

%         [xi, yi, but] = evalin('base','myginput(fig_9465748, 1);');      % get a point
%         [xi, yi, but] = myginput(fig, 1)      % get a point
        figure(fig);
        [xi, yi, but] = ginput(1);
        if ~isequal(but, 1)             % stop if not button 1
            break
        end
        xi = round(xi); 
        yi = round(yi);
        [xi yi]

        for n=1:Ntile
            xx = xi + sz(2)/Ntile*(n-1);
            xx = mod(xx-1,sz(2))+1;
            ydata(:,n) = squeeze( I(yi,xx,:));% 
        end
        if nargin<3 ||isempty(xdata)
            xdata = [1:size(ydata,1)]';
        end
        if isempty(plt)
%             figure(1000);
%             plt = plot(xdata,ydata);
        else
            for k=1:Ntile
            set(plt(k),'xdata',xdata,'ydata',ydata(:,k),'LineStyle',linestyle{k});
            end
            set(plt(1).Parent,'ylim',ylim);
        end
  

%          ang = angle(ifft(ydata))
%          mag = abs(ifft(ydata))
    polarPrint( fft(ydata,[],1),'deg','.5f');
    disp(' ');
    end