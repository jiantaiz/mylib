function show_time_series2(fig,ylim,xdata)
% SHOW_TIME_SERIES2 shows the plot time series in a figure
%    SHOW_TIME_SERIES2(fig,ylim,xdata) specifies the figure handle, y-lim
%    and x-axis data.
%
%    Example:
%    show_time_series2(gcf,[-1 1],[1 2 3 4 5])
%    pick a point in current figure, and plot the time series at that
%    point. The current figure is assumed to have 5 images in this example.
%
%    See also: 

% AUTHOR    : Yi Sui
% DATE      : 05/16/2017
%%
% show_time_series(fig,ylim)
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
    plt = plot([1:9],magic(9),'-x');
    grid on;
    linestyle={'-','-.','--','--'};
%     fig = gcf;
%     assignin('base','fig_9465748',fig)
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

%         ydata(:,1) =squeeze( I(yi,xi,:));% 
        y=squeeze( I(yi,xi,:))% 

        for t=-4:4;
            ydata(:,t+5) = y + t.*2.*pi;
        end
        
        if nargin<3
            xdata = [1:size(ydata,1)]';
        end
        if isempty(plt)
%             figure(1000);
%             plt = plot(xdata,ydata);
        else
            for k=1:9

                set(plt(k),'xdata',xdata,'ydata',ydata(:,k),'LineStyle',linestyle{1});
            end
            set(plt(1).Parent,'ylim',ylim);
        end
  

%          ang = angle(ifft(ydata))
%          mag = abs(ifft(ydata))
    end
