function plot_moment(ax,m)
if nargin<1
    ax = gca;
end
if nargin<2
    m = 0;
end
tag = sprintf('moment%d',m);
obj = findobj(ax,'Tag',tag);
if ~isempty(obj)
    delete(obj);
    return;
end


rf1 = findobj(gcf,'Tag', 'gzrf1');
% rf2 = findobj(gcf,'Tag', 'gzrf2');
rf2l =findobj(gcf,'Tag', 'gzrf2l1d');
rf2r =findobj(gcf,'Tag', 'gzrf2r1a');
Gmax = 5 ;%Gs/cm

wave = findobj(ax,'Tag','waveform');

t = wave.XData*1000; %usec
G = wave.YData / 2^31 * Gmax;

t_rf1 =  mean(rf1.XData);
if isempty(rf2l)
    t_rf2 = t(end);
else
    t_rf2 =  (rf2l.XData(end) + rf2r.XData(1))/2;
end

[tt,GG,mm] = calc_moment([t' G'], m ,[t_rf1, t_rf2]);
hold(ax,'on')
plot(ax,tt,mm./max(mm(:))*2^30,'linewidth',2,'tag',sprintf('moment%d',m));

