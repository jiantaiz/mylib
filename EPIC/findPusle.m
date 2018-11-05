function findPusle(varargin)
if numel(varargin)<2
    pulsename = varargin{1};
    ax = gca;
else
    h=varargin{1};
    pulsename = h.String{h.Value};
    ax = h.UserData;
end


% get(ax);
if strcmpi(pulsename, '0th moment')
    plot_moment(ax,0)
    return;
end
if strcmpi(pulsename, '1st moment')
    plot_moment(ax,1)
    return;
end

obj = findobj(ax,'Tag',pulsename);
if isempty(obj)
    fprintf('Pulse: %s is not found', pulsename);
    return;
end

cursorLine = findobj(ax,'Tag', 'cursorLine');

x = mean(obj(1).XData);
x= [x x];
y = get(ax,'YLim');

if isempty(cursorLine)
    
    line( x,y,'Color',[0 0 0], 'Linewidth',1, 'Tag', 'cursorLine','parent',ax);
else
    cursorLine.XData = x;
    cursorLine.YData = y;
end
