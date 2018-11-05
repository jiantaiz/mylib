function addPulseList(handler)
% ADDPUSLELIST 
%    ADDPUSLELIST(handler) ...
%
%    Example:
%    ADDPUSLELIST(gca)
%
%    Subfunctions: 
%    See also: 

% AUTHOR    : Yi Sui
% DATE      : 06/20/2017
%%
lines = findobj(handler,'type','line');
tags = get(lines,'tag');
tags = unique(tags);

tags = [{'0th moment'; '1st moment'} ; tags];
% 
pos = get(handler,'Position');

pos(1) = pos(1)+pos(3)+0.004;
pos(3) = pos(3)/9;
h = uicontrol('units','norm','style','listbox',...
  'string',tags,'position',pos,'Callback',@findPusle,'UserData',handler); %a listbox
