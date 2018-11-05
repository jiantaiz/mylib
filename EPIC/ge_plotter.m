function ge_plotter(xmlfile,sequencers)
% GE_PLOTTER ...
%    GE_PLOTTER(xmlfile,sequencers) ...
%
%    Example:
%    ... 
%
%    Subfunctions: 
%    See also: 

% AUTHOR    : Yi Sui
% DATE      : 07/02/2018
%%
if nargin<2
    sequencers=[2:5];
end
xml = xmlread(xmlfile);
seqs = xml.getElementsByTagName('sequencer');

for k=1:seqs.getLength
    
   seq = seqs.item(k-1);
   attr = seq.getAttributes;
   da = seq.getChildNodes.item(1).getTextContent;
   
   Seq(k).data = str2num(da);
   Seq(k).id = char(attr.item(0).getValue);
   Seq(k).ttl = char(attr.item(1).getValue);
   Seq(k).xttl = char(attr.item(2).getValue);
   Seq(k).yttl = char(attr.item(3).getValue);
   
end

N = numel(sequencers);
ax = tight_subplot(N,1,0.006,0.1,0.1);
for k =1:N
%     ax(k)=subplot(N,1,k);
    s = sequencers(k);
    plot(ax(k),Seq(s).data(:,1),Seq(s).data(:,2));
%     title(Seq(s).ttl);
    ylabel(ax(k),Seq(s).ttl(9:end));
    ylim=get(ax(k),'YLim');
    set(ax(k),'YLim',ylim.*1.1);
end
% set(ax(2:2:end),'YAxisLocation','right');
set(ax(1:N-1),'XTickLabel',''); 
set(ax,'XGrid','on');
xlabel(ax(end),Seq(sequencers(end)).xttl);
linkaxes(ax,'x');
