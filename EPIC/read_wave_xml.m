function wv =  read_wave_xml(xmlfile,fig,M)
xml = xmlread(xmlfile);
elem = xml.getElementsByTagName('sequencer');
N = elem.getLength;
for k=1:N
    data = elem.item(k-1).getTextContent;
    wv{k} = str2num(data);
    sequencer{k} = char(elem.item(k-1).getAttributes.item(1).getNodeValue);
end

%%
if fig == 1
    figure;
    
    if numel(M) == 1
        M=[1:M];
    end
    
    for k=1:numel(M);
        ax(k)=subplot(numel(M),1,k);
        x=wv{M(k)};plot(x(:,1),x(:,2));
        title(sequencer{M(k)})
    end
    linkaxes(ax,'x')
end