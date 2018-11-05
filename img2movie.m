function img2movie(folder,output_filename,varargin)
% MAT2GIF converts a 3D matrix to an animated gif file
%    MAT2GIF(A,output_filename,varargin)
%        Options:
%        'cmap': colormap, gray(256) (default)| 256x3 colormap 
%        'clim': color limits, [min max](default) | two-element vector
%        'DelayTime': Delay before displaying next image, 0.5(default) | scalar 
%        'LoopCount': Number of times to repeat animation, Inf(default) | integer
%    Example:
%    img2movie('./','test.avi','DelayTime',0.2)
%
%    See also: 

% AUTHOR    : Yi Sui
% DATE      : 04/01/2017

p = inputParser;
p.addRequired('folder',@isstr);
p.addRequired('output_filename',@isstr);
p.addParameter('cmap',gray(256));
p.addParameter('DelayTime',0.5);
p.addParameter('LoopCount',Inf);
p.addParameter('clim',[]);
p.parse(folder,output_filename,varargin{:});

DelayTime = p.Results.DelayTime;
clim = p.Results.clim;
cmap = p.Results.cmap;
LoopCount = p.Results.LoopCount;

writerObj = VideoWriter(output_filename);
writerObj.FrameRate = 1./DelayTime;
open(writerObj);

IMs = dir(folder);
IMs = IMs(~[IMs.isdir]);

for k = 1:numel(IMs)
    %progress bar
    IM = imread(fullfile(folder,IMs(k).name));
    f = im2frame(IM);
    writeVideo(writerObj,f);
end
close(writerObj);
return
