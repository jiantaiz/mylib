function mat2movie(A,output_filename,varargin)
% MAT2GIF converts a 3D matrix to an animated gif file
%    MAT2GIF(A,output_filename,varargin)
%        Options:
%        'cmap': colormap, gray(256) (default)| 256x3 colormap 
%        'clim': color limits, [min max](default) | two-element vector
%        'DelayTime': Delay before displaying next image, 0.5(default) | scalar 
%        'LoopCount': Number of times to repeat animation, Inf(default) | integer
%    Example:
%    mat2movie(A,'test.gif','clim',[0 7],'cmap',hsv(256),'DelayTime',1)
%
%    See also: 

% AUTHOR    : Yi Sui
% DATE      : 04/01/2017

p = inputParser;
p.addRequired('A',@isnumeric);
p.addRequired('output_filename',@isstr);
p.addParameter('cmap',gray(256));
p.addParameter('DelayTime',0.5);
p.addParameter('LoopCount',Inf);
p.addParameter('clim',[]);
p.parse(A,output_filename,varargin{:});

DelayTime = p.Results.DelayTime;
clim = p.Results.clim;
cmap = p.Results.cmap;
LoopCount = p.Results.LoopCount;

if isempty(clim)
    clim = minmax(A(:)');
end
[x, y, n] = size(A);
A = double(A);
writerObj = VideoWriter(output_filename);
writerObj.FrameRate = 1./DelayTime;
open(writerObj);
for k = 1:n
    %progress bar
    IM = A(:,:,k);
    IM(3,1:n) = clim(2);
    IM(1:5,k) = clim(2); 
    IM = mat2gray(IM, clim)*length(cmap);
    IM(IM<1)=1;
    f = im2frame(IM,cmap);
    writeVideo(writerObj,f);
end
close(writerObj);
return
