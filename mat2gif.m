function mat2gif(A,output_filename,varargin)
% MAT2GIF converts a 3D matrix to an animated gif file
%    MAT2GIF(A,output_filename,varargin)
%        Options:
%        'cmap': colormap, gray(256) (default)| 256x3 colormap 
%        'clim': color limits, [min max](default) | two-element vector
%        'DelayTime': Delay before displaying next image, 0.5(default) | scalar 
%        'LoopCount': Number of times to repeat animation, Inf(default) | integer
%    Example:
%    mat2gif(A,'test.gif','clim',[0 7],'cmap',hsv(256),'DelayTime',1,'LoopCount',Inf)
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
p.addParameter('ProgressBar',0);
p.addParameter('Resize',1);
p.addParameter('Method','nearest');
p.addParameter('Mask',[]);
p.parse(A,output_filename,varargin{:});

DelayTime = p.Results.DelayTime;
clim = p.Results.clim;
cmap = p.Results.cmap;
LoopCount = p.Results.LoopCount;
ProgressBar = p.Results.ProgressBar;
Resize =  p.Results.Resize;
Method =  p.Results.Method;
Mask = p.Results.Mask;
if isempty(clim)
    clim = minmax(A(:)');
end
[x, y, n] = size(A);
A = double(A);
for k = 1:n
    
    IM = A(:,:,k);
    
    if ~isempty(Mask)
        IM(Mask) = clim(1);
    end
    IM(isnan(IM)) = clim(1);
    if ProgressBar == 1
    IM(3,1:n) = clim(2);
    IM(1:5,k) = clim(2);
    end

    A(:,:,k) = mat2gray(IM, clim)*256;
end
if Resize>1
    A = imresize(A,Resize,Method);
    [x, y, n] = size(A);
end

A = reshape(A, x, y, 1, n);
% Save as a gif
imwrite(A, cmap, output_filename, 'LoopCount', LoopCount, 'DelayTime', DelayTime);
return
