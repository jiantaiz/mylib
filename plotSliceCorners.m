function plotSliceCorners(corners, varargin)
% PLOTSLICECORNERS ...
%    PLOTSLICECORNERS(corners) ...
%    corners: 3 x 3 x numslices, UpperLeft = corners(:,1,:), UpperRight = corners(:,2,:), LowerLeft = corners(:,3,:); 
%    Example:
%    ... 
%
%    Subfunctions: 
%    See also: 

% AUTHOR    : Yi Sui
% DATE      : 05/04/2018
%%
sz = size(corners);
if sz(1) == 9
    cor = reshape(corners,3,3,sz(end));
else
    cor=corners;
end
cor(:,4,:) = cor(:,2,:)+cor(:,3,:)-cor(:,1,:);
cor = cor(:,[1 2 4 3 1],:);% reorder points
plot3( squeeze(cor(1,:,:)),squeeze(cor(2,:,:)),squeeze(cor(3,:,:)),varargin{:});
hold on;
ULmarkers = plot3( squeeze(cor(1,1,:)),squeeze(cor(2,1,:)),squeeze(cor(3,1,:)),'x',varargin{:});
URmarkers = plot3( squeeze(cor(1,2,:)),squeeze(cor(2,2,:)),squeeze(cor(3,2,:)),'o',varargin{:});
LLmarkers = plot3( squeeze(cor(1,4,:)),squeeze(cor(2,4,:)),squeeze(cor(3,4,:)),'d',varargin{:});
legend([ULmarkers,URmarkers,LLmarkers],{'UpperLeft','UpperRight','LowerLeft'});
axis equal
xlabel('x (mm)');
ylabel('y (mm)');
zlabel('z (mm)');
grid on;
box on;
