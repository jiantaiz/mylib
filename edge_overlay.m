function IM = edge_overlay(I,magn,mask,method)
% IM = EDGE_OVERLAY(I,magn)  takes edge of magn and overlay it on I
% IM = EDGE_OVERLAY(I,magn,mask,method) masks the image
% method = 'Canny' 'log' 'Prewitt' 'Roberts' 'Sobel' or 'zerocross'
% by Yi Sui
% 4/1/2017

if nargin<3
    mask = [];
end
if nargin<4
    method = 'Sobel';
end
for k = 1:size(magn,3)
    myedges(:,:,k) = edge(magn(:,:,k),method);
end
mx = max(I(:));
IM = bsxfun(@plus, I, myedges.*mx*5);
if ~isempty(mask)
    IM = bsxfun(@times, IM,mask);
end