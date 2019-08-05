function ov2(img,voxel_size,clim,mask,func,varargin)
if nargin <2
    voxel_size=[];
end
if nargin <3
    clim = [];
end
if nargin<4
    mask = [];
end
if nargin<5
    func = [];
end
if ~iscell(img)
    img={img};
end

%
% warning('flip and rot nii cordinate image for display')
% for k=1:numel(img)
%     img{k} = rot90(flip(img{k},1));
% end
% if ~isempty(mask)
%    mask = rot90(flip(mask,1));
% end
%
for k=1:numel(img)
    if ~isempty(mask)
        sz1 = size(img{k});
        sz2 = size(mask);
        if any(sz1(1:3) ~= sz2(1:3))
            mask2 = imresize_nd(mask,sz1(1:3));
            img{k} = bsxfun(@times,img{k},mask2);
        else
            mask = single(mask);
%             mask(mask==0)=nan;
            
            img{k} = bsxfun(@times,img{k},mask);
        end
    end
    if ~isempty(func)
        img{k} = func(img{k},varargin{:});
    end
end
ov(img,clim,voxel_size)


