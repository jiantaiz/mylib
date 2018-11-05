% for Yi's animation

% mat.mat should be a 3D array nx*ny*nz 
% load('mag.mat');
load('cimgs.mat');
mag = abs(cimgs(:,:,:,1));
mask = mag > 700;
mask = imopen(mask, strel('disk',3));
mask = imclose(mask, strel('disk',5));

nx=size(mag,1);ny=size(mag,2);nz=size(mag,3);
% slice index
slice_index = zeros(size(mag));%disp_digital (nx, ny, nz);
save(fullfile(pwd,'slice_index.mat'),'slice_index');

load('phs_3duw_bmf.mat');


 load(fullfile(pwd,'MRE_ZY','disp_absmu.mat'));
% mu_disp = 
if (exist(fullfile(pwd,'output_images'),'dir'))
else mkdir(fullfile(pwd,'output_images'));
end

big_array1=[mag;slice_index];
for i = 1:size(big_array1,3)
    %     imdisp(big_array1(:,:,11),[0 65536/10]);
    %     frame = getframe(gcf);
    %     im = frame2im(frame);
    %     [X,map] = rgb2ind(im,256);
    
    im = mat2gray(big_array1(:,:,i),[0 65536/10]);
    X = gray2ind(im, 256);
    map = gray(256);
    
    
    imwrite(X,map,fullfile(pwd,'output_images',['mag_s',int2str(i),'.png']));
    close;
end
close;
mask = repmat(mask,[1 1 1 size(x_uw,4)]);
big_array1=[x_uw.*mask,y_uw.*mask,z_uw.*mask];
for j = 1:size(big_array1,4)
    for i = 1:size(big_array1,3)
%         imdisp(big_array1(:,:,i,j),[-pi/2,pi/2],'map','awave'); %[-3.3 4.3] [-0.2 0.2]
%         frame = getframe(gcf);
%         im = frame2im(frame);
%         [X,map] = rgb2ind(im,256);
        im = mat2gray(big_array1(:,:,i,j),[-pi/2,pi/2]);
        X = gray2ind(im, 256);
        map = awave(256);
        if j == 1
            imwrite(X,map,fullfile(pwd,'output_images',['phs_s',int2str(i),'.gif']),'DelayTime',.5, 'LoopCount',3000); %3000, 0.1
            % delay time: default 0.1, 0.5 for 4 offsets
        else
            imwrite(X,map,fullfile(pwd,'output_images',['phs_s',int2str(i),'.gif']),'WriteMode','Append','DelayTime',.5);
        end
        close;
    end
end
close;

big_array3=[mu_disp.*mask(:,:,:,1)];
for i = 1:size(big_array3,3)
    %     imdisp(big_array3(:,:,i),[0,8],aaasmo); %[-3.3 4.3] [-0.2 0.2]
    %     axis equal
    % %     set(gca,'unit','pixels')
    % %     frame = getframe(gca,get(gca,'position'));
    %     im = frame2im(frame);g
    %     [X,map] = rgb2ind(im,256);
    
    
    im = mat2gray(big_array3(:,:,i),[0,8]);
    X = gray2ind(im, 256);
    map = aaasmo(256);
    imwrite(X,map,fullfile(pwd,'output_images',['abscmu_s',int2str(i),'.png']));
    close;
end
close;

% export to pptx
disp('export to pptx')
ppt_title = 'MRE_summary';
isOpen  = exportToPPTX();
if ~isempty(isOpen),
    % If PowerPoint already started, then close first and then open a new one
    exportToPPTX('close');
end
exportToPPTX('new','Dimensions',[10 7.5], ... %13.3 7.5
    'Title',ppt_title, ...
    'Author','MatLab', ...
    'HorizontalAlignment','center', ...
    'FontWeight','bold',...
    'FontSize',40);

for m = 1:size(mag,3)
    slideNum = exportToPPTX('addslide');
    exportToPPTX('addpicture',fullfile(pwd,'output_images',['mag_s',int2str(m),'.png']),'Scale','maxfixed','Position',[0,0,2.5,5]);
    %     exportToPPTX('addpicture',fullfile(filepath,['phs_amp_s',int2str(m),'.png']),'Scale','maxfixed','Position',[2.5,0,7.5,5]);
    exportToPPTX('addpicture',fullfile(pwd,'output_images',['phs_s',int2str(m),'.gif']),'Scale','maxfixed','Position',[2.5,0,7.5,2.5]);
%     exportToPPTX('addpicture',fullfile(filepath,['4phs_s',int2str(m),'.gif']),'Scale','maxfixed','Position',[3.75,1.88,5.63,5.63]);
     exportToPPTX('addpicture',fullfile(pwd,'output_images',['abscmu_s',int2str(m),'.png']),'Scale','maxfixed','Position',[2.5,2.5,2.5,2.5]);
    close;
end
newFile = exportToPPTX('save',fullfile(pwd,'output_images',ppt_title));
exportToPPTX('close');