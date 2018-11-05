
Nt = 16;
mask = repmat( magn > 400,[1 1 1 Nt]);
phs_x_uw2 = first_harmonic(phs_x_uw(:,:,:,:),Nt).*mask;
phs_y_uw2 = first_harmonic(phs_y_uw(:,:,:,:),Nt).*mask;
phs_z_uw2 = first_harmonic(phs_z_uw(:,:,:,:),Nt).*mask;



% [Nx Ny Nz Nt] = size(phs_x_uw2);
% x = 1:Nx;
% y = 1:Ny;
% z = 1:Nz;
% [X Y Z] = meshgrid(x,y,z);
% mask = repmat( magn > 0,[1 1 1 Nt]);
% 
% k=1;
% X2 = bsxfun(@plus,  X,k*mask.*phs_x_uw2);
% Y2 =  bsxfun(@plus, Y,k*mask.*phs_y_uw2);
% Z2 =  bsxfun(@plus, Z,k*mask.*phs_z_uw2);
% X2 = double(X2) ;
% Y2 = double(Y2);
% Z2 = double(Z2);

mask = repmat( magn > 2000,[1 1 1 Nt]);
%%
d=0;
for t=1:Nt;
    
    msk = mask>2;
    msk(:,:,:,t) = mask(:,:,:,t);
    ux = phs_x_uw2(:,:,:,t) - d*mean(phs_x_uw2(msk)) ;
    uy = phs_y_uw2(:,:,:,t) - d*mean(phs_y_uw2(msk));
    uz = phs_z_uw2(:,:,:,t) - d*mean(phs_z_uw2(msk));
    
    D = cat(4,ux,uy,uz);
    M(:,:,:,t) = imwarp(magn,D);
end