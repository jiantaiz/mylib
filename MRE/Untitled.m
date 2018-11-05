sli =[6];
Mx =cat(2,imresize(s6.x_ux(:,:,sli,:),2),s10.x_ux(:,:,sli,:),s10.x_uw(:,:,sli,:)); 
Mx = permute(Mx,[1,2,4,3]);

My =cat(2,imresize(s6.y_ux(:,:,sli,:),2),s10.y_ux(:,:,sli,:),s10.y_uw(:,:,sli,:)); 
My = permute(My,[1,2,4,3]);

Mz =cat(2,imresize(s6.z_ux(:,:,sli,:),2),s10.z_ux(:,:,sli,:),s10.z_uw(:,:,sli,:)); 
Mz = permute(Mz,[1,2,4,3]);

M = cat(1,Mx,My,Mz);
imdisp(M(:,:,:))
%%
sli =[6];
Mx =cat(2,imresize(s6.x_ux(:,:,sli,:),2),s9.x_uw(:,:,sli,:)); 
Mx = permute(Mx,[1,2,4,3]);

My =cat(2,imresize(s6.y_ux(:,:,sli,:),2),s9.y_uw(:,:,sli,:)); 
My = permute(My,[1,2,4,3]);

Mz =cat(2,imresize(s6.z_ux(:,:,sli,:),2),s9.z_uw(:,:,sli,:)); 
Mz = permute(Mz,[1,2,4,3]);

M = cat(1,Mx,My,Mz);
imdisp(M(:,:,:))
%%
M=cat(2,s6.cimgs(:,:,sli,1,1),s10.cimgs(1:2:end,1:2:end,sli,1,1));
imdisp(abs(M))
%%
M=cat(2,s6.x_ux(:,:,sli,1,1),s10.x_ux(1:2:end,1:2:end,sli,1,1));
imdisp((M))
