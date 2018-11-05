numofslice = 8;
numoftimeoff = 8;
phs_x2 = angle(exp(1i*phs_x).* exp(1i*0));
phs_x2 = reshape(phs_x2,[128,128,numofslice,4,numoftimeoff]);
imdisp(phs_x2(:,:,1,:,1))

phs_y2 = angle(exp(1i*phs_y).* exp(1i*pi));
phs_y2 = reshape(phs_y2,[128,128,numofslice,4,numoftimeoff]);
imdisp(phs_y2(:,:,1,:,1))

phs_z2 = angle(exp(1i*phs_z).* exp(1i*0));
phs_z2 = reshape(phs_z2,[128,128,numofslice,4,numoftimeoff]);
imdisp(phs_z2(:,:,1,:,1))



phs_x_fft = ifft(phs_x2,[],4);
phs_y_fft = ifft(phs_y2,[],4);
phs_z_fft = ifft(phs_z2,[],4);

imdisp(angle(phs_x_fft(:,:,1,2,:)))
imdisp(angle(phs_y_fft(:,:,1,2,:)))
imdisp(angle(phs_z_fft(:,:,1,2,:)))
%%
sl=8;
f = 2;
% imdisp(mag)
img = cat(2, phs_x_fft(:,:,sl,f,:),phs_y_fft(:,:,sl,f,:),phs_z_fft(:,:,sl,f,:));
% imdisp(img)
% imdisp(angle(img))
% imdisp(angle(bsxfun(@times, (img(:,:,:,:,1:end-1)), (conj(img(:,:,:,:,end))))))
% imdisp(bsxfun(@minus, abs(img(:,:,:,:,1:end-1)), abs(img(:,:,:,:,end))))

imdisp(angle(img(:,:,:,:,1:end)))
