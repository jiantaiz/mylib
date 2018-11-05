N=2;
fft_points = 8;


v = [-N:N];
clear in;
for k = 1:fft_points
    in{k} =v;% {v,v,v,v};
end
uw = combvec(in{:}).*2*pi;
%%
cmb_fft = ifft(uw,[],1);
cmb_fft_1 = cmb_fft(2,:);
cmb_fft_1 = unique(cmb_fft_1);
%%
figure;plot(real(cmb_fft_1),'.');hold on; plot(imag(cmb_fft_1),'.');
 


%%
yd2_fft = ifft(yd(:,2));
M=size(uw,2);
% for j = 1:size(uw,2);
%     delta = yd(:,1) + cmb(:,j) - yd(:,2) ;
%     delta_fft = ifft(delta);
    
    y_uw = repmat(yd(:,1),[1,M]) + uw; 
    y_uw_fft = ifft(y_uw,[],1);
    y_uw_fft(1,:) = 0;
%     y_uw2 = fft(y_uw_fft);
%     y_uw_mean = mean(y_uw);
    
%     delta = y_uw_fft(2) - yd2_fft(2) ;
    delta =  fft(y_uw_fft,[],1) - repmat(yd(:,2),[1,M]);
      
    D = sum(delta.*conj(delta),1);
% end
[minD, idx] = min(D);