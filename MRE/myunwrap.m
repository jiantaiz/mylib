function uw = myunwrap(wph_big_menc,ph_small_menc, menc_ratio, uw_kernal,tol)
% uw = myunwrap(wph_big_menc,ph_small_menc, menc_ratio, uw_kernal)
% The motion is encoded with big and small menc. Big menc gives large and high
% SNR phase which however may be phase-wrapped. small menc gives small
% but no phase wrapp signal. We can use the small menc signal to get a
% rough idea of actural displacement and use it to guide the unwrapping of
% big menc phase signal.
%INPUT:
% wph_big_menc: phase signal from big menc
% ph_small_menc: phase signal from small menc
% ratio between big and small menc menc_ratio = big_menc / small_menc
%OUTPUT:
% uw: unwrapped phase

if nargin<3
    menc_ratio = 1;
end

if nargin <4
    uw_kernal = [];
end
if nargin <5
    tol=1;
end


if isempty(uw_kernal)
    
    ph_est = real(ph_small_menc * menc_ratio);
    
    idx = abs(ph_est) >0.*pi; 
    
    D = abs(ph_est - wph_big_menc)./(2.*pi.*tol);
    
     N = round(D);
%     N = floor(D);
    N = sign(ph_est-wph_big_menc).*N;
%     uw = wph_big_menc + N.*2*pi;
    uw = wph_big_menc;
    uw(idx) = wph_big_menc(idx) + N(idx).*2*pi;
    
%     uw_fft =ifft(uw);
%     
%     ph_est = ifft(uw);
%     
else
    %
 
    sz = size(wph_big_menc);
    uw = zeros(sz);
    M=size(uw_kernal,2);
    for i = 1:sz(1)
        i
        for j = 1:sz(2)
            yd1=wph_big_menc(i,j,:);
            yd2=ph_small_menc(i,j,:);
            y_uw = repmat(yd1(:),[1,M]) + uw_kernal ;
            
            delta = y_uw - repmat(yd2(:),[1,M]);
            
            delta_fft = ifft(delta,[],1);
            delta_fft(1,:) = 0;
            delta = fft(delta_fft,[],1);
            
%             y_uw_fft = ifft(y_uw,[],1);
%             y_uw_fft([1,3:end-1],:) = 0;
            
%             delta =  fft(y_uw_fft,[],1) - repmat(yd2(:),[1,M]);
            D = sum(delta.*conj(delta),1);
            [minD, idx] = min(D);
            uw(i,j,:) = y_uw(:,idx);
            
        end
    end
    
end



