function ph_bmf = remove_bulk_motion(ph_uw,mask,method,varargin)
% REMOVE_BULK_MOTION removes MRE bulk motion
%    ph_bmf = REMOVE_BULK_MOTION(ph_uw,mask,method,varargin) ...
%
%    Example:
%    ... 
%
%    Subfunctions: 
%    See also: 

% AUTHOR    : Yi Sui
% DATE      : 05/16/2017
%%
%ph_bmf = remove_bulk_motion(ph_uw,mask,method,vargin)
%
% ph_uw: unwrapped phase map
% mask: binaray mask map


[x y n] = size(ph_uw);


% switch method
%     case 'average'
%     case 'gaussian'
%
% end

ph_bmf = ph_uw; % init ph_bmf to ph_uw;

for k = 1:n
    ph = ph_uw(:,:,k);
    msk = mask(:,:,k);
    %     sigma = 0.5;
    switch method
        case 'average'
            ph(isnan(ph)) = 0;
            ph(~msk) = 0;
            m = mean(ph(ph~=0));
            ph_bmf(:,:,k) = ph_bmf(:,:,k) - m ;
        case 'gaussian'
            %             if numel(varargin) > 0
            %                 sigma = varargin(1);
            %             else
            %                 sigma = 0.5;
            %             end
            m = imgaussfilt(ph,varargin{:},'FilterDomain','frequency');
            ph_bmf(:,:,k) = ph_bmf(:,:,k) - m ;
        case 'bhp'
            %             thresh = varargin(1);
            %             n = varargin(2);
            ph(isnan(ph)) = 0;
            ph(~msk) = 0;
            ph_fft = fftshift(fft2(ph));
            ph_fft = bhp(ph_fft,varargin{:});
            ph_bmf(:,:,k)  = (ifft2(ifftshift(ph_fft))).*msk;
        case 'blpf'
            ph(isnan(ph)) = 0;
            ph(~msk) = 0;
            ph_fft = fftshift(fft2(ph));
            ph_fft = blpf(ph_fft,varargin{:});
            ph_bmf(:,:,k)  = (ifft2(ifftshift(ph_fft))).*msk;
        case 'poly'
            if k ==1
                order =varargin{1};
                [X, Y] = meshgrid(1:x,1:y);
                 A = poly_nd(order,X(:),Y(:));
                 %[ones(numel(X), 1), X, Y, X.^2, Y.^2, X.*Y];
            end
%             X = X(msk);
%             Y = Y(msk);
            s = ph(msk);
            A2 = A(msk(:), :);

            coeff = pinv(A2)*s;
            ph(msk) = s - A2*coeff;            
            ph_bmf(:,:,k) = ph.*msk;
            
        otherwise
            error('The filter method "%s" not defined',method);
    end
    
end
