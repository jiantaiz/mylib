function [IMC, IMw, IM, EPI, SE] = PSF_correction(PSF_data, MRE_data,pe_dir, acc)
if nargin<3
    pe_dir=1;
end
if nargin<4
    acc=4;
end
ro_dir=mod(pe_dir,2)+1;
%PSF_data: [nx ny nz ns]
% acc = 4;
% K = permute(PSF_data,[ro_dir pe_dir 4 3]);
K = permute(PSF_data,[ro_dir pe_dir 3 4]);
%K : [nx ny ns nz]
% K = flip(K,3);
[nx ny ns nz] = size(K);
w = hanning(ns);
w = reshape(w,[1 1 ns 1]);
w=1
K= bsxfun(@times, K , w);
IM=K;
dim = 3;
imsize = floor(nx/acc);
IM(:,:,end+1:imsize,:)=0; %zero padding;
IM = circshift(IM,[0 0 ceil((imsize-ns)/2) 0]);
IM = ifftshift(ifft(fftshift(IM,dim),[],dim),dim);
IM = repmat(IM, [1 1 acc 1]);
% IM = circshift(IM,[0, 0, floor(size(IM,3)/(2*acc)), 0]);

IM = flip(IM,3);
%% water mask
[nx ny ns nz] = size(IM);
mask = zeros([ny ns]);
a=floor(ns/acc/2*1);%mask width
for k=1:ny
    idx = (round(k/ny*ns)-a) :(round(k/ny*ns)+a) ;
    %     idx = (round((ny-k+1)/ny*ns)-a) :(round((ny-k+1)/ny*ns)+a) ;
    idx = mod( idx-1 + ns, ns )+1;
    %     mask(k,idx)=gausswin(2*a+1,5);
    %     mask(k,idx)=gausswin(2*a+1,4)./max(gausswin(2*a+1,4));
    mask(k,idx)=1;
end
IMw = bsxfun(@times, IM, reshape(mask, [1 ny ns 1]));
% [IMwmax idx_m] = max(abs(IMw),[],2);
% mask2 = bsxfun(@eq, abs(IMw), IMwmax);
%% correction
if ~isempty(MRE_data)
    cimgs = MRE_data;
    IM_raw= permute(cimgs,[ro_dir pe_dir 3 4 5]);
    % IM_raw = flip(IM_raw,2);
    [nx ny nz nt]= size(IM_raw);
    sz = size(IM_raw);
    sz(2) =ns;
    IMC = zeros(sz);
    for t=1:nt
        t
        IM_y = squeeze(IM_raw(:,:,:,t));%distorted
        IM_c = zeros([nx, ns, nz]);
        %     win= gausswin(2*a+1,2);
        win=1;
        parfor sl =1:nz
            for x = 1:nx
                tmp=0;
                for y=1:ny
                    p0 = IMw(x,y,:,sl);
                    p=p0(:);
                    %p^2 gives more weighting on PSF peak, image is sharper
                    %             IM_c(x,:,sl) = IM_c(x,:,sl) + (IM_y(x,y,sl)) .*  abs(p).^2 ./ sqrt(sum(abs(p).^2)) ;
                    %  IM_y^2 gives more weighting on high intensity pixels, less noise is spread to other pixesl, this increases SNR
                    %             IM_c(x,:,sl) = IM_c(x,:,sl) + (IM_y(x,y,sl))*abs(IM_y(x,y,sl)) .*  abs(p) ./ (sum(abs(p))) ;
                    %
                    % both IM_y^2 and p^2
                    if abs(IM_y(x,y,sl))>10 && max(abs(p))>0
                        tmp = tmp + (IM_y(x,y,sl))*abs(IM_y(x,y,sl)) .*  abs(p).^2 ./ (sum(abs(p).^2)) ;
                        %                     tmp = tmp + (IM_y(x,y,sl))*abs(IM_y(x,y,sl)) .*  abs(p)./ (sum(abs(p))) ;
                        %                     tmp = tmp + (IM_y(x,y,sl)) .*  abs(p)./ (sum(abs(p))) ;
                    end
                    
                end
                IM_c(x,:,sl) = tmp;
            end
            
        end
        IMC(:,:,:,t)= sqrt(abs(IM_c)).*exp(1i*angle(IM_c));
    end
    IMC = permute(IMC,[ro_dir pe_dir 3 4 5]);
    
else
    IMC=[];
end
IMw2 = abs(IMw(:,:,:,:)).^2 .*exp(1i*angle(IMw(:,:,:,:)));

if nargout>3
%     EPI = squeeze(sqrt(sum(abs(IMw(:,:,:,:)).^2,3))); %EPI
    EPI = squeeze(sum(IMw2,3));
    EPI = sqrt(abs(EPI)).*exp(1i*angle(EPI));
    EPI = permute(EPI,[ro_dir,pe_dir,3]);
end
if nargout>4
%     SE =  squeeze(sqrt(sum(abs(IMw(:,:,:,:)).^2,2)));
    SE = squeeze(sum(IMw2,2));
    SE = sqrt(abs(SE)).*exp(1i*angle(SE));
    SE = permute(SE,[ro_dir,pe_dir,3]);
end



