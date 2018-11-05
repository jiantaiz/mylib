function [ph_uw, Ds, jmps]= myunwrap_multi_menc(ph, menc,N,dim,weight,flag_qual_guide,magn,mask,flag_plot)
% MYUNWRAP_MULTI_MENC unwraps phases with multi-menc data
%    ph_uw = MYUNWRAP_MULTI_MENC(ph, menc,N,dim,weight) ...
%        ph: 6-dimention matrix ph(x,y,z,dir,timeoffsets,menc)
%        menc: ascending with the the first one normalized to 1. e.g. menc = [1 2.6 2.7 4.3]
%            the menc is supposed to be a non-integer number for this algrithm to work
%        N: max 2pi jump, +-N*2pi.
%        dim: specifies the dimension of menc
%        weight: weightings of each menc
%
%    See also: unwrap_dual_menc

% AUTHOR    : Yi Sui
% DATE      : 03/24/2017


if nargin < 4
    dim = 6;
end
if nargin < 5
    weight = ones(1,numel(menc));
end
if nargin < 6
    flag_qual_guide=0;
end
if nargin <9
    flag_plot = 0;
end
[phs perm,nshifts]= shiftdata(ph, dim);

nmenc = numel(menc);
npoints = numel(phs(1,:));
d=zeros(nmenc-1, npoints);
D =zeros(2*N+1,npoints);
for k=-N:N
    
    x=phs;
    x(end,:) = phs(end,:) +2*pi*k;
    
    %I'm assuming the smallest phase should be less than 2pi. Not really needed.
    %     ix1 = x(end,:) > 2*pi*menc(end)/menc(1);
    %     ix2 = x(end,:) < -2*pi*menc(end)/menc(1);
    %
    %     x(end,ix1&ix2) = 0;
    
    for m = 1:numel(menc)-1;
        x(m,:) = x(end,:) * (menc(m) ./ menc(end));
        %         x(m,:) = angle(exp(1i.*x(m,:)));
        %         d(m,:) = abs(x(m,:) - phs(m,:)); %not working well, should do it in complex domain
        %         d(m,:) = abs(exp(1i.*x(m,:)) - exp(1i.*phs(m,:))  );
        %            d(m,:) = abs(angle( exp( 1i.*x(m,:) - 1i.*phs(m,:) )  )); %no difference with above
        %           d(m,:) = weight(m).* abs(exp(1i.*x(m,:)) - exp(1i.*phs(m,:))  ).^2;
        %             d(m,:) = weight(m).* angle( exp(1i.*x(m,:) - 1i.*phs(m,:))  ).^2 .*(menc(m) ./ menc(end)).^2;
        %             d(m,:) = weight(m).* abs(angle( exp(1i.*x(m,:) - 1i.*phs(m,:))  )) .* (menc(m) ./ menc(end));
%         d(m,:) = weight(m).* abs(angle( exp(1i.*x(m,:) - 1i.*phs(m,:))))  .* (1+ abs(x(m,:))/10  );
        d(m,:) = 1./menc(m)*weight(m).* abs(angle( exp(1i.*x(m,:) - 1i.*phs(m,:))));
        
        
    end
%     rm = (abs(x(1,:)) >pi);
    
    D(k+N+1,:) = sum(d,1);
%     D(k+N+1,rm) = inf; 
    %    D(k+N+1,:) = sum(d,1).* (abs(k)+1);
    %     D(k+N+1,:) = sqrt(abs(k)+1)*sum(d,1);
    %     D(k+N+1,:) = (10+abs(x(end,:)).^1).*sum(d,1);
    %     D(k+N+1,:) = (abs(x(end,:)).^1).*sum(d,1);
    
    %      D(k+N+1,:) = (abs(k)+1)*20000;
    
end
% D

% reorder D so that the jmp is [0 +1 -1 +2 -2 ... +N -N]*2pi. if +1 and -N
% both have minimal error, the one with smaller abs (i.e. +1) will be
% chosen.

neworder = [N+2:2*N+1; N:-1:1];
neworder = [N+1;neworder(:)];
D = D(neworder,:);
[M,idx] = min(D,[],1);
jmp = neworder(idx)'-1-N;
% ph_uw = phs(end,:)+2*pi*jmp;
phs(end,:) = phs(end,:)+2*pi*jmp;


[Ds,idx2] = sort(D,1);
sz = size(phs);
sz(1) = size(Ds,1);
Ds  = reshape(Ds, sz);
Ds = unshiftdata(Ds, perm, nshifts);
idx2 = reshape(idx2,sz);
idx2  = unshiftdata(idx2, perm, nshifts);
jmps = neworder(idx2)-1-N;
%
% jmp1 = idx2(1,:)-1-N;
% jmp2 = idx2(2,:)-1-N;
% ph_uw1 = phs(end,:)+2*pi*jmp1;
% ph_uw2 = phs(end,:)+2*pi*jmp2;
%
% Ds1 = Ds(1,:);
% Ds2 = Ds(2,:);
% % DDs = (Ds2-Ds1) < 1 ;
%
% c1=abs(ph_uw1);
% c2=abs(ph_uw2);
%
% keep1 = (Ds2-Ds1 >=0.5) | (Ds2-Ds1 < 0.5 & c1 <= c2 );
% keep2 = ~keep1;
%
% % keep1 = 1;
% % keep2 = 0;
% % [~, idx3] = min(abs(ph_uw),[],1);
%
%
% phs(end,:) =   ph_uw1.* keep1 + ph_uw2 .* keep2;

% ph_uw = unshiftdata(ph_uw, perm, nshifts);
phs = unshiftdata(phs, perm, nshifts);
ph_uw = phs(:,:,:,:,:,end);



if flag_qual_guide ==1
    [nx, ny, nz, nt,nd] = size(ph_uw);
    
    if ~exist('magn','var');
        magn = ones([nx,ny,nz]);
    end
    if ~exist('mask','var');
        mask=ones([nx,ny,nz]);
    end
    
    
    for t = 1:nt;
        for d = 1:nd
            for sl = 1:nz
                ph = ph_uw(:,:,sl,t,d);
                % mask = mag > 200;
                
                mag = magn(:,:,sl);
                msk = mask(:,:,sl);
                pickstart = 0;
                % im_phase_quality = Ds(:,:,sl,t,1,1);
                im_phase_quality=Ds(:,:,sl,t,d,1)./Ds(:,:,sl,t,d,2);
                %                 im_phase_quality=[];
                jump = squeeze( jmps(:,:,sl,t,d,1:4));
                jump=bsxfun(@minus,jump,jump(:,:,1));
                % [ph_uw, ph_qual]=unwrap2D(ph,mask,pickstart,im_phase_quality,mag);
                
                [uw, ph_qual]=unwrap2D_sui(ph,msk,pickstart,im_phase_quality,mag,jump,flag_plot);
                ph_uw(:,:,sl,t,d) =uw ;
            end
        end
    end
end










