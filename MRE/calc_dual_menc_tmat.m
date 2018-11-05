function [T,menc,menc_mtx]=calc_dual_menc_tmat(cp,meg,freq_motion,te,meg_num)
% CALC_DUAL_MENC_TMAT calculate menc, and transform matrix for dual menc
%    [T,menc,menc_mtx]=CALC_DUAL_MENC_TMAT(cp,meg,freq_motion,te) ...
%       INPUT:
%       cp: cornerpoint matrix (Nx4) or cornerpoint file
%       meg: relative amplitude of the negtive MEG ( 0.0-1.0 )
%       freq_motion: frequence of vibration
%       te: echo time in msec
%       OUTPUT:
%       T: transform matrix from low menc to high menc
%       menc: 3x1 vector, menc values in rad/um for x,y,and z
%       menc_mtx: 3x3x4 complex menc matrix at 4 different combination: pos+neg, neg, pos, and pos-neg 
%       
%    See also: calc_menc, unwrap_dual_menc

% AUTHOR    : Yi Sui
% DATE      : 07/05/2017
%%
if nargin < 5
    meg_num=1;
end
gz=1;%0.9 in psd, shows in cornerpoints already, so gz=1 here.
cme_gamp = [1,0,0;
    0,1,0;
    0,0,gz]; %positive amplitude for x-, y-, z-MEG

cme_gamp2 = cme_gamp.* ...
    [-meg, 0, 0;
    0, -meg, 0;
    0, 0, -1]; %negtive amplitude for x-, y-, z-MEG

if ischar(cp)
    cp = read_cp(cp);
end


[cp1,t180]=set_meg(cp, 1, meg_num);
t_rf=[t180-te/2,t180];

%caclulate menc matrix
for k = 1:size(cme_gamp,1)
    cp1=set_meg(cp, cme_gamp(k,:),meg_num);
    [menc_pos(k,:),phase_pos,menc_pos_cplx(k,:)] =calc_menc(cp1,t_rf,t180+te/2,freq_motion);
    
    cp2=set_meg(cp, cme_gamp2(k,:),meg_num);
    [menc_neg(k,:),phase_neg,menc_neg_cplx(k,:)] =calc_menc(cp2,t_rf,t180+te/2,freq_motion);
end

ang= angle(menc_pos_cplx(1)-menc_neg_cplx(1));
menc_pos_cplx = menc_pos_cplx.*exp(-1i*ang);
menc_neg_cplx = menc_neg_cplx.*exp(-1i*ang);
menc_add = menc_pos_cplx+menc_neg_cplx;
menc_sub = menc_pos_cplx-menc_neg_cplx;

menc_mtx = cat(3,menc_add,menc_neg_cplx,menc_pos_cplx,menc_sub);
menc = abs( diag(menc_mtx(:,:,end)));
T = menc_sub*(inv(menc_add));

if nargout<1
    disp('menc_pos:'); polarPrint(menc_pos_cplx,'deg','7.4f')
    disp('menc_neg:'); polarPrint(menc_neg_cplx,'deg','7.4f')
    disp('+:');polarPrint(menc_add,'deg','7.4f')
    disp('-:');polarPrint(menc_sub,'deg','7.4f')
    disp('Transform matrix:');polarPrint(T,'deg','7.4f');
    
    % animating PSD
    figure;
    
    for k = 1:size(cme_gamp,1)
        cp1=set_meg(cp, cme_gamp(k,:),meg_num);
        plot(cp1(:,1),cp1(:,2:end));
        hold on;
        plot([t_rf t180+te/2]*1000,[0 0 0],'x');
        hold off;
        legend({'X','Y','Z'})
        
        pause(0.4)
        cp2=set_meg(cp, cme_gamp2(k,:),meg_num);
        plot(cp2(:,1),cp2(:,2:end));
        hold on;
        plot([t_rf t180+te/2]*1000,[0 0 0],'x');
        hold off;
        
        legend({'X','Y','Z'})
        pause(.4)
    end
end
