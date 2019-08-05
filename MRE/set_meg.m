function [cp,t180] = set_meg(cp,meg_scale,num_meg_pairs)
% SET_MEG  set motion encoding gradient amplitude
%    cp = SET_MEG(cp,meg_scale) will find the location of the MEG in the given cornerpoint 'cp' and scale its amplitude with meg_scale
%    [cp,t180] = SET_MEG(cp,meg_scale) will also return the location of 180 RF pulse in ms
%    See also:

% AUTHOR    : Yi Sui
% DATE      : 07/05/2017
%%
if nargin<3
    num_meg_pairs=1;
end
m = max(cp(:,2));
m=1;

onesin =[0 m m 0 -m -m];
seq1 = [repmat(onesin,[1 num_meg_pairs]), 0 m m 0];
seq1_len = length(seq1);
% seq1 = [0 m m 0 -m -m 0 m m 0];% MEG sequence
seq1 = seq1./sqrt(sum(seq1.^2));
% seq2 = seq(end:-1:1);% reverse
seq2 = -seq1;

if numel(meg_scale) == 1
    meg_scale = ones(1,3).*meg_scale;
elseif numel(meg_scale)~=3
    error('meg_scale must be scaler or 1x3 vector');
end

for grad = 2:4 % x, y, z-gradient
    if grad==2
        g1 = conv(cp(:,grad), seq1,'same');
        g2 = conv(cp(:,grad), seq2,'same');
        gm1 = max(g1);
        gm2 = max(g2);
        gm2 = -100;
        g2=0;
        sgn = sign(gm1-gm2);
        gm = max([gm1, gm2]);
        
        df = ones(size(g1));
        for k = seq1_len/2 : numel(g1)-seq1_len/2;
            df(k) = sum(abs(g1(k).*seq1' - cp(k-seq1_len/2+1:k+seq1_len/2,grad)));
        end
        idx = g1 >1 & df <0.01;
        gm = max(g1(idx));
        idx = find(g1 ==gm);
    end
    %     idx = [find(g1 == gm), find(g2 == gm)]
    
    for k = 1:numel(idx)
        
        %         cp(idx(k)-4:idx(k)+5, grad) = sgn.*seq1.*meg_scale(grad - 1);
        cp(idx(k)-seq1_len/2+1:idx(k)+seq1_len/2, grad) = cp(idx(k)-seq1_len/2+1:idx(k)+seq1_len/2, grad).*meg_scale(grad - 1);
    end
end
t180 = (cp(idx(1)-seq1_len/2+1,1) + cp(idx(2)+seq1_len/2,1))/2/1000; %ms
