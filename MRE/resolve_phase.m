function ph_small =  resolve_phase(ph_small,phs_big, k)
% RESOLVE_PHASE experimental...
%    ph_small =  RESOLVE_PHASE(ph_small,phs_big, k) ...
%
%    Example:
%    ... 
%
%    Subfunctions: 
%    See also: 

% AUTHOR    : Yi Sui
% DATE      : 05/16/2017

% if ph_big =7*phs_small, we can guess if the small phase has 2pi jump or
% not using the code below, k needs to be half integer e.g 3.5

idx = find(ph_small >0);
alpha = pi - ph_small(idx); %
beta =angle(exp(1i* (k*pi -k*alpha ))); % to see if it is really positive
ph = phs_big;
delta = abs(ph(idx) - beta);
idx2neg = idx(delta > pi/2.5 );


idx = find(ph_small <0);
alpha = ph_small(idx)-pi; %
beta =angle(exp(1i* (-k*pi +k*alpha )));% to see if it is really negtive
ph = phs_big;
delta = abs(ph(idx) - beta);
idxpos = idx(delta > pi/2.5 );

ph_small(idx2neg) = ph_small(idx2neg)-2*pi;
ph_small(idxpos) = ph_small(idxpos)+2*pi;
