function w = fermiwin(n,radius,width)
% FERMIWIN Fermi Filter
%    w = FERMIWIN(n,radius,width)
%    n: length of window function
%    radius: fermi radius in pixel, the frequencies less than which are in the pass-band; the freq greater than which are attenuated smoothly to zero
%    width: fermi width in pixel, larger number gives stronger low-pass filtering
%    Example:
%    w = FERMIWIN(128,32,2.5) 
%
%    Subfunctions: 
%    See also: 

% AUTHOR    : Yi Sui
% DATE      : 03/20/2019
%%
%
if numel(n)==1 % 1D window
    x = abs([0:n-1] - n/2);
    w = 1./(1+exp((x-radius)/width));
elseif numel(n) == 2 % 2D window
    y = [0:n(1)-1] - n(1)/2;
    x = [0:n(2)-1] - n(2)/2;
    [X,Y] = meshgrid(x,y);
    R =sqrt(X.^2+Y.^2);
    w = 1./(1+exp((R-radius)./width));
end
w = w./max(w(:)); %normalize
    