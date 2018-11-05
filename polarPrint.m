function polarPrint(A,rad_or_deg,num_fmt)
% POLARPRINT ...
%    POLARPRINT(A,rad_or_deg,num_fmt) ...
%
%    Example:
%    ... 
%
%    Subfunctions: 
%    See also: 

% AUTHOR    : Yi Sui
% DATE      : 05/16/2017
%%%polarPrint(A,rad_or_deg,num_fmt)if nargin < 2;    rad_or_deg = 'deg';endif nargin < 3    g = 'g';else    g = num_fmt;endfmt = ['%' g  '%c' '%' g '%c' '\t'];mag = abs(A);ang = rad2deg(angle(A));[nx, ny] = size(A);angsign = ones(1,ny).*8736;% /_switch rad_or_deg    case 'rad'        degsign = ones(1,ny).*' ';        ang = angle(A);    case 'deg'        degsign = ones(1,ny).*186; %º        ang = rad2deg(angle(A));endfor k = 1:nx    fprintf(fmt, [mag(k,:);angsign;ang(k,:);degsign] );    fprintf('\n')end
