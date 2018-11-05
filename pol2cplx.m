function cplx = pol2cplx(mag, deg)
% POL2CPLX returns a complex number from its polar form
%    cplx = POL2CPLX(mag, deg) returns a complex number given magnitude
%    and angle
%

% AUTHOR    : Yi Sui
% DATE      : 05/16/2017
%%

rad =  deg2rad(deg);
cplx = mag.*exp(1i.*rad);
