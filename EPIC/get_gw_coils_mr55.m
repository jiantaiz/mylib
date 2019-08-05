function gwcoefs = get_gw_coils_mr55()
% #{sdc@mr55}[101] cat /w/config/gw_coils.dat
% #Gradwarp Coefficient file for a XRMW Body Gradient Coil - Spherical Harmonic Correction
% GRADWARPTYPE            1
SCALEX(1) =                 0.0;
SCALEX(2) =                -3.2192e-3;
SCALEX(3) =                 -6.8787e-4;
SCALEX(4) =                 1.3123e-5;
SCALEX(5) =                 -2.9957e-7;
SCALEX(6) =                 -4.2178e-8;
SCALEX(7) =                 2.5713e-10;
SCALEX(8) =                 1.1517e-11;
SCALEX(9) =                 -3.0381e-13;
SCALEX(10) =                -8.4245e-15;
SCALEY(1) =                 0.0;
SCALEY(2) =                 -2.3675e-3;
SCALEY(3) =                 -6.2607e-4;
SCALEY(4) =                 1.3177e-5;
SCALEY(5) =                 -6.4110e-7;
SCALEY(6) =                 -4.5042e-8;
SCALEY(7) =                 4.1225e-10;
SCALEY(8) =                 1.0885e-11;
SCALEY(9) =                 -5.0161e-13;
SCALEY(10) =                -4.8482e-15;
SCALEZ(1) =                 0.0;
SCALEZ(2) =                 9.4428e-6;
SCALEZ(3) =                 -7.7195e-5;
SCALEZ(4) =                 -6.7632e-9;
SCALEZ(5) =                 -6.5554e-7;
SCALEZ(6) =                 7.1369e-12;
SCALEZ(7) =                 7.3912e-10;
SCALEZ(8) =                 -9.6227e-16;
SCALEZ(9) =                 -9.1711e-14;
SCALEZ(10) =                -2.5266e-19;
% DELTA                   0.0
% #created 10 June 2009, includes 7th order correction to Z5 term
% # SCALEX(1 - SCALEX(5, SCALEY(1- SCALEY(5, SCALEZ(1 - SCALEZ(5 (The spherical
% # harmonic coefficients) = need to be computed and recorded in this file.
% # However, the numbers entred in the file should have the multiplicative
% # factors for each of the term already incorporated.
% # Therefore,
% # SCALEX(1 *= 1, SCALEX(2 *= 3,   SCALEX(3 *= 3/2, SCALEX(4 *= 5/2, SCALEX(5 *= 15/8
% # SCALEY(1 *= 1, SCALEY(2 *= 3,   SCALEY(3 *= 3/2, SCALEY(4 *= 5/2, SCALEY(5 *= 15/8
% # SCALEZ(1 *= 1, SCALEZ(2 *= 1/2, SCALEZ(3 *= 1/2, SCALEZ(4 *= 1/8, SCALEZ(5 *= 1/8.
% #
% #
% #       X:
% #               hc31
% #               hc51
% #       Y:
% #               hc31
% #               hc51
% #       Z:
% #               hc30
% #               hc50

gwcoefs = [SCALEX',SCALEY',SCALEZ'];