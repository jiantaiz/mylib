function [theta, t_out, G_out, r_out] =  calc_mre_phase(cp,t_rf,f_motion,phi0,disp_in_meter)
% CALC_MRE_PHASE calculates MRE phase
%    [theta, t_out, G_out, r_out] =  CALC_MRE_PHASE(cp,t_rf,f_motion,phi0,disp_in_meter) ...
%    INPUT:
%       cp - cornerpoints (or cornerpoints file name)of mre pusle sequence in Gs/cm (GE convention). 
%       t_rf - iso center of rf pulses in msec
%       f_motion - frequency of motion in Hz
%       phi0 - initial phase in radians (phase offset) of motion
%       disp_in_meter - displacement of motion in meter
%    OUTPUT:
%       theta - mre phase evolution along time on each gradient axes 
%       t - time 
%       G - Gradient waveform in mT/m (SI)
%       r - motion in meter
%    See also: calc_menc, calc_moment, calc_concomitant

% AUTHOR    : Yi Sui
% DATE      : 05/16/2017

if isstr(cp)
    cp = read_cp(cp);
end

if numel(t_rf) <2
    t90 = 0;
    t180 = t_rf;
else
    t90 = t_rf(1);
    t180 = t_rf(2);
end

if nargin <5
    A = 10e-6;%default 10um
else
    A= disp_in_meter;
end


t0 = cp (:,1)/1e3;% convert usec to msec
G0 = cp(:,2:end).*10 ; % convert Gs/cm to mT/m

dt = 8e-3; %msec

N= ceil (t0(end)/dt);

t = linspace(t0(1), t0(end), N); % in msec

G = interp1(t0,G0, t);% mT/m

gamma = 2*pi*42.576; %rad/msec/mT
% A = 10e-6; % in m 
% phi =pi/180 * phi0; %initial phase of motion

r = A*cos(2*pi*f_motion.*t/1e3 + phi0); % f_motion in Hz, r in m
r2 = A*sin(2*pi*f_motion.*t/1e3 + phi0); % f_motion in Hz, r in m

% v=1; % velocity m/msec
% r = v*t; % first moment
% 
% a=1; %acceleration m/msec^2;
% r = 0.5.*a.*t.^2; 


w = gamma.*G.*[r' r' r']; % spin angular freq in rad/usec

theta0 = zeros(sum(t<t90),3);

theta1 = cumtrapz(t(t>=t90 & t<t180), w(t>=t90 & t<t180,:)); % in rad
theta2_init = repmat(-theta1(end,:), [ sum(t>=t180), 1]);

theta2 = cumtrapz(t(t>=t180), w(t>=t180, :)) + theta2_init;
theta =[theta0;theta1; theta2];


if nargout >1
    t_out = t;
end

if nargout >2
    G_out = G;
end

if nargout >3
    r_out = r;
end
