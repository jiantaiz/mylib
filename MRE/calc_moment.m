function [t, G, m] =  calc_moment(cp, order, t_rf)
% CALC_MOMENT calculates gradient moments
%    [t, G, m] =  CALC_MOMENT(cp, order, t_rf) returns time t, gradient G,
%    and its gradient moment m. cp is cornerpoint of gradient or
%    cornerpoint file, order specifies the order of moments, t_rf gives the
%    time location (in msec) of RF pulses [t90 t180].
%
%    Example:
%    [t, G, m] =  calc_moment(cp, 1, [1.4 50.1])
%
%    See also: calc_mre_phase, calc_menc, calc_concomitant

% AUTHOR    : Yi Sui
% DATE      : 05/16/2017
%%

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


t0 = cp (:,1)/1e3;% convert usec to msec
G0 = cp(:,2:end).*10 ; % convert Gs/cm to mT/m

dt = 1e-4; %msec

N= ceil (t0(end)/dt);

t = linspace(t0(1), t0(end), N+1)'; % in msec

G = interp1(t0,G0, t);% mT/m

gamma = 2*pi*42.576/100; %rad/msec/mT

if order>0
    r = (t/1000).^order; %in m
else
    r =(1e-3).*t.^0;%in m
end
% w = gamma.*G.*[r' r' r']; % spin angular freq in rad/usec
w = gamma.* bsxfun(@times, G, r);% spin angular freq in rad/usec


theta0 = zeros(sum(t<t90),size(G,2));

theta1 = cumtrapz(t(t>=t90 & t<t180), w(t>=t90 & t<t180,:)); % in rad
theta2_init = repmat(-theta1(end,:), [ sum(t>=t180), 1]);

theta2 = cumtrapz(t(t>=t180), w(t>=t180, :)) + theta2_init;
m =[theta0;theta1; theta2];
