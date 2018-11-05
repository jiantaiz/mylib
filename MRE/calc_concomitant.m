function [t, G, m] =  calc_concomitant(cp,r,t_rf,B0)
%[t, G, m] =  calc_concomitant(cp, r , t_rf)
% calculate comcomitant field of givin waveform (cp)
% cp: time and gradients
%     cp(:,1) - time in usec
%     cp(:,2:4) - X Y Z gradients in mT/m
% r: 1x3 location vector for one point [x,y,z], 
% t_rf: 1x2 vector time position of 90 and 180 RF
% B0: main magnetic field in mT. e.g. 3000 for 3T magnet
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


dt = 0.02; %msec

N= ceil (t0(end)/dt)

t = linspace(t0(1), t0(end), N); % in msec

G = interp1(t0,G0, t);% mT/m

gamma = 2*pi*42.576; %rad/msec/mT

Gx = G(:,1);
Gy = G(:,2);
Gz = G(:,3);

x = r(1);
y = r(2);
z = r(3);

if nargin < 4
    B0 = 3000 %mT, 3T = 3000mT
end

Bc = 1/2/B0 .* (Gx.^2 * z.^2 + Gy.^2 .* z.^2 + Gz.^2 .* (x.^2 + y.^2)/4 - Gx.*Gz.*x.*z - Gy.*Gz.*y.*z );
w = gamma.*Bc;
theta0 = zeros(sum(t<t90),1);
theta1 = cumtrapz(t(t>=t90 & t<t180), w(t>=t90 & t<t180)); % in rad
theta2_init = repmat(-theta1(end), [ sum(t>=t180), 1]);
theta2 = cumtrapz(t(t>=t180), w(t>=t180)) + theta2_init;
m =[theta0;theta1; theta2];
