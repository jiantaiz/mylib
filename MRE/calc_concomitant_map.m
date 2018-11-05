function m =  calc_concomitant_map(cp,r,t_rf,TE,B0)
% m =  calc_concomitantmap(cp,r,t_rf,TE,B0)
% calculate comcomitant field of givin waveform (cp)
% cp: time and gradients
%     cp(:,1) - time in usec
%     cp(:,2:4) - X Y Z gradients in mT/m
% r: a cell of {x, y, z} deines the grid on the intersection of x,y,z.
% e.g r = {[1 2 3], [1 2 3], 4} defines a 3x3x1=9 points grid.
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


dt = 0.01; %msec

N= ceil (t0(end)/dt);

t = linspace(t0(1), t0(end), N); % in msec

G = interp1(t0,G0, t);% mT/m

gamma = 2*pi*42.576; %rad/msec/mT

Gx = G(:,1); 
Gy = G(:,2); 
Gz = G(:,3); 
if nargin < 5
    
    B0 = 3000 %mT, 3T = 3000mT
end
if iscell(r)
    xv = r{1}; nx = numel(xv);
    yv = r{2}; ny = numel(yv);
    zv = r{3}; nz = numel(zv);
    m = zeros(nx,ny,nz);
else
    error ('r should be a cell')
end
for kx = 1:nx;
        x = xv(kx);
    for ky = 1:ny;
        y = yv(ky);
        for kz = 1:nz
            z = zv(kz);
            Bc = 1/2/B0 .* (Gx.^2 * z.^2 + Gy.^2 .* z.^2 + Gz.^2 .* (x.^2 + y.^2)/4 - Gx.*Gz.*x.*z - Gy.*Gz.*y.*z );
            w = gamma.*Bc;
            theta1 = trapz(t(t>=t90 & t<t180), w(t>=t90 & t<t180)); % in rad
            m(kx,ky,kz) = -theta1 + trapz(t(t>=t180 & t<=TE), w(t>=t180 & t<=TE)) ;

        end
    end
end
