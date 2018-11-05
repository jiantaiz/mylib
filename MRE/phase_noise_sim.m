%%!!!!!NOT WORK
%% denoising with High SNR magnitude, not possible.
N=1024;
x = linspace(-20,20,N);
[X,Y] = meshgrid(x);
x0=5;
y0=5;
R2 = (X-x0).^2 + (Y-y0).^2;
C1 = exp(-0.3*R2);

R = abs(X + 1i*Y);
r0 = abs(x0 + 1i*y0);
a=1.4;
r0 = a*r0;

C2 = exp(-(R-r0).^2);
%
% imagesc(x,x,C1+C2)
set(0,'DefaultFigureWindowStyle','docked')

figure;
contour(x,x,(C1.*C2),120);
hold on;
plot([0 x0],[0,y0],'-xr')
plot(a*[0 x0],a*[0,y0],'-xr')