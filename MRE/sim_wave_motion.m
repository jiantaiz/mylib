
%simulate the motion of combined sheer wave and compressed wave.
As = 9;
Ac =4.5;
w =2*pi*60;%
t=linspace(0,1/60*2,100);

wave_num = 4;
r = linspace(0,0.4,100);%meters

phi_s = 2*pi*r.*wave_num;

% phi_s = 0;
phi_c = 0;
figure;
t0=linspace(0,1/12,100);
D=[]
plt = plot([1 2 3]);
for t=t0;   


ds = As* cos(w*t + phi_s);
dc = Ac* cos(w/4*t + phi_c);

d=ds+dc;
D=[D,d'];
 plot(r,d); 
% set(plt,'xdata',r,'ydata',d);
ylim([-(As+Ac),(As+Ac)]); xlim([0 r(end)]);
 hold on
plot(r(6),d(6),'o',r(30),d(30),'o')
plot(r(15),d(15),'x',r(45),d(45),'x')
% plot(t,d)
xlabel('r')
grid on

pause(0.2)
end
grid on
%%
imdisp(D);
xlabel('t');
ylabel('r')
figure;
Dfft= ifft(D,[],2);
plot(2*abs(Dfft(:,6)));hold on;plot(abs(max(D,[],2) ),'--')
