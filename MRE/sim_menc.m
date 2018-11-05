
N = 8;
f_motion =60; % Hz
% t_rf = [7.661 29.698]; % msec for fcomp3
% t_rf = [5.592 38.6]; % msec  'cp_70_fcomp2.all'

% t_rf =[89.730 117.4]; % msec for cp_cmre_140_fc2.all
% 
% t_rf = [92.830 126.500 ]; % 'cp_cmre_140_fc3.all';
% 
% 
% t_rf = [192.8 225.80];%cp_cmre_75_140_fc3.all
% t_rf = [89.75 117.9];%cp_cmre_140_asymeg.all
% t_rf =[90.850 118.0];%cp_cmre_140_asymeg_noear.all

% TE= 51.360;% for fcomp3
% TE=72.45;%fcomp2

% TE = 145.050; %msec cp_cmre_140_fc2.all
% 

% TE = 161.500; %'cp_cmre_140_fc3.all';
% 
% TE = 258.5;%cp_cmre_75_140_fc3.all
% 
% TE = 146.7;%cp_cmre_140_asymeg
% TE = 145.3; %cp_cmre_140_asymeg_noear.all

% fname = 'cp_70_fcomp3.all'
% fname = 'cp_70_fcomp2.all';
% fname = 'cp_cmre_140_fc2.all';
% fname = 'cp_cmre_140_fc3.all';
% fname = 'cp_cmre_75_140_fc3.all';
% fname = 'cp_cmre_140_asymeg.all';
% fname = 'cp_cmre_140_asymeg_noear.all'
% fname = 'cp_cmre_140_asym_2drf.all';
% fname = 'cp_cmre_140_conv.all';


fname = 'C:\Users\m165355\Documents\Matlab\Data\BrainMRE\cpfiles_zy_brainmre\cornerPointsRot.0.all';
%  fname = 'C:\Users\m165355\Documents\Matlab\Data\BrainMRE\cpfiles_zy_brainmre\cornerPoints.all';
% fname = 'C:\Users\m165355\Documents\Matlab\Data\BrainMRE\cpfiles_zy_brainmre\cornerPointsRot2.0.all'; % fc2 meg55 se6
fname = 'C:\Users\m165355\Documents\Matlab\Data\BrainMRE\cpfiles_zy_brainmre\cornerPointsRot4.0.all';

TE= 67.97;%cpfiles_zy_brainmre\cornerPointsRot.0.all
t_rf = [5.360,36.46];%cpfiles_zy_brainmre\cornerPointsRot.0.all



% fname = 'C:\Users\m165355\Documents\Matlab\Data\BrainMRE\cpfiles_zy_brainmre\cornerPoints.fc3.meg68.se10.all';
% TE = 48.7;
% t_rf = [5.360 26.865];%cornerPoints.fc3.meg68.se10.all
% 

% fname = 'C:\Users\m165355\Documents\Matlab\Data\BrainMRE\cpfiles_zy_brainmre\cornerPoints.fc3.meg60.se9.all';
% TE = 52.1;
% t_rf = [5.360 28.56];%cornerPoints.fc3.meg68.se10.all

cp = read_cp(fname);
% 
% time_scale=1;
% cp(:,1) = cp(:,1).*time_scale;
% t_rf = t_rf.*time_scale;
% TE = TE.*time_scale;


amp = [];
ang =[];
f_motions = 60;
for f_motion = f_motions
    ph=[];
    for phi0 = [0:2*pi/N:2*pi-2*pi/N];
      [theta,t] = calc_mre_phase(cp,t_rf,f_motion,phi0);
      TE_idx = find(t>=TE);
      TE_idx = TE_idx(1);
      ph  = [ph ;theta(TE_idx,:)];
      
      if 0
          plot(t,rad2deg(theta(:,1)), t, r*1e6);
          hold on
%           pause
      end
    end
    ph_fft = fft(ph)/N;
    
    amp = [amp; 2*abs(ph_fft(2,:))];
    ang = [ang;angle(ph_fft(2,:))];
    
    
%     ang_init= -rad2deg(mean(ang(end,1)));
    ang_init= -mean(ang(1,3));
    [theta,t,G, r] = calc_mre_phase(cp, t_rf,f_motion,ang_init);
    
end
%%
figure;
 [t, G, m] =  calc_moment(cp,2, t_rf);
for k=1:3
    ax(k) = subplot(3,1,k);
    plot(t,G(:,k),t,m(:,k)./max(m(:,k)) *30 )
     hold on
    plot(t,rad2deg(theta(:,k)), t, r*1e6);
    grid on
end
linkaxes(ax);

%%
figure
plot(t,[G(:,:)]);
grid on;
% plot(f_motions,amp)
% plot(ph)
% 
%%
cp = 'C:\Users\m165355\Documents\Matlab\Data\BrainMRE\cpfiles_zy_brainmre\cornerPointsRot.0.all';

 [menc_pos,phase_pos] =  calc_menc(cp,t_rf,TE,f_motion);
 
 MENC_pos = menc_pos.*exp(1i.*phase_pos)
 
cp = 'C:\Users\m165355\Documents\Matlab\Data\BrainMRE\cpfiles_zy_brainmre\cornerPointsRot2.0.all'; % fc2 meg55 se6
[menc_neg,phase_neg] =  calc_menc(cp,t_rf,TE,f_motion);

 MENC_neg = menc_neg.*exp(1i.*phase_neg)
 
cp = 'C:\Users\m165355\Documents\Matlab\Data\BrainMRE\cpfiles_zy_brainmre\cornerPointsRot3.0.all'; % fc2 meg55 se6
[menc_neg2,phase_neg2] =  calc_menc(cp,t_rf,TE,f_motion);

 MENC_neg2 = menc_neg2.*exp(1i.*phase_neg2)
 
 cp = 'C:\Users\m165355\Documents\Matlab\Data\BrainMRE\cpfiles_zy_brainmre\cornerPointsRot4.0.all'; % fc2 meg55 se6
[menc_img,phase_img] =  calc_menc(cp,t_rf,TE,f_motion);

 MENC_img = menc_img.*exp(1i.*phase_img)
 