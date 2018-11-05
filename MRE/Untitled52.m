fname = 'C:\Users\m165355\Documents\Matlab\Data\BrainMRE\cpfiles_zy_brainmre\cornerPointsRot2.0.all';

fname = 'epimre.tgt.pdx.3';
TE= 66;
t_rf = [5.42,35.765];
%%
fname = 'C:\Users\m165355\Documents\Matlab\Data\epic\epimre.tgt.pd';
TE= 82;%cpfiles_zy_brainmre\cornerPointsRot.0.all
t_rf = [7.864,44.675];%cpfiles_zy_brainmre\cornerPointsRot.0.all
%%
 [t_out,wave_out] = plotPdFile(fname);
%  figure; plot(t_out,wave_out);
cp =[t_out,wave_out*4];
%%
 r=[1,1,-1.5]/10;
 [t, G, m] =  calc_concomitant(cp, r , t_rf);
 figure;plot(t,G,t,m)
%%
x= [-20:20]/100; 
y= [-20:20]/100; 
z= 1/100;
r={x,y,z};
m =  calc_concomitant_map(cp,r,t_rf,TE,3000);

figure;surf(x,y,squeeze(m))

