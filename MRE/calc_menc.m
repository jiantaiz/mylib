function [menc,ph_offset,menc_cplx] = calc_menc(cp,t_rf,echo_pos,f_motion)
% CALC_MENC     Calculate MRE motion encoding sensitivity (menc) using actural waveform
%
%   [menc,phase] =  CALC_MENC(cp,t_rf,TE,f_motion)
%   Calculate MRE motion encoding sensitivity (menc) using actural waveform
%INPUT:
% cp - cornerpoints (or cornerpoints file) of gradient waveform
% t_rf - iso center time of rf pulses (90 and/or 180)
% echo_pos - echo position (note: not actural echo time, but where the echo is on the time axes)
%OUTPUT:
% menc : motion encoding sensitivity in rad/um
% phase: the offset of motion that gives the maxium accumulated mre phase.

N=4;

ph=[];
for phi0 = [0:2*pi/N:2*pi-2*pi/N];
    [theta,t] = calc_mre_phase(cp,t_rf,f_motion,phi0,1e-6);
    TE_idx = find(t>=echo_pos);
    TE_idx = TE_idx(1);
    ph  = [ph ;theta(TE_idx,:)];
end

ph_fft = fft(ph)/N;

menc = 2*abs(ph_fft(2,:));%rad/um
if nargout >1
    ph_offset = angle(ph_fft(2,:));
end

if nargout >2
    menc_cplx = menc.*exp(1i*ph_offset);
end
