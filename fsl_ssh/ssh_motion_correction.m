function [vol_corr, omat] = ssh_motion_correction(vol_4d,Interp,ssh2_conn)
% SSH_MOTION_CORRECTION run FSL eddy_correct to correct subject movements. 
%    vol_corr = SSH_MOTION_CORRECTION(vol_4d,Interp,ssh2_conn) 
%    INPUTS:
%       vol_4d: time series 4D volume data to correct.
%       Interp: interpolation method {'trilinear'|'spline'} 
%       ssh2_conn: ssh2 connection object (see ssh2_config)
%    OUTPUTS:
%       vol_corr: motion corrected volume data
%       omat: transformation matrices of all time points.
%    Example:
%    ssh2_conn = ssh2_login;
%    Icorr = ssh_motion_correction(I,'spline',ssh2_conn); 
%
%    Subfunctions: 
%    See also: ssh2_config, ssh2_login, fsl_wrapper, fsl_bet, bet_ssh

% AUTHOR    : Yi Sui
% DATE      : 09/12/2017
%%

if nargin<3
    ssh2_conn = ssh2_config_publickey('mr-cim', 'm165355', 'C:\Users\m165355\OneDrive\Mayo\Matlab\mr_cim_id_rsa', '');
end
if nargin <2 || isempty(Interp)    
    Interp = 'trilinear';%{trilinear,spline} 
end



tmp_dir = [tempname,'_tmp']; %get temp dir name
mkdir(tmp_dir);
tmp_infile = 'infile.nii';
tmp_outfile = 'outfile';

%save to nii file
vol_4d = single(vol_4d);
nii=make_nii(vol_4d);
save_nii(nii,fullfile(tmp_dir, tmp_infile));


%files to be sent to server
infiles = {fullfile(tmp_dir,tmp_infile)};
p = fileparts( mfilename('fullpath'));


if isreal(vol_4d)
    prg = 'eddy_correct';
else %complex number
    prg = 'eddy_correct_cplx';
end
infiles {end+1} = fullfile(p,prg);
%files to be fetched from server
outfiles = {fullfile(tmp_dir,[tmp_outfile, '.nii.gz'])};
if nargout>1
%     outfiles{end+1} =fullfile(tmp_dir,[tmp_outfile, '*.dat']);
    outfiles{end+1} =fullfile(tmp_dir,'*.dat');
end

verbose = 0;
run_in_bg=0;

cmd = sprintf('sh %s %s %s 0 %s', prg, tmp_infile,tmp_outfile,Interp);

[ssh2_conn,COMMAND_RESULT] = fsl_wrapper(cmd, infiles, outfiles, ssh2_conn,verbose, run_in_bg);

nii = load_nii(outfiles{1});
vol_corr = nii.img;

omat_files = dir(outfiles{end});
for k=1:numel(omat_files)
    fname = fullfile(tmp_dir,omat_files(k).name);    
    omat(:,:,k) = reshape(textread(fname,'%f'),[4,4])';
end
