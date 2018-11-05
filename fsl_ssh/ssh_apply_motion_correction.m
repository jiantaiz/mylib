function vol_corr = ssh_apply_motion_correction(vol_4d,omat,Interp,ssh2_conn)
% SSH_MOTION_CORRECTION run FSL eddy_correct to correct subject movements. 
%    vol_corr = SSH_MOTION_CORRECTION(vol_4d,Interp,ssh2_conn) 
%    INPUTS:
%       vol_4d: time series 4D volume data to correct.
%       omat: transformation matrices of all time points.
%       Interp: interpolation method {'trilinear'|'spline'} 
%       ssh2_conn: ssh2 connection object (see ssh2_config)
%    OUTPUTS:
%       vol_corr: motion corrected volume data
%       
%    Example:
%     % 4D MRE data subject movement correction
%     ssh2_conn = ssh2_login;
%     vol = cimgs(:,:,:,:);
%     [Vcorr, mat]  = ssh_motion_correction(abs(vol),'spline'); %get trans matrices using magnitude images     
%     [Vcorr_r ]  = ssh_apply_motion_correction(real(vol),mat,'spline'); % apply to real and imaginary parts seperately
%     [Vcorr_i ]  = ssh_apply_motion_correction(imag(vol),mat,'spline');
%     cimgs_corr = Vcorr_r + 1i*Vcorr_i;
%     cimgs_corr = reshape(cimgs_corr, size(cimgs));
%
%    Subfunctions: 
%    See also: ssh_motion_correction, ssh2_config, ssh2_login, fsl_wrapper, fsl_bet, bet_ssh

% AUTHOR    : Yi Sui
% DATE      : 09/12/2017
%%

if nargin<4
    ssh2_conn = ssh2_config_publickey('mr-cim', 'm165355', 'C:\Users\m165355\OneDrive\Mayo\Matlab\mr_cim_id_rsa', '');
end
if nargin <3 || isempty(Interp)    
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
infiles {end+1} = fullfile(p,'apply_eddy_correct');

for k=1:size(omat,3)
    omat_file  = fullfile(tmp_dir, sprintf('%s_tmp%04d_mat.dat',tmp_outfile,k-1));
    fid = fopen(omat_file,'w');
    fprintf(fid,'%g %g %g %g\n',omat(:,:,k)');
    fclose(fid);
    
    infiles{end+1} = omat_file;
end

%files to be fetched from server
outfiles = {fullfile(tmp_dir,[tmp_outfile, '.nii.gz'])};

verbose = 0;
run_in_bg=0;

cmd = sprintf('sh apply_eddy_correct %s %s 0 %s', tmp_infile,tmp_outfile,Interp);
[ssh2_conn,COMMAND_RESULT] = fsl_wrapper(cmd, infiles, outfiles, ssh2_conn, verbose, run_in_bg);

nii = load_nii(outfiles{1});
vol_corr = nii.img;

