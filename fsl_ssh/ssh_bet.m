function [brain,brain_mask,skull_mask,scalp_mask] = ssh_bet(mat,voxsize,ssh2_conn,g,f,dorobust,betopt)
% SSH_BET run FSL bet (brain extraction tool) on remote server via ssh
%    [brain,brain_mask] = SSH_BET(mat,voxsize,ssh2_conn,g,f,dorobust,betopt)
%    [brain,brain_mask,skull_mask,scalp_mask] = SSH_BET(...); %returns
%                                               additional skull and sclp masks.
%    Examples:
%    [brain,brain_mask] = SSH_BET(mat,voxsize); % mat is a 3D brain volume.
%
%    % Config a ssh2 connection
%    ssh2_conn = ssh2_config(HOSTNAME,USERNAME,PASSWORD, PORT);
%    % or use publickey
%    ssh2_conn = ssh2_config_publickey('mr-cim', 'm165355', 'C:\Users\m165355\OneDrive\Mayo\Matlab\mr_cim_id_rsa3', '');

%    % Or simply use the login window.
%    ssh2_conn = ssh2_login;
%    [brain,brain_mask] = ssh_bet(mat,voxsize,ssh2_conn,0,0.4,1,'-S'); % -S: eye & optic nerve cleanup
%
%    

%
% TIPS:
% Set up FSL enviroment on remote server first and make sure 'bet' can be
% excuted.
% Add the fowllowing lines to the end of ~/.bashrc file on the remote server:
% export FSLDIR=/usr/local/fsl
% export PATH=$PATH:$FSLDIR/bin/
% source $FSLDIR/etc/fslconf/fsl.sh
%
% Main bet2 options (just for your reference, not all options work with ssh_bet):
%   -o          generate brain surface outline overlaid onto original image
%   -m          generate binary brain mask
%   -s          generate approximate skull image
%   -n          don't generate segmented brain image output
%   -f <f>      fractional intensity threshold (0->1); default=0.5; smaller values give larger brain outline estimates
%   -g <g>      vertical gradient in fractional intensity threshold (-1->1); default=0; positive values give larger brain outline at bottom, smaller at top
%   -r <r>      head radius (mm not voxels); initial surface sphere is set to half of this
%   -c <x y z>  centre-of-gravity (voxels not mm) of initial mesh surface.
%   -t          apply thresholding to segmented brain image and mask
%   -e          generates brain surface as mesh in .vtk format
%
% Variations on default bet2 functionality (mutually exclusive options):
%   (default)   just run bet20
%   -R          robust brain centre estimation (iterates BET several times)
%   -S          eye & optic nerve cleanup (can be useful in SIENA)
%   -B          bias field & neck cleanup (can be useful in SIENA)
%   -Z          improve BET if FOV is very small in Z (by temporarily padding end slices)
%   -F          apply to 4D FMRI data (uses -f 0.3 and dilates brain mask slightly)
%   -A          run bet2 and then betsurf to get additional skull and scalp surfaces (includes registrations)
%   -A2 <T2>    as with -A, when also feeding in non-brain-extracted T2 (includes registrations)
%
% Miscellaneous options:
%   -v          verbose (switch on diagnostic messages)
%   -h          display this help, then exits
%   -d          debug (don't delete temporary intermediate images)
% See also: bet_wrapper, fsl_wrapper, ssh2_config

% AUTHOR    : Yi Sui
% DATE      : 06/23/2017
%%

if ~exist('g','var')|| isempty(g), g = 0; end
if ~exist('f','var')|| isempty(f), f = 0.5; end
if ~exist('dorobust','var') || isempty(dorobust) , dorobust = 0; end
if ~exist('betopt','var')|| isempty(betopt), betopt = ''; end
if nargin<1, mfile_showhelp; return; end
if ~exist('ssh2_conn','var'), ssh2_conn=[]; end
if dorobust
    betopt =[betopt,' -R']; %robust
end
if nargout > 1
    betopt =[betopt,' -m']; % get brain mask
end
if nargout > 2
    betopt =[betopt,' -A']; % get additional skull and scalp surfaces.
end

tmp_dir = [tempname,'_tmp']; %get temp dir name
mkdir(tmp_dir);

tmp_infile = 'infile.nii';
tmp_outfile = 'outfile';

%save to nii file
nii=make_nii(mat,voxsize);
save_nii(nii,fullfile(tmp_dir, tmp_infile));

betopt = sprintf('-f %2.2f -g %2.2f %s', f, g, betopt); %bet options
%%

try
    % run bet
    bet_wrapper(fullfile(tmp_dir,tmp_infile), fullfile(tmp_dir,tmp_outfile),betopt,ssh2_conn);
    
    
    %read in .nii files.
    nii = load_nii( fullfile(tmp_dir, [tmp_outfile,'.nii.gz']));
    brain = nii.img;
    if nargout>1
        nii = load_nii( fullfile(tmp_dir, [tmp_outfile,'_mask.nii.gz']));
        brain_mask = nii.img;
    end
    if nargout>2
        nii = load_nii(fullfile(tmp_dir, [tmp_outfile,'_skull_mask.nii.gz']));
        skull_mask = nii.img;
        nii = load_nii(fullfile(tmp_dir, [tmp_outfile,'_outskin_mesh.nii.gz']));
        scalp_mask = nii.img;
    end
catch err
    ssh2_conn = ssh2_close(ssh2_conn);
    %delete the tmp files on local PC
    if isempty( strfind(tmp_dir,'*'))&& ~isempty( strfind(tmp_dir,'_tmp'))
        rmdir(tmp_dir,'s');
    end
    rethrow(err);
end
%delete the tmp files on local PC
if isempty( strfind(tmp_dir,'*'))&& ~isempty( strfind(tmp_dir,'_tmp'))
    rmdir(tmp_dir,'s');
end
