function [ssh2_conn, COMMAND_RESULT] = bet_wrapper(infile, outfile, opt, ssh2_conn,verbose,run_in_bg)
% BET_WRAPPER runs FSL bet (brain extraction tool) on a remote server
%    
%    BET_WRAPPER(infile, outfile, opt, ssh2_conn,verbose,run_in_bg)
%        infile   : input file.
%        outfile  : base name of output file(s).
%        opt      : bet2 options see below.
%        ssh2_conn: ssh connection created by ssh2_config. If empty or not
%                   specified, login dialog will show.
%        verbose  : diagnostic messages
%        run_in_bg: run the command in background via parfeval, so your matlab command
%                   window will not be blocked.
%    Examples:
%        BET_WRAPPER example/epi.nii  example/epi_brain;
%        BET_WRAPPER('example/epi.nii', 'example/epi_brain');
%
%        BET_WRAPPER example/epi.nii example/epi_brain '-f 0.4 -g -0.1 -m';
%        BET_WRAPPER('example/epi.nii', 'example/epi_brain', '-f 0.4 -g -0.1 -m');
%
%        ssh2_conn = ssh2_login;
%        [ssh2_conn, COMMAND_RESULT] = BET_WRAPPER('example/epi.nii', 'example/epi_brain', '-f 0.4 -g -0.1 -m -A',ssh2_conn,1);
%
% Usage:    bet <input> <output> [options]
% Main bet2 options (just for your reference, not all options work with bet_ssh):
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
%   (default)   just run bet2
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
%
% TIP:
% Set up FSL enviroment on remote server first and make sure 'bet' can be
% excuted.
% Add the fowllowing lines to the end of ~/.bashrc file on the remote server:
% export FSLDIR=/usr/local/fsl
% export PATH=$PATH:$FSLDIR/bin/
% source $FSLDIR/etc/fslconf/fsl.sh
%
% See also: bet_ssh, fsl_wrapper, ssh2_config, ssh2_command

% AUTHOR    : Yi Sui
% DATE      : 06/23/2017
%%
if ~exist('opt','var'), opt = ''; end
if  ~exist('ssh2_conn','var'), ssh2_conn=[]; end
if ~exist('verbose','var') || isempty(verbose)
    verbose = 0;
end
if  ~exist('run_in_bg','var'), run_in_bg=0; end
% if run_in_bg
%     if isempty(ssh2_conn)
%         ssh2_conn = ssh2_login;
%     end
%     COMMAND_RESULT=parfeval(@bet_wrapper,1,infile, outfile, opt, ssh2_conn,verbose,0);
%     return
% end


[infile_path,infile_name,infile_ext] = fileparts(infile);
[outfile_path,outfile_name,outfile_ext] = fileparts(outfile);
try    
    cmd = sprintf('bet %s %s %s',[infile_name,infile_ext],outfile_name, opt);
    [ssh2_conn COMMAND_RESULT]= fsl_wrapper(cmd, infile,[outfile,'*'], ssh2_conn,verbose, run_in_bg);
catch err
    ssh2_conn = ssh2_close(ssh2_conn);
    rethrow(err);
end
