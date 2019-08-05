function [ssh2_conn,COMMAND_RESULT] = fsl_wrapper(cmd, infiles, outfiles, ssh2_conn,verbose, run_in_bg)
% FSL_WRAPPER run FSL command on a remote Linux server.
%    FSL_WRAPPER(cmd, infiles, outfiles, ssh2_conn,verbose)
%         cmd      : FSL commond
%         infiles  : input files that will be sent to the remote server,
%                    it can be either a single string or a cell array of strings for
%                    multiple files.
%         outfiles : output files that will be downloaded from the remote server.
%                    it is either a single string or a cell array of strings for
%                    multiple files. '*' wildcard is supported. Note the
%                    path of outfiles is local, the remote path is always a
%                    temp folder.
%         ssh2_conn: ssh connection created by ssh2_config. If it's empty or not
%                    specified, login dialog will show.
%         verbose  : print diagnostic messages
%
%    Example:
%    fsl_wrapper('bet epi.nii epi_brain','example/epi.nii','example/epi_brain*');
%       %note that example is a local path that the output files will be downloaded to.
%    See also: bet_wrapper, flirt_wrapper

% AUTHOR    : Yi Sui
% DATE      : 06/23/2017
%%

%
% TIP:
% Set up FSL enviroment on remote server first and make sure 'bet' can be
% excuted.
% Add the fowllowing lines to the end of ~/.bashrc file on the remote server:
% export FSLDIR=/usr/local/fsl
% export PATH=$PATH:$FSLDIR/bin/
% source $FSLDIR/etc/fslconf/fsl.sh
%
% See also: ssh2_config, ssh2_command

% AUTHOR    : Yi Sui
% DATE      : 06/23/2017
%%

if ~iscell(infiles)
    infiles= {infiles};
end

if ~iscell(outfiles)
    outfiles = {outfiles};
end
flag_close_ssh = false;
if  ~exist('ssh2_conn','var')|| isempty(ssh2_conn)
    if nargout <1
        flag_close_ssh = true;
    end
    private_key = 'C:\Users\m165355\OneDrive\Mayo\Matlab\mr_cim_id_rsa';
    if exist(private_key,'file') == 2
        ssh2_conn = ssh2_config_publickey('mr-cim', 'm165355',private_key, '');
    else
        ssh2_conn = ssh2_login('m165355',[],'mr-cim');
    end
        
    if isempty(ssh2_conn) % cancelled
        COMMAND_RESULT=[];
        return
    end
end

if ~exist('verbose','var') || isempty(verbose)
    verbose = 0;
end

if  ~exist('run_in_bg','var'), run_in_bg=0; end

if run_in_bg
    if isstruct(ssh2_conn)
        ssh2_conn = ssh2_close(ssh2_conn);%error occurs in workers if connection is still open.
        ssh2_conn.connection=[];
        ssh2_conn.command_session=[];
    end
    COMMAND_RESULT=parfeval(@fsl_wrapper,1,cmd, infiles, outfiles, ssh2_conn,verbose);
    return
end

%
if verbose == 0 ;
    verbose = @(varargin)fprintf('');
else
    verbose = @(varargin)fprintf(varargin{:});
end

if ~isempty (regexpi(cmd,'rm .*(-.*r|*)')); % safety check
    warning('The command %s is not allowed, it will delete massive files.',cmd);
    COMMAND_RESULT=[];
    return;
end

remote_tmp_dir = ['tmp_fslwrapper_',num2str(matlab.internal.timing.timing('cpucount'))];
ssh2_setup;
try
    %send file to server
    ssh2_conn = ssh2_command(ssh2_conn,sprintf('mkdir %s',remote_tmp_dir));
    
    for k=1:numel(infiles)
        if ~isempty(infiles{k})
            verbose('uploading %s\n', infiles{k});
            [infile_path,infile_name,infile_ext] = fileparts(infiles{k});
            ssh2_conn = scp_put(ssh2_conn,[infile_name,infile_ext],remote_tmp_dir,infile_path);
        end
    end
    
    %execute cmd on server
    verbose('running on %s ...\n',ssh2_conn.hostname);
    cmd = sprintf('cd %s && %s',remote_tmp_dir,cmd);
    verbose('%s\n',cmd);
    tic
    
    [ssh2_conn, COMMAND_RESULT]= ssh2_command(ssh2_conn,cmd);
    if nargout<2
       fprintf('%s\n',COMMAND_RESULT{:});
    end
    verbose('%s\n',COMMAND_RESULT{:});
    verbose('Done! (%g sec)\n',toc);
    
    %retrieve files from server
    for k = 1:numel(outfiles)
        if ~isempty(outfiles{k})
            verbose('downloading %s\n', outfiles{k});
            [outfile_path,outfile_name,outfile_ext] = fileparts(outfiles{k});
            [ssh2_conn, files ]= ssh2_command(ssh2_conn,sprintf('ls %s/%s', remote_tmp_dir, [outfile_name,outfile_ext]));
            ssh2_conn = scp_get(ssh2_conn,files,outfile_path);
        end
    end
    %delete the tmp files on server
    if isempty( strfind(remote_tmp_dir,'*')) && ~isempty(strfind(remote_tmp_dir,'tmp_fslwrapper_')) % to be safe, check no '*' in file name
        verbose('delete tmp files in %s\n',remote_tmp_dir);
        ssh2_conn = ssh2_command(ssh2_conn,sprintf('rm -rf %s',remote_tmp_dir));
    end
    
    if flag_close_ssh        
        ssh2_conn = ssh2_close(ssh2_conn);
        verbose('SSH connection closed.\n');
    end
catch err
%     ssh2_conn = ssh2_close(ssh2_conn);
    rethrow(err);
end
