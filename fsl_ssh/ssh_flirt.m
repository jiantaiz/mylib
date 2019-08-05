function [regvol,omat] = ssh_flirt(inputvol,inputvoxsize,refvol,refvoxsize,ssh2_conn,options,verbose)
% ssh_flirt linear registration
%    [regvol,omat] = ssh_flirt(inputvol,inputvoxsize,refvol,refvoxsize)
%    linear registrate inputvol to refvol, the voxsize is in mm.
%    [regvol,omat] = ssh_flirt(...,ssh2_conn,options,verbose) specifies
%    optional ssh2_conn, additional FLIRT options, and verbose flag.
%
%    Examples:
%    [regvol omat] = ssh_flirt(inputvol,[1.875 1.875 3], refvol,[2 2 2]);
%
%    % Config a ssh2 connection
%    ssh2_conn = ssh2_config(HOSTNAME,USERNAME,PASSWORD, PORT);
%    % Or simply use the login window.
%    ssh2_conn = ssh2_login;
%    [regvol omat] = ssh_flirt(inputvol,[1.875 1.875 3], refvol,[2 2 2],ssh2_conn);
%
%    [regvol omat] = ssh_flirt(inputvol,[1.875 1.875 3], refvol,[2 2 2],ssh2_conn,'-dof 6',1);% set dof=6, and verbose=1
%    See also: flirt_wrapper 

% AUTHOR    : Yi Sui
% DATE      : 06/26/2017
%%
if ~exist('ssh2_conn','var'), ssh2_conn = []; end
if ~exist('options','var'), options =''; end
if ~exist('verbose','var'), verbose =false; end

tmp_dir = [tempname,'_tmp']; %get temp dir name
mkdir(tmp_dir);
tmp_inputvol = fullfile(tmp_dir, 'inputvol.nii');
tmp_refvol = fullfile(tmp_dir,'refvol.nii');
tmp_outputvol = fullfile(tmp_dir,'outputvol');
tmp_omat = fullfile(tmp_dir,'omat.mat');

%save to nii file
nii=make_nii(inputvol,inputvoxsize);
save_nii(nii,tmp_inputvol);

nii=make_nii(refvol,refvoxsize);
save_nii(nii,tmp_refvol);

addition_args={};
if nargout >1 %output affine matrix
    addition_args=[addition_args,{'omat',tmp_omat}];
end

try
    [ssh2_conn,command_result]=flirt_wrapper(tmp_inputvol,tmp_refvol,tmp_outputvol,'ssh2_conn',ssh2_conn,'option',options,addition_args{:},'verbose',verbose);
    
    %read in .nii files.
    nii = load_nii( [tmp_outputvol,'.nii.gz']);
    regvol = nii.img;
    
    if nargout >1
        omat = reshape(textread( tmp_omat,'%f'),[4,4])';
    end
    
catch err
    ssh2_conn = ssh2_close(ssh2_conn);
    %delete the tmp files on local PC
    if isempty( strfind(tmp_dir,'*'))&& ~isempty( strfind(tmp_dir,'_tmp'))
        rmdir(tmp_dir,'s');
    end
    rethrow(err);
end

if isempty( strfind(tmp_dir,'*')) && ~isempty( strfind(tmp_dir,'_tmp'))
    disp(tmp_dir)
    rmdir(tmp_dir,'s');
end
