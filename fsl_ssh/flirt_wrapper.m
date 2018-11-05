function [ssh2_conn, command_result] = flirt_wrapper(varargin)
% FLIRT_WRAPPER FLIRT wrapper
%    FLIRT_WRAPPER(inputvol,refvol,outputvol)
%    FLIRT_WRAPPER(...,'omat', <matrix-filname>) saves affine matrix file
%    FLIRT_WRAPPER(...,'applyxfm', [0|1], 'init', <matrix-filname> ) applies affine transform
%    FLIRT_WRAPPER(...,'option',<flirt options>, 'ssh2_conn', <ssh2_conn>, 'verbose', [0|1])
%    []FLIRT_WRAPPER(...,'run_in_bg', [0|1])  run command in background via parfeval()
%    [ssh2_conn, command_result] = FLIRT_WRAPPER(...) outputs ssh2 connection and command results
%
%    Examples:
%    FLIRT_WRAPPER('example/epi.nii','example/T2_2mm.nii','example/epi_to_T2_2mm');
%
%    ssh2_conn  = ssh2_login;
%    [ssh2_conn,command_result]=FLIRT_WRAPPER('example/epi.nii','example/T2_2mm.nii','example/epi_to_T2_2mm','verbose',1,'omat','example/epi_to_T2_2mm.mat','ssh2_conn',ssh2_conn);
%    
% 
% Usage: flirt [options] -in <inputvol> -ref <refvol> -out <outputvol>
%        flirt [options] -in <inputvol> -ref <refvol> -omat <outputmatrix>
%        flirt [options] -in <inputvol> -ref <refvol> -applyxfm -init <matrix> -out <outputvol>
% 
%   Available options are:
%         -in  <inputvol>                    (no default)
%         -ref <refvol>                      (no default)
%         -init <matrix-filname>             (input 4x4 affine matrix)
%         -omat <matrix-filename>            (output in 4x4 ascii format)
%         -out, -o <outputvol>               (default is none)
%         -datatype {char,short,int,float,double}                    (force output data type)
%         -cost {mutualinfo,corratio,normcorr,normmi,leastsq,labeldiff,bbr}        (default is corratio)
%         -searchcost {mutualinfo,corratio,normcorr,normmi,leastsq,labeldiff,bbr}  (default is corratio)
%         -usesqform                         (initialise using appropriate sform or qform)
%         -displayinit                       (display initial matrix)
%         -anglerep {quaternion,euler}       (default is euler)
%         -interp {trilinear,nearestneighbour,sinc,spline}  (final interpolation: def - trilinear)
%         -sincwidth <full-width in voxels>  (default is 7)
%         -sincwindow {rectangular,hanning,blackman}
%         -bins <number of histogram bins>   (default is 256)
%         -dof  <number of transform dofs>   (default is 12)
%         -noresample                        (do not change input sampling)
%         -forcescaling                      (force rescaling even for low-res images)
%         -minsampling <vox_dim>             (set minimum voxel dimension for sampling (in mm))
%         -applyxfm                          (applies transform (no optimisation) - requires -init)
%         -applyisoxfm <scale>               (as applyxfm but forces isotropic resampling)
%         -paddingsize <number of voxels>    (for applyxfm: interpolates outside image by size)
%         -searchrx <min_angle> <max_angle>  (angles in degrees: default is -90 90)
%         -searchry <min_angle> <max_angle>  (angles in degrees: default is -90 90)
%         -searchrz <min_angle> <max_angle>  (angles in degrees: default is -90 90)
%         -nosearch                          (sets all angular search ranges to 0 0)
%         -coarsesearch <delta_angle>        (angle in degrees: default is 60)
%         -finesearch <delta_angle>          (angle in degrees: default is 18)
%         -schedule <schedule-file>          (replaces default schedule)
%         -refweight <volume>                (use weights for reference volume)
%         -inweight <volume>                 (use weights for input volume)
%         -wmseg <volume>                    (white matter segmentation volume needed by BBR cost function)
%         -wmcoords <text matrix>            (white matter boundary coordinates for BBR cost function)
%         -wmnorms <text matrix>             (white matter boundary normals for BBR cost function)
%         -fieldmap <volume>                 (fieldmap image in rads/s - must be already registered to the reference image)
%         -fieldmapmask <volume>             (mask for fieldmap image)
%         -pedir <index>                     (phase encode direction of EPI - 1/2/3=x/y/z & -1/-2/-3=-x/-y/-z)
%         -echospacing <value>               (value of EPI echo spacing - units of seconds)
%         -bbrtype <value>                   (type of bbr cost function: signed [default], global_abs, local_abs)
%         -bbrslope <value>                  (value of bbr slope)
%         -setbackground <value>             (use specified background value for points outside FOV)
%         -noclamp                           (do not use intensity clamping)
%         -noresampblur                      (do not use blurring on downsampling)
%         -2D                                (use 2D rigid body mode - ignores dof)
%         -verbose <num>                     (0 is least and default)
%         -v                                 (same as -verbose 1)
%         -i                                 (pauses at each stage: default is off)
%         -version                           (prints version number)
%         -help
%    See also: ssh_flirt, bet_wrapper, fsl_wrapper


% AUTHOR    : Yi Sui
% DATE      : 06/26/2017
%%
OptionsParser = inputParser;
OptionsParser.KeepUnmatched = true;
OptionsParser.addRequired('inputvol',@ischar);
OptionsParser.addRequired('refvol',@ischar);
OptionsParser.addRequired('outputvol',@ischar);
OptionsParser.addParameter('omat','',@ischar);
OptionsParser.addParameter('applyxfm',0);
OptionsParser.addParameter('init','',@ischar);
OptionsParser.addParameter('option','');
OptionsParser.addParameter('verbose',0);
OptionsParser.addParameter('ssh2_conn',[]);
OptionsParser.addParameter('run_in_bg',0);

% Parse Input Arguments
try
    OptionsParser.parse(varargin{:});
catch Error
    OptionsParser.parse;
    if strcmpi(Error.identifier, 'MATLAB:InputParser:ArgumentFailedValidation')
        error(Error.identifier, Error.message);
    end;
end;
Options = OptionsParser.Results;

[inputvol.pathstr,inputvol.name,inputvol.ext]  = fileparts(Options.inputvol);
infiles = {Options.inputvol};

[refvol.pathstr,refvol.name,refvol.ext]  = fileparts(Options.refvol);
if ~strcmpi( refvol.pathstr ,'fsl_template')
    infiles {end+1} = Options.refvol;
else
    refvol.name = ['$FSLDIR/data/standard/', refvol.name, refvol.ext];
end

[outputvol.pathstr,outputvol.name,outputvol.ext]  = fileparts(Options.outputvol);
outfiles = {[Options.outputvol,'*']};

cmd = sprintf('flirt -in %s -ref %s -out %s', inputvol.name, refvol.name,outputvol.name);
if ~isempty(Options.omat)
    [omat.pathstr,omat.name,omat.ext]  = fileparts(Options.omat);
    outfiles{end+1} = Options.omat;
    cmd = sprintf('%s -omat %s',cmd,[omat.name,omat.ext]);
end

if Options.applyxfm
    [init.pathstr,init.name,init.ext]  = fileparts(Options.init);
    infiles{end+1} = Options.init;
    cmd = sprintf('%s -applyxfm -init %s',cmd, [init.name,init.ext]);
end

cmd = [cmd ' ' Options.option];
[ssh2_conn, command_result] = fsl_wrapper(cmd,infiles,outfiles,Options.ssh2_conn,Options.verbose,Options.run_in_bg);
