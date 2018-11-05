function [ssh2_conn, COMMAND_RESULT] = fsl_bet(infile, outfile, varargin)
% FSL_BET run FSL bet on remote server, this command can use the same
% syntax as in Linux, see example below
%    [ssh2_conn, COMMAND_RESULT] = FSL_BET(infile, outfile, varargin)
%
%    Example:
%    % run bet in background by giving a '&' at the end
%    fsl_bet magn.nii magn_brain -m -f 0.6 -R -A & ;
%
%    See also: fsl_wrapper, bet_wrapper, bet_ssh

% AUTHOR    : Yi Sui
% DATE      : 07/07/2017
%%
global fsl_ssh2_conn
opt = sprintf('%s ',varargin{:});
opt=strtrim(opt);
if opt(end) == '&'
    run_in_bg = 1;
    opt=opt(1:end-1);
else
    run_in_bg = 0;
end

if isstruct(fsl_ssh2_conn)
    ssh2_conn=fsl_ssh2_conn;
else
    warning('Tip: you can create a global variable ''fsl_ssh2_conn'' for ssh2 connection.')
    ssh2_conn=[];
end

verbose = 0;
[ssh2_conn, COMMAND_RESULT] = bet_wrapper(infile, outfile, opt, ssh2_conn,verbose,run_in_bg);
