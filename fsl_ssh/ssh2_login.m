function ssh2_conn = ssh2_login(username, password, host)
% SSH2_LOGIN displays Login dialog and create a ssh2 config structure.
%    ssh2_conn = SSH2_LOGIN displays Login dialog and returns a ssh2 config
%    structure ssh2_conn for further use.
%    ssh2_conn = SSH2_LOGIN(username, password, host) specifies the default
%    value for username, password and host.
%
%    See also: logindlg2, ssh2_config

% AUTHOR    : Yi Sui
% DATE      : 06/26/2017
%%

if ~exist('username','var') || isempty(username)
    username = '';
end
if ~exist('password','var') || isempty(password)
    password = '';
end
if ~exist('host','var') || isempty(host)
    host = '';
end

[password,username,host] = logindlg2('enterUserName', true,'DefaultUserName',username,'enterHostName', true,'DefaultHostName',host,'WindowName',host);
pause(0.001);
if password ==-1 %cancelled
    ssh2_conn=[];
    return;
end

ssh2_conn = ssh2_config(host,username,password);

end
