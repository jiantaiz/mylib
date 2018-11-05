function M = unwrap3D_ssh(mat,DO_MEAN_SHIFT,ssh2_conn,print_output)
% unwrap3D_ssh remotely runs Josh's unix program GC_UNWRAP on mre-cim (graph cut unwrap)
%    M = unwrap3D_ssh(mat, DO_MEAN_SHIFT,username,password,print_output)
%    See also: unwrap3D

% AUTHOR    : Yi Sui
% DATE      : 05/16/2017

host = 'mr-cim';
if nargin<2
    DO_MEAN_SHIFT = 0;
end
if nargin<3 || isempty(ssh2_conn)
    
    private_key = 'C:\Users\m165355\OneDrive\Mayo\Matlab\mr_cim_id_rsa';
    if exist(private_key,'file') == 2
        ssh2_conn = ssh2_config_publickey('mr-cim', 'm165355',private_key, '');
    else
        ssh2_conn=ssh2_login('m165355','','mr-cim');
    end
end
    
if isempty(ssh2_conn) % Cancelled
    disp('Login cancelled!!');
    M=-1;
    return ;
end
pause(0.1);
if nargin<4
    print_output = 0;
end

tmp_file_in = tempname();
[tmp_file_in_pathstr,tmp_file_in_name,ext]=fileparts(tmp_file_in);
tmp_file_out = tempname();
[pathstr,tmp_file_out_name,ext]=fileparts(tmp_file_out);

%save to binary file
mat(isnan(mat))=0;
mat2bin('',tmp_file_in,mat);

% GC_UNWRAP IN_FILENAME OUT_FILENAME Nx Ny Nz Nt DO_MEAN_SHIFT
[Nx, Ny, Nz, Nt]=size(mat);
cmd = sprintf('~/bin/GC_UNWRAP %s %s %d %d %d %d %d > GC_UNWRAP.output& ',tmp_file_in_name,tmp_file_out_name,Nx, Ny, Nz, Nt, DO_MEAN_SHIFT );


try
    %config ssh connection to server
%     ssh2_conn = ssh2_config(host,username,password);
    %send file to server
    ssh2_conn = scp_put(ssh2_conn,tmp_file_in_name,'',tmp_file_in_pathstr);
    %execute cmd on server
    [ssh2_conn COMMAND_RESULT]= ssh2_command(ssh2_conn,cmd,1);
        
    % waiting for the job done
    fprintf('3D unwrapping on %s (%d ts). ',host,Nt);
    cnt = 0;
    
    % unicode
    %     9474?
    %     9585?
    %     9472?
    %     9586?
    %     9587?
    %     9532?
    pinwheel =[9474 9585 9472 9586 ];
    
    % check if GC_UNWRAP process exists
    [ssh2_conn, r]= ssh2_command(ssh2_conn,'pgrep GC_UNWRAP');
    tic
    last_step = -1;
    while  ~isempty(r{1})
        if (mod(cnt,10) == 0 )
            %              fprintf('%c',bagua(mod(cnt,8)+1 ));% 8=backspace
            
            try
                %grep GC_UNWRAP output to see the progress.
                [ssh2_conn, txt]= ssh2_command(ssh2_conn,'grep -a unwrapping GC_UNWRAP.output | tail -1');
                cur_step = sscanf(txt{1},'%*s %*s %*s %*s %d %*s %*s %d');
                if cur_step(1) == last_step
                    fprintf('%c.%c',8,pinwheel(mod(cnt,4)+1 ));% 8=backspace
                else
                    fprintf('%c%d.%c',8,cur_step(1),pinwheel(mod(cnt,4)+1 ));% 8=backspace
                    last_step = cur_step(1);
                end
            end
            
            % check if GC_UNWRAP is still running
            [ssh2_conn, r]= ssh2_command(ssh2_conn,'pgrep GC_UNWRAP');
            
        else
            % spin the pinwheel
            fprintf('%c%c',8,pinwheel(mod(cnt,4)+1 ));% 8=backspace
            %              fprintf('%c%c',8,bagua(mod(cnt,8)+1 ));% 8=backspace
            
            pause(2/10);
        end
        cnt = cnt+1;
    end
    run_time = toc;
    fprintf('%c',8);% 8=backspace
    fprintf(' done! (%.1f sec)',run_time);
    
    if print_output
        [ssh2_conn, ~]= ssh2_command(ssh2_conn,'cat GC_UNWRAP.output',1);
    end
    
    %retrieve file from server
    ssh2_conn = scp_get(ssh2_conn,tmp_file_out_name,tmp_file_in_pathstr);
    %delete the tmp files on server
    if isempty( strfind(tmp_file_in_name,'*')) % to be safe, check no '*' in file name
        ssh2_conn = ssh2_command(ssh2_conn,sprintf('rm %s',tmp_file_in_name));
    end
    if isempty( strfind(tmp_file_out_name,'*')) % to be safe, check no '*' in file name
        ssh2_conn =  ssh2_command(ssh2_conn,sprintf('rm %s',tmp_file_out_name));
    end
    %close connection
    ssh2_conn = ssh2_close(ssh2_conn);
catch err
    ssh2_conn = ssh2_close(ssh2_conn);
    delete(tmp_file_in,tmp_file_out);
    rethrow(err);
end

dims=size(mat);
fid = fopen(tmp_file_out,'r');
M = fread(fid,'single');
M = single(reshape(M,dims));
fclose(fid);

%delete the tmp files on local PC
if isempty( strfind(tmp_file_in_name,'*')) && isempty( strfind(tmp_file_out_name,'*'))
    delete(tmp_file_in,tmp_file_out);
end

