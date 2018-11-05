function [job,t] = bgbatch(cmd)
% BGRUN running command in background
%    [job,t] = BGRUN(cmd) running command string in background, cmd is a
%    string
%
%    Example:
%    [job,t] = bgrun('[ brain,brain_mask skull_mask ,scalp_mask ] = bet_ssh(M,[240/256,240/256,2],ssh)');
%
%    Subfunctions: timer_fcn
%    See also: batch

% AUTHOR    : Yi Sui
% DATE      : 07/06/2017
%%
%parse output arguments
[str, remain] = strtok(cmd, '=');
argout={};
if ~isempty(remain)%with output
    delimiter = ' ,[]';
    [str, remain] = strtok(str, delimiter);
    argout(end+1) = {str};
    while true
        [str, remain] = strtok(remain, delimiter);
        if isempty(str),  break;  end
        argout(end+1) = {str};
    end
end

cmd = strrep(cmd,'''','''''');

expr = sprintf('batch(''%s'',''Profile'', ''local'', ''Pool'', 0,''CaptureDiary'',false);',cmd);
warning off
job = evalin('base',expr);
warning on
disp('Job is running in backgroud.')

%set up a timer to monitor the job
t = timer;
t.TimerFcn = @timer_fcn;
t.Period = 1;
t.TasksToExecute = inf;
t.ExecutionMode = 'fixedRate';

ud.job = job;
ud.argout = argout;
t.UserData = ud;
start(t);


function timer_fcn(mTimer,~)
% mTimer
ud = mTimer.UserData;
jb = ud.job;
% disp(jb.State)
msg='';
if strcmpi(jb.State,'finished');
    if ~isempty(ud.argout)
        S = load(jb,ud.argout{:});
        fields = fieldnames(S);
        for k=1:numel(fields)
            assignin('base',fields{k},S.(fields{k}));
        end
        msg = sprintf('%s, ',fields{:});
        msg = sprintf('%s saved in the work space!!',msg);
    end
    stop(mTimer);
    delete(mTimer);

    msg = sprintf('%s\n%s\n %s\n',jb.Name,'Job finished!!',msg);
    disp(msg);
    msgbox(msg,'bgrun');
    
end
% disp(mTimer.UserData.State);
