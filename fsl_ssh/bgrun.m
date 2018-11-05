function job = bgrun(varargin)
% BGRUN run command in background
%    BGRUN(cmd) cmd is a command string
%    job = BGRUN(cmd) return the job handle
%    you can add BGRUN in front of any command to put it in background
%
%    Example:
%    BGRUN a = rand(3)
%    job = BGRUN ('a = rand(3)')
%
%    See also: parfeval

% AUTHOR    : Yi Sui
% DATE      : 07/13/2017
%%
cmd = [varargin{:}];
[argout,funcname,argin] = parse_cmd_str(cmd);
for k = 1: numel(argin)
    argin{k} = evalin('base',argin{k});
end
job = parfeval(funcname,numel(argout),argin{:});
fprintf('%s is running in backgroud.\n',job.Function);

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
job = ud.job;
% disp(jb.State)
msg='';
if strcmpi(job.State,'finished');
    
    if ~isempty(ud.argout) &&  isempty(job.Error)
        out = cell(1,numel(ud.argout));
        [out{:}] = job.fetchOutputs;
        for k=1:numel(out)
            assignin('base',ud.argout{k},out{k});
        end
        msg = sprintf('%s, ',ud.argout{:});
        if numel(msg) > 2
            msg=msg(1:end-2); %remove ',' and space at the end
        end
        
    end
    stop(mTimer);
    delete(mTimer);
    if isempty(msg)
        msg = sprintf('%s finished!!\n',job.Function);
    else
        msg = sprintf('%s finished!!\nThe output(s): ''%s'' saved in work space\n',job.Function,msg);
    end
    fprintf(msg);
    if ~isempty(job.Error)
        fprintf('ERROR in %s: %s\n',job.Function, job.Error.message);
        disp(job)
    end
    %     msgbox(msg,'bgrun');
    
end
% disp(mTimer.UserData.State);

function [argout,funcname,argin]=parse_cmd_str(cmd)
[outputstr, bodystr] = strtok(cmd, '=');
argout={};
if ~isempty(bodystr)%with output
    delimiter = ' ,[]';
    [outputstr, remain] = strtok(outputstr, delimiter);
    argout(end+1) = {outputstr};
    while true
        [outputstr, remain] = strtok(remain, delimiter);
        if isempty(outputstr),  break;  end
        argout(end+1) = {strtrim(outputstr)};
    end
else % no output
    bodystr = outputstr;
end

[funcname, inputstr] = strtok(bodystr, '=()');
funcname = strtrim(funcname);

argin={};
if ~isempty(inputstr)%with input
    delimiter = ' ,()';
    [inputstr, remain] = strtok(inputstr, delimiter);
    argin(end+1) = {inputstr};
    while true
        [inputstr, remain] = strtok(remain, delimiter);
        if isempty(inputstr),  break;  end
        argin(end+1) = {strtrim(inputstr)};
    end
end


