function M = bin2mat(filepath,filename,dims)
% % 
% fid = fopen(fullfile(filepath,filename),'r');
% M = fread(fid,'single'); % single
% fclose(fid);
% M = reshape(M,dims);
%%
%
% filepath = 'E:\ZY\Sugarsync\MRE Data\Brain Data\meningioma\09-096-835\';
cd(uigetdir(filepath)); 
folder = dir(pwd);
folder = folder(~strncmpi('.', {folder.name}, 1)); % folders and files
% dims=[128 128 48 8];
% files1_to_look_for = {'*BMF1.0_all.phs'}; 
% files1_to_look_for = {'*3duw.phs'}; 
files1_to_look_for = {filename};
N_file1_to_look_for = length(files1_to_look_for);


n1 = 1; n2 = 1; n3 = 1; n4 = 1;

dir_list = {pwd}; Ndirs = length(dir_list); done_looking = 0; dir_cnt = 1; file_cnt = 0;
while (done_looking == 0)
    if (dir_cnt > Ndirs), break; end
    cur_dir = dir_list{dir_cnt};
    for k_file_search = 1:N_file1_to_look_for % 8 offsets
        D = dir(fullfile(cur_dir,files1_to_look_for{k_file_search}));
        if (isempty(D)), continue; end
        N_files_found = length(D);
        for k_files = 1:N_files_found
            cur_file = D(k_files).name; if (strcmp(cur_file,'.') || strcmp(cur_file,'..')), continue; end
            if (D(k_files).isdir == 1), continue; end
            disp(fullfile(cur_dir,cur_file));
            fid = fopen(fullfile(cur_dir,cur_file),'r');
            M = fread(fid,'single');
            M = single(reshape(M,dims));
            fclose(fid); 
            save(fullfile(cur_dir,[cur_file '.mat']),'M');
            disp(fullfile(cur_dir,cur_file));
            file_cnt = file_cnt + 1;
        end
    end
    D = dir(fullfile(cur_dir,'*')); N_obj = length(D);
    for k_obj = 1:N_obj
        cur_obj = D(k_obj).name; if (strcmp(cur_obj,'.') || strcmp(cur_obj,'..')), continue; end
        if (D(k_obj).isdir == 0), continue; end
        dir_list = cat(1,dir_list,[cur_dir,'\',cur_obj,'\']); Ndirs = Ndirs+1;
    end
    dir_cnt = dir_cnt+1;
    if (dir_cnt > Ndirs), done_looking = 1; end
end