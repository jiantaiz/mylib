function move_cal_h5_file(root)
if nargin<1
    root='.';
end
cals=fullfile(root,'CALS');
mkdir(cals);

if exist(cals,'dir')
    f=rdir(fullfile(root,'*','CALS','*.h5'));
    for k=1:numel(f);
        movefile(f(k).name,cals);
    end
    
    d=rdir(fullfile(root,'*','CALS*'));
    for k=1:numel(d);
        rmdir(d(k).name)
    end
else
    error('CALS folder not created successfully');
end