function se = find_series_by_desc (pat_dir,pattern)
% FIND_SERIES_BY_DESC returns series ID by its descriptions
%    se = FIND_SERIES_BY_DESC (pat_dir,pattern) searches "pat_dir" folder and
%    find series with descrition having the regexp "pattern".
%
%    Example:
%    se = FIND_SERIES_BY_DESC ('.','^CMRE.*') will search current folder and 
%    returns the series ID with description beginning with 'CMRE'.
%
%    See also: 

% AUTHOR    : Yi Sui
% DATE      : 05/16/2017

se=[];
folders = dir(fullfile(pat_dir,'s*'));
tic
parfor d = 1:numel(folders)
    if folders(d).isdir
        se_folder = fullfile(pat_dir, folders(d).name);
        dcm = dir(fullfile(se_folder,'*.dcm'));
        dcm = dcm(~[dcm.isdir]);
%         fullfile(se_folder,dcm(1).name)
        dinfo = dicominfo(fullfile(se_folder,dcm(1).name));
        se_desc = dinfo.SeriesDescription;
        if ~isempty(regexpi(se_desc,pattern))
            se = [se dinfo.SeriesNumber];
        end
    end
end
toc
