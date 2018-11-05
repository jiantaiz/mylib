
function [ results, cancel ] = Standard_MRE_Pipeline(directory, directory_nm, t1_directory, rlist, DI_option, ...
                                                     no_results, reg_reslice_only, n_erode_roi, save_cmu_files, ...
                                                     select_roi_n, sub_directory)
%results = Standard_MRE_Pipeline
%          (
%                      directory,   %     DICOM "motion" file directory (NOTE 1)
%                   directory_nm,   %     DICOM "no-motion" file directory
%                   t1_directory,   %     DICOM "T1" file directory
%                          rlist,   % [1:20,22,37] Analysis region range
%                      DI_option,   % [1] 1 = "Normal" Direct Inversion (amp^2 weighting of XYZ)
%                                         2 = Josh Trzasko's least-squares(XYZ) Direct Inversion
%                                         3 = Josh Trzasko's minimum absolute deviation (MAD) DI
%                     no_results,   % [0] 1 = do not create .csv and .mat results files
%                                             (and return empty 'results' from this function)
%               reg_reslice_only,   % [0] 1 = register and reslice data sets, but do not process
%                    n_erode_roi,   % [0] NET erosion of (user-defined) ROI regions (0 or 1 only);
%                                        (NOTE that regular brain atlas regions are always eroded by 1)
%                 save_cmu_files,   % [0] 1 = save binary files of complex-mu images (each region)
%                   select_roi_n,   % [0] 0  = interactive user selection of each ROI file
%                                         1+ = select nth newest (modified date) ROI file, if exists
%                                        (1 = newest, 2 = next newest, etc.)
%                  sub_directory    %  [] New or existing sub-directory (of folder pipeline_results) for all results
%          )
%
%(NOTE 1: If ONLY 'directory' is supplied, it will be used as the 
%starting directory for interactive (user) selection of the 'motion', 
%'no-motion' and 'T1' directories.)
%
%Modified version of Matt Murphy's original "Brain MRE pipeline" script:
%David Lake, beginning Jan 2013.

%Modified version of Brain_MRE_Pipeline.m (3-mm pixels), to process 2.5-mm data.
%(Should process "any" given input data, AT ITS NATIVE ACQUISITION RESOLUTION;
% slice- (z-) dimension is NOT interpolated.)
%David Lake, beginning Apr 21, 2015.

%This version is the "standard pipeline" BEFORE proposed revision of 
%pipeline processing (median-filtering of full-volume DI result, followed 
%by post-masking of desired regions).  That is, this version still uses 
%Matt's edge-adaptive DI and curl filters to reduce edge effects.

clearvars -except directory directory_nm t1_directory rlist DI_option no_results reg_reslice_only ...
                  n_erode_roi save_cmu_files select_roi_n sub_directory
              
% Reset Atlas regions for each new run...              
RESET = 1;
define_region(0, 0, 0, 0, RESET);

% Set special CSF fraction limit for masking pixels
% (Set to 0 for normal operation)
special_csf_fraction = 0;

temp = 0;
CALC_CMEDIAN_MAG = temp;    % 1 = calculate magnitude of iterative-complex-median 
                            % of complex-mu, and include result in spreadsheet.
         
temp = 0;                            
CALC_MAG_FROM_RI_MEANS = temp; % 1 = calculate magnitude of complex-mu, from
                               % independent MEANS of Real/Imag parts

% Define default parameter values here...
DEFAULT_rlist      =  [1:20,22,37];
DEFAULT_rlist_str  = '[1:20,22,37]';
DEFAULT_DI_option         = 1;
DEFAULT_no_results        = 0;
DEFAULT_reg_reslice_only  = 0;
DEFAULT_n_erode_roi       = 0;
DEFAULT_save_cmu_files    = 0;
DEFAULT_select_roi_n = 0;

% try
    
    display_selected_dirs = 1;
    f = 0; %#ok<NASGU>

    if (~exist('sub_directory','var') || length(sub_directory) < 1)
        sub_directory = [];
    end
    
    if ~exist('DI_option','var')
        DI_option = DEFAULT_DI_option;
    end
    
    if ~exist('no_results','var')
       no_results = DEFAULT_no_results; 
    end
    
    if ~exist('reg_reslice_only','var')
       reg_reslice_only = DEFAULT_reg_reslice_only; 
    end

    if ~exist('n_erode_roi','var')
        n_erode_roi = DEFAULT_n_erode_roi; 
    end
    
    if ~exist('save_cmu_files','var')
       save_cmu_files = DEFAULT_save_cmu_files; 
    end

    if ~exist('select_roi_n','var')
       select_roi_n = DEFAULT_select_roi_n; 
    end
    
    if (n_erode_roi < 0 || 1 < n_erode_roi)
       fprintf('\n### ERROR: n_erode_roi may be only 0 or 1 ###\n\n');
       return;
    end
    
    if ~exist('directory_nm','var')

        pause_flag = 1;
        if (exist('directory','var'))
            %%%fprintf('%s\n', directory);
            
            default_dir = dir_fullfile(fullfile(directory, '*_motion'));
            if (~isempty(default_dir) && length(default_dir) == 1)
                directory = uigetdir(default_dir.name, 'Select "motion" data directory...');
            else
                directory = uigetdir(directory, 'Select "motion" data directory...');
            end
        else
            directory = uigetdir('/mnt/huston_data', 'Select "motion" data directory...');
            if (~isempty(directory))
                fprintf('%s%s\n', '                      ', directory);
            end
        end
        
        if (directory == 0)
            fprintf('(Cancelled)\n\n');
            return
        else
            % Simple directory validity check (contains DICOM file)...
            dfilename = fullfile(directory, 'i001.dcm');
            if (exist(dfilename, 'file') ~= 2)
                fprintf('\nDirectory %s\n', directory);
                fprintf('is not a valid "motion" data directory.\n');
                fprintf('(Processing Cancelled)\n\n');
                return
            end
        end
        
        directory = fullfile(directory, filesep);   % MRE images with motion

        motion_imgInfo = dicominfo(fullfile(directory, 'i001.dcm'));
        
        if (ispc)
            directory_nm = uigetdir(fullfile(directory,''), 'Select "NO motion" data directory...');
        else
            default_dir = dir_fullfile(fullfile(directory, '..', '*_NOmotion'));
            if (~isempty(default_dir) && length(default_dir) == 1)
                directory_nm = uigetdir(default_dir.name, 'Select "NO motion" data directory...');
            else
                directory_nm = uigetdir(fullfile(directory,'..'), 'Select "NO motion" data directory...');
            end
        end
        
        if (directory_nm == 0)
            fprintf('(Cancelled)\n\n');
            return
        else
            % Simple directory validity check (contains DICOM file)...
            dfilename = fullfile(directory_nm, 'i001.dcm');
            if (exist(dfilename, 'file') ~= 2)
                fprintf('\nDirectory %s\n', directory_nm);
                fprintf('is not a valid "NO motion" data directory.\n');
                fprintf('(Processing Cancelled)\n\n');
                return
            end
        end
        
        directory_nm = fullfile(directory_nm, filesep); % MRE images with NO motion

        if (ispc)
            t1_directory = uigetdir(fullfile(directory_nm,''), 'Select T1 data directory...');
        else
            default_dir = dir_fullfile(fullfile(directory, '..', '*_T1'));
            if (~isempty(default_dir) && length(default_dir) == 1)
                t1_directory = uigetdir(default_dir.name, 'Select T1 data directory...');
            else
                t1_directory = uigetdir(fullfile(directory_nm,'..'), 'Select T1 data directory...');
            end
        end
        
        if (t1_directory == 0)
            fprintf('(Cancelled)\n\n');
            return
        else
            % Simple directory validity check (contains DICOM file)...
            dfilename = fullfile(t1_directory, 'i001.dcm');
            if (exist(dfilename, 'file') ~= 2)
                fprintf('\nDirectory %s\n', t1_directory);
                fprintf('is not a valid T1 data directory.\n');
                fprintf('(Processing Cancelled)\n\n');
                return
            end
        end
        
        t1_directory = fullfile(t1_directory, filesep); % T1-weighted images

        t1_imgInfo = dicominfo(fullfile(t1_directory, 'i001.dcm')); % get dicom info from the T1 image, this will let us reconstruct the names of the atlases, masks, segmentation images, etc. that come from grinder/ATLAS toolbox
        
        try
            f = motion_imgInfo.Private_0019_10bc;
            if (f < 0)
                f = motion_imgInfo.Private_0019_10b4;
            end
        catch
            f = 0;
        end
        
        if (f <= 0)
            % Get driver frequency from user...
            f = 60;   % Default driver frequency, Hz
            instr = input('\nEnter driver frequency, Hz (<Enter> for 60Hz) : ','s');
            if (~isempty(instr))
                fin = str2double(instr);
                if (fin > 0)
                    f = fin;
                end
            end
        end
        
        % Display selected directories...
        disp( '---------------------');
        disp(['    Driver frequency: ', num2str(f), ' Hz']);
        try
            proc_zspacing_mm = single(motion_imgInfo.SpacingBetweenSlices);
            proc_xspacing_mm = proc_zspacing_mm;
            proc_yspacing_mm = proc_zspacing_mm;
            if (abs(proc_xspacing_mm - 3.0) > 0.01 || abs(proc_yspacing_mm - 3.0) > 0.01 || abs(proc_zspacing_mm - 3.0) > 0.01)
                fprintf('### Note pixel spacing for processing: %.3f, %.3f, %.3f mm ###\n', proc_xspacing_mm, proc_yspacing_mm, proc_zspacing_mm);
            end
        catch
        end
        disp(['    Motion directory: ', directory]);
        disp([' No-motion directory: ', directory_nm]);
        disp(['        T1 directory: ', t1_directory]);
        disp( '---------------------');
        display_selected_dirs = 0;  % No need to display again (below)...

    else
        motion_imgInfo = dicominfo(fullfile(directory, 'i001.dcm'));
        try
            f = motion_imgInfo.Private_0019_10bc;
            if (f < 0)
                f = motion_imgInfo.Private_0019_10b4;
            end
        catch
            f = 0;
        end
        
        t1_imgInfo = dicominfo(fullfile(t1_directory, 'i001.dcm')); % get dicom info from the T1 image, this will let us reconstruct the names of the atlases, masks, segmentation images, etc. that come from grinder/ATLAS toolbox
        pause_flag = 0;
    end

    if ~exist('rlist','var')
        % Get region range from user...
        rlist = [];
        while (isempty(rlist))
            try
                instr = input(['Enter regions (<Enter> for ' DEFAULT_rlist_str ', ''help'' for list): '],'s');
                if (length(instr) < 1)
                    rlist = DEFAULT_rlist;
                    disp('Please wait... ');
                elseif (strcmp(instr, 'help'))
                    % List valid regions...
                    fprintf('\n');
                    for rr = 1:37
                        fprintf('   %2d %s\n', rr, region_name(rr));
                        if (rr == 20 || rr == 21 || rr == 37)
                            fprintf('\n');
                        end
                    end
                    fprintf('   41 ROI (user will select .nii file)\n');
                    fprintf('\n');
                    fprintf('   51 cerebrum MINUS whole ROI/tumor (user will select .nii file)\n');
                    fprintf('\n');
                    rlist = [];
                else
                    rlist = eval(instr);
                end
            catch %#ok<*CTCH>
                rlist = [];
            end
        end
    end
    
    % Determine acquisition/DICOM-data/processing parameters...
    
    % Note that acquired X/Y sizes are often different than the DICOM X/Y sizes 
    % of each slice; the scanner often interpolates reconstructed slices to a 
    % common FFT size (e.g., 128x128 pixels).
    % =========================================
    
    % Determine acquisition X- and Y-sizes...
    acq_nx = -1;
    for am_index = 1:size(motion_imgInfo.AcquisitionMatrix)
        if (motion_imgInfo.AcquisitionMatrix(am_index) > 0)
            if (acq_nx == -1)
                % First non-zero value in AcquisitionMatrix...
                acq_nx = motion_imgInfo.AcquisitionMatrix(am_index);
            else
                if (motion_imgInfo.AcquisitionMatrix(am_index) ~= acq_nx)
                    acq_nx = -1;
                end
            end
        end
    end

    if (acq_nx <= 0)
        fprintf('\n### Invalid AcquisitionMatrix size(s): %d, %d, %d, %d ###', ...
                motion_imgInfo.AcquisitionMatrix(1), motion_imgInfo.AcquisitionMatrix(2), ...
                motion_imgInfo.AcquisitionMatrix(3), motion_imgInfo.AcquisitionMatrix(4));
        results = [];
        cancel = 0;
        return
    end
    
    % Note: acquisition X- and Y-sizes are equal...
%     acq_ny = acq_nx;
%     acq_nz = motion_imgInfo.Private_0021_104f;
    acq_nt = motion_imgInfo.NumberOfTemporalPositions;
    
    % Determine acquisition FOV (same in X and Y)...
    acq_fovx_mm = motion_imgInfo.ReconstructionDiameter;
    acq_fovy_mm = motion_imgInfo.ReconstructionDiameter;
    
%     % Determine acquisition pixel spacings...
%     acq_xspacing_mm = acq_fovx_mm / acq_nx;
%     acq_yspacing_mm = acq_fovy_mm / acq_nx;
%     acq_zspacing_mm = motion_imgInfo.SpacingBetweenSlices;
    
    % Determine related parameters of DICOM data...
    dicom_nx = motion_imgInfo.Columns;
    dicom_ny = motion_imgInfo.Rows;
    dicom_nz = motion_imgInfo.Private_0021_104f;
    
%     dicom_xspacing_mm = motion_imgInfo.PixelSpacing(1);
%     dicom_yspacing_mm = motion_imgInfo.PixelSpacing(2);
%     dicom_zspacing_mm = motion_imgInfo.SpacingBetweenSlices;

    % Calculate related parameters for processing...
    proc_zspacing_mm = single(motion_imgInfo.SpacingBetweenSlices);
    
    % Interpolate data (in X and Y) to isotropic pixels for processing...
    proc_xspacing_mm = proc_zspacing_mm;
    proc_yspacing_mm = proc_zspacing_mm;
    
    proc_xfov_meter = single(acq_fovx_mm) / 1000.0;
    proc_yfov_meter = single(acq_fovy_mm) / 1000.0;
    
    proc_nx = round(acq_fovx_mm / proc_xspacing_mm);
    proc_ny = round(acq_fovy_mm / proc_yspacing_mm);
    
    proc_nz = dicom_nz + 2;
    proc_zfov_meter = single(proc_nz) * single(proc_zspacing_mm) / 1000.0;

    % =========================================
    
    nx = int32(dicom_nx);  % DICOM xsize (usually, nx > inx)
    ny = int32(dicom_ny);  % DICOM ysize
    nz = int32(dicom_nz);  % number of acquired (and DICOM) slices
    nt = int32(acq_nt);    % number of time offsets in the motion data

    inx = int32(proc_nx);  % xsize at which to process data
    
    iny = int32(proc_ny);  % ysize at which to process data
    if (iny == 0)
        iny = inx;
    end
    
    if (inx == 0)
        % Dimension error...
        fprintf('\n    Dimension Error (inx == 0)!');
        fprintf('\n### Standard_MRE_Pipeline ABORTED ###\n');
        results = [];
        return
    end
    
%     % Calculate field of view, in METERS {x y z]...
%     fovx = single(motion_imgInfo.PixelSpacing(1)) * single(nx) / 1000.0;
%     fovy = single(motion_imgInfo.PixelSpacing(2)) * single(ny) / 1000.0;
%     fovz = single(motion_imgInfo.SpacingBetweenSlices) * single(nz + 2) / 1000.0;
%         % Error corrected 2015_0623: "nz + 2" was just "nz,"
%         % resulting in erroneous processing FOVz of 0.144 meter!

    % Define FOV in meters, for di_JDT2()...
    % fov = [fovx fovy fovz];   // NOTE axis order
    fov = [proc_xfov_meter proc_yfov_meter proc_zfov_meter];

    if (f == 0)
        % Get driver frequency from user...
        f = 60;   % Default driver frequency, Hz
        instr = input('\nEnter driver frequency, Hz (<Enter> for 60Hz) : ','s');
        if (~isempty(instr))
            fin = str2double(instr);
            if (fin > 0)
                f = fin;
            end
        end
    end
    
    resultsdir = fullfile(directory,'pipeline_results');
    resultsdir2 = resultsdir;
    if ~isempty(sub_directory)
        resultsdir2 = fullfile(resultsdir,sub_directory);
    end
    
    n_roi = 0;
    cancel = 0;
    roi_filename = [];
    if (~isempty(find(rlist == 41, 1)) || ~isempty(find(rlist == 51, 1)))
        % Have user select ROI (.nii) file here...
        % (Need to have valid ROI file, and need to
        %  know number of regions, to rebuild rlist)...
        
        starting_dir = resultsdir;
        try
            imp = strfind(resultsdir, '_motion/pipeline_results');
            if (isempty(imp))
                % This is not the 60Hz results directory;
                % Change starting directory to be the 60Hz results directory...
                [success,file_path,] = fileattrib([resultsdir '/../..']);
                find_dir = dir([resultsdir '/../../*_motion']);
                if (success && length(find_dir) == 1)
                    starting_dir = fullfile(file_path.Name, [find_dir(1).name '/pipeline_results']);
                end
            end
        catch
            starting_dir = resultsdir;
        end

        file = 0; %#ok<NASGU>
        
        % Instead of just displaying the ROI file names (in uigetfile window, below), 
        % list them by date last modified (most-recent first), and list the region 
        % numbers they each contain, so that the user can select more easily...
        
        % Former code...
        %       [file, path, filter_index] = uigetfile('*NF.nii', 'Select input ROI (.nii) file', starting_dir);   
        %       if (filter_index == 0)
        %           % User selected Cancel, or closed selection menu...
        %           results = [];
        %           cancel = 1;
        %           return
        %       end
        
        roi_files = [dir(fullfile(starting_dir, '*NF.nii')); ...
                     dir(fullfile(starting_dir, '*MS.nii')); ...
                     dir(fullfile(starting_dir, '*KMP.nii')); ...
                     dir(fullfile(starting_dir, '*KP.nii')); ...
                     dir(fullfile(starting_dir, '*ROI.nii'))];
                 
        if isempty(roi_files)
            % No ROI files found; make user find and select file...
            [file, path, filter_index] = uigetfile('*NF.nii; *MS.nii; *KMP.nii; *KP.nii; *ROI.nii', 'Select input ROI (.nii) file', starting_dir);   
            if (filter_index == 0)
                % User selected Cancel, or closed selection menu...
                results = [];
                cancel = 1;
                return
            end
        else
            path = starting_dir;
            nroi = numel(roi_files);
            if (nroi > 1)
                % Sort files by modification date; newest first...
                for i = 1:nroi - 1
                    for j = i + 1:nroi
                        if (roi_files(i).datenum < roi_files(j).datenum)
                            % Swap more recent file(j) with file(i)...
                            temp = roi_files(i);
                            roi_files(i) = roi_files(j);
                            roi_files(j) = temp;
                        end
                    end
                end
            end
            
            % Determine max region number contained in each ROI file...
            max_roi = zeros(nroi);
            for i = 1:nroi
                temp_filename = fullfile(path, roi_files(i).name);
                roi_array = spm_read_vols(spm_vol(temp_filename));
                max_roi(i) = max_value(roi_array);
            end
            
            if (0 < nroi && nroi < select_roi_n)
                % select_roi_n is invalid; cancel processing of this subject...
                disp(['Motion directory: ', directory]);
                fprintf('  (Specified ROI file [%d] does not exist)\n', select_roi_n);
                results = [];
                cancel = 1;
                return
            end
            
            % Display ROI file(s) for user selection...
            fprintf('  Select input ROI (.nii) file to use...\n');
            for i = 1:nroi
                if (select_roi_n > 0 && select_roi_n == i && nroi > 1)
                    % Display auto-selected ROI name (of two-or-more) in bold text...
                    fprintf('\n  [%d] ', i);
                    cprintf('*text', '%-60s', roi_files(i).name);
                    fprintf('  %s   (max region = %d)', roi_files(i).date, max_roi(i));
                else
                    fprintf('\n  [%d] %-60s  %s   (max region = %d)', i, roi_files(i).name, roi_files(i).date, max_roi(i));
                end
            end
            fprintf('\n\n');
            
            if (n_erode_roi == 0)
                fprintf('  [[ Net ROI erosion = %d ]]\n', n_erode_roi);
            else
                cprintf('*text', '  [[ Net ROI erosion = %d ]]\n', n_erode_roi);
            end
            
            if (select_roi_n > 0)
                selected_roi = select_roi_n;
                fprintf('  (Roi file [%d] auto-selected)\n', selected_roi);
            else
                selected_roi = -1;
            end
            
            while (selected_roi < 0 || selected_roi > nroi)
                selected_roi = input('  Enter number of ROI file to use (0 to select a different file): ');
            end
            
            fprintf('\n');
            if (selected_roi == 0)
                [file, path, filter_index] = uigetfile('*NF.nii; *MS.nii; *KMP.nii; *KP.nii; *ROI.nii', 'Select input ROI (.nii) file', starting_dir);   
                if (filter_index == 0)
                    % User selected Cancel, or closed selection menu...
                    results = [];
                    cancel = 1;
                    return
                end
            else
                file = roi_files(selected_roi).name;
            end
        end

        if (file ~= 0)
            roi_filename = fullfile(path, file);
            if (~exist(roi_filename, 'file') == 2)
                fprintf('\n    ROI file not found!');
                fprintf('\n### Standard_MRE_Pipeline ABORTED ###\n');
                results = [];
                return
            end

            % Read ROI array from file...
            roi_array = spm_read_vols(spm_vol(roi_filename));
            max_roi = max_value(roi_array);
            if (max_roi < 1)
                fprintf('\n    ROI array contains no regions!');
                fprintf('\n### Standard_MRE_Pipeline ABORTED ###\n');
                results = [];
                return
            end
            
            % Build list of roi numbers (not necessarily contiguous integers!)...
            roi_list = [];
            for i = 1:max_roi
                if ~isempty(find(roi_array == i, 1))
                    % Value 'i' exists (at least once) in roi_array...
                    % (Add 'i' to list of ROI region numbers)
                    roi_list = [roi_list i]; %#ok<AGROW>
                end
            end
            
            if (isempty(find(roi_list == 1, 1)) && max_roi > 5)
                % roi_array does not contain value '1', but DOES
                % contain sub-volume numbers, so add region '1'...
                roi_list = [1 roi_list];
            end
            
            n_roi = length(roi_list);
            
            sizes = size(roi_array);
            if sizes(3) == nz
                % Add two slices of zeros; one to each side of roi_array volume...
                roi_array = cat(3,zeros(iny,inx),roi_array,zeros(iny,inx));
            end
            
        end
        
        if (~isempty(find(rlist == 41, 1)))       
            % Rebuild rlist, to include all ROI regions...
            % n_roi = max_value(roi_array);
            new_rlist = zeros(1, length(rlist) - 1 + n_roi);
            j = 1;
            for i = 1:length(rlist)
                if (rlist(i) == 41)
                    for r = 1:n_roi
                        roi_num = 40 + roi_list(r);
                        new_rlist(j) = roi_num;
                        j = j + 1;
                    end
                else
                    new_rlist(j) = rlist(i);
                    j = j + 1;
                end
            end

            rlist = new_rlist;
        end
        
    end
    
    if (display_selected_dirs)
        disp( '---------------------');
        disp(['    Driver frequency: ', num2str(f), ' Hz']);
        try
            if (abs(proc_xspacing_mm - 3.0) > 0.01 || abs(proc_yspacing_mm - 3.0) > 0.01 || abs(proc_zspacing_mm - 3.0) > 0.01)
                fprintf('### Note pixel spacing for processing: %.3f, %.3f, %.3f mm ###\n', proc_xspacing_mm, proc_yspacing_mm, proc_zspacing_mm);
            end
        catch
        end
        disp(['    Motion directory: ', directory]);
        disp([' No-motion directory: ', directory_nm]);
        disp(['        T1 directory: ', t1_directory]);
        disp( '---------------------');
    end

    if ~exist('rlist','var')
        rlist = DEFAULT_rlist;
    end

    reg0_index = find(rlist == 0, 1);
    if (~isempty(reg0_index) && reg0_index > 1)
        % Move region '0' to BEGINNING of list...
        rlist = [0, rlist(1:reg0_index-1) rlist(reg0_index+1:end)];
    end
    
    if (pause_flag == 1)
        if (~isequal(rlist,DEFAULT_rlist) && ~isequal(rlist,[0,DEFAULT_rlist]))
            fprintf('\nSelected regions %2d: %s\n', rlist(1), region_name(rlist(1)));
            for rnum = 2:size(rlist, 2)
                fprintf('                 %2d: %s\n', rlist(rnum), region_name(rlist(rnum)));
            end
            disp( '---------------------');
            disp('Press any key to begin processing (Ctrl-C to abort)... ');
            pause;
            disp('Please wait... ');
        end
    end

    original_rlist = rlist;
    
    % Ensure that 'rlist' includes (begins with) region '0'...
    no_reg0 = 0;
    reg0_index = find(rlist == 0, 1);
    if isempty(reg0_index)
%         % Add region '0' to beginning of rlist...
%         rlist = [0, rlist];
        no_reg0 = 1;    % rlist did not originally include region '0'
    end
    
    if (special_csf_fraction ~= 0)
        fprintf('\n### NOTE: using special CSF fraction (%.3f) for brain_mask ###\n', special_csf_fraction);
    end
    
%     if (pause_flag == 1)
%         disp('Press any key to begin processing (Ctrl-C to abort)... ');
%         pause;
%     end

    % Disable interpreter in figures (waitbars)...
    % (This avoids special characters in waitbar text)
    set(0, 'DefaulttextInterpreter', 'none');

    elapsed = tic;
    currentDir=pwd;

    % create temporary processing directory and cd into it
    procdir =  get_procdir('MRE','brainPipeline');
    cd(procdir);

    inz = int32(proc_nz); % number of slices after adding a slice of zeros on the top and bottom of the volume
                          % (this is done so that erosion will account for edge artifacts in the top and bottom slices)

    patientid_nospaces = t1_imgInfo.PatientID; % get the patient ID
    patientid_nospaces(patientid_nospaces==' ')=[]; % remove spaces
    patientid_nospaces(patientid_nospaces=='-')=[]; % remove dashes
    filename_idToSeries = [patientid_nospaces '_' num2str(t1_imgInfo.StudyDate) num2str(t1_imgInfo.SeriesTime) '_' num2str(t1_imgInfo.SeriesNumber)]; % concatentate with study date, series time, and series number

    % calculate the total number of images
    megDirns=6;
    noImgs = nz*nt*megDirns*2;

    % load complex images
    if fopen(fullfile(directory,'cimgs.mat'))==-1 % if a mat file does NOT exist that contains all the images
        % load the dicom images
        h=waitbar(0,'loading motion images...','Name','Brain MRE Pipeline');
        % initialize a variable to hold the dicom images
        data=zeros(ny,nx,noImgs);
        % read in each image
        for ii=1:noImgs
            data(:,:,ii) = dicomread([directory 'i' sprintf('%03d',ii) '.dcm']);
            waitbar(single(ii)/single(noImgs),h)
        end
        delete(h);  % Close the waitbar window
        pause(0.1);

        % cast data to floating point
        data=single(data);
        % reshape data to put the real and imaginary parts in the 3rd dimension
        data=reshape(data,[ny,nx,2,nz*nt*megDirns]);
        % create complex images
        cimgs = squeeze(data(:,:,1,:)+1i*data(:,:,2,:));
        % reshape complex images so the slices are in the 3rd
        % dimension, time offsets in the 4th and motion encoding
        % direction in the 5th
        cimgs = reshape(cimgs,[ny,nx,megDirns, nz,nt]);
        cimgs = permute(cimgs,[1 2 4 5 3]);

        % get the physical position of the top and bottom slice
        info_slice1=dicominfo([directory 'i' sprintf('%03d',1) '.dcm']);
        info_sliceEnd=dicominfo([directory 'i' sprintf('%03d',noImgs) '.dcm']);

        % if the volume starts of the top of the head then flip the
        % slices
        if info_slice1.SliceLocation > info_sliceEnd.SliceLocation
            cimgs=flipdim(cimgs,3);
        end

        % move the data into nii orientation
        cimgs = permute(cimgs,[2 1 3 4 5]);
        cimgs = flipdim(cimgs,2);
        cimgs = flipdim(cimgs,1);

        % save the images to a mat file
        save([directory 'cimgs.mat'],'cimgs')

    else % if this was already done, just load the mat file
        load(fullfile(directory,'cimgs.mat'),'cimgs')
    end

    % calculate phase-difference complex images between the positive
    % and negative motion encoding images
    % the magnitude is the geometric mean of pos. and neg. encoded
    % images
    % the phase is the phase of the complex conjugate product
    pdcimgs(:,:,:,:,1) = sqrt(abs(cimgs(:,:,:,:,1)).*abs(cimgs(:,:,:,:,2))).*exp(1i*angle(cimgs(:,:,:,:,1).*conj(cimgs(:,:,:,:,2)))); % +x and -x
    pdcimgs(:,:,:,:,2) = sqrt(abs(cimgs(:,:,:,:,3)).*abs(cimgs(:,:,:,:,4))).*exp(1i*angle(cimgs(:,:,:,:,3).*conj(cimgs(:,:,:,:,4)))); % +y and -y
    pdcimgs(:,:,:,:,3) = sqrt(abs(cimgs(:,:,:,:,5)).*abs(cimgs(:,:,:,:,6))).*exp(1i*angle(cimgs(:,:,:,:,5).*conj(cimgs(:,:,:,:,6)))); % +z and -z

    % Interpolate slices to the desired X/Y resolution for PROCESSING...
    if inx~=nx || iny~=ny
        % take the 2D FFT
        ksp = ifftshift(ifftshift(ifft(ifft(pdcimgs,[],1),[],2),1),2);
        % grab the center of k-space
        pdcimgs_interp = ksp(((ny-iny)/2+1):((ny-iny)/2+iny),((nx-inx)/2+1):((nx-inx)/2+inx),:,:,:);
        % go back to image space
        pdcimgs_interp = fft(fft(fftshift(fftshift(pdcimgs_interp,1),2),[],1),[],2);
    else
        pdcimgs_interp=pdcimgs;
    end

    laDir = fullfile(filesep,'mnt','la','mc',patientid_nospaces(1:3),patientid_nospaces,t1_imgInfo.StudyDate,'Protocol_UNKNOWN');
    
    if (~exist(laDir, 'dir') && (~exist(resultsdir,'dir') || length(dir(fullfile(resultsdir,'r*.nii'))) < 1))
        disp(['NOT FOUND: ' laDir]);
        laDir = fullfile(filesep,'mnt','la','mc',patientid_nospaces(1:3),patientid_nospaces,t1_imgInfo.StudyDate,'');
        laDir = uigetdir(laDir, 'Select data directory containing "ATLAS" (Cancel to abort)...');
        if (laDir == 0)
            fprintf('\n### Standard_MRE_Pipeline ABORTED ###\n');
            fprintf('(%s)\n\n', directory);
            results = [];
            return
        end
    end
    
    process_data = 'yes';
    
    gm_fn  = ['c1M_' filename_idToSeries '_MT1_GW_N3m_HMRF.nii']; % name of the Gray Matter segmentation image
    wm_fn  = ['c2M_' filename_idToSeries '_MT1_GW_N3m_HMRF.nii']; % White Matter image
    csf_fn = ['c3M_' filename_idToSeries '_MT1_GW_N3m_HMRF.nii']; % CSF image
    wm22_fn  = ['M_' filename_idToSeries '_MT1_GW_N3m_STAND400_WM_LOBAR_22_ATLAS.nii'];
    aal35_fn = ['M_' filename_idToSeries '_MT1_GW_N3m_STAND400_AAL_35_ATLAS.nii'];
    tiv_fn   = ['M_' filename_idToSeries '_MT1_GW_N3_STAND400_TIV_TIVMASK_ERODED.nii'];
    
    proc_gm_fn    = fullfile(procdir,gm_fn);
    proc_wm_fn    = fullfile(procdir,wm_fn);
    proc_csf_fn   = fullfile(procdir,csf_fn);
    proc_wm22_fn  = fullfile(procdir,wm22_fn);
    proc_aal35_fn = fullfile(procdir,aal35_fn);
    proc_tiv_fn   = fullfile(procdir,tiv_fn);

    if ~exist(resultsdir,'dir') || length(dir(fullfile(resultsdir,'r*.nii'))) < 1

        atlasDir = fullfile(laDir,'ATLAS','STAND400',[num2str(t1_imgInfo.StudyDate) num2str(t1_imgInfo.SeriesTime) '_' num2str(t1_imgInfo.SeriesNumber) '_MT1_GW_N3m']);
        volsDir = fullfile(laDir,'vols');

        t1Filename = fullfile(atlasDir,['M_' filename_idToSeries '_MT1_GW_N3m.nii']);

        [~,fnroot,fnext]=fileparts(t1Filename);
        proc_t1_fn = fullfile(procdir,[fnroot fnext]);

        aal35Filename = fullfile(atlasDir,aal35_fn);
        aal120Filename = fullfile(atlasDir,['M_' filename_idToSeries '_MT1_GW_N3m_STAND400_AAL_120_ATLAS.nii']);
        wm22Filename = fullfile(atlasDir,wm22_fn);
        wm26Filename = fullfile(atlasDir,['M_' filename_idToSeries '_MT1_GW_N3m_STAND400_WM_LOBAR_26_ATLAS.nii']);

        [~,fnroot,fnext]=fileparts(aal120Filename);
        proc_aal120_fn = fullfile(procdir,[fnroot fnext]);

        [~,fnroot,fnext]=fileparts(wm26Filename);
        proc_wm26_fn = fullfile(procdir,[fnroot fnext]);

        gmFilename = fullfile(atlasDir,gm_fn); % name of the resliced GM segmentation image
        wmFilename = fullfile(atlasDir,wm_fn); % WM image
        csfFilename = fullfile(atlasDir,csf_fn); % CSF image

        bmFilename = fullfile(volsDir,['M_' filename_idToSeries '_MT1_GW_N3_brainmask.nii']);
        ventFilename = fullfile(volsDir,['M_' filename_idToSeries '_MT1_GW_N3_ventmask2.nii']);
        tivFilename = fullfile(volsDir,tiv_fn);

        [~,fnroot,fnext]=fileparts(bmFilename);
        proc_bm_fn = fullfile(procdir,[fnroot fnext]);

        [~,fnroot,fnext]=fileparts(ventFilename);
        proc_vent_fn = fullfile(procdir,[fnroot fnext]);

        disp(' ');
        if ~exist(resultsdir,'dir')
            disp('(Preprocessed images not found)')
        elseif length(dir(fullfile(resultsdir,'r*.nii'))) < 1
            disp('(Resliced images not found)')
        end
        disp('Performing registration and reslicing steps...')
        disp(' ');

        if ~exist(atlasDir, 'dir')
            cprintf('*text', '    Grinder Atlas directory does not (yet) exist for this exam.');
            fprintf('\n\n### Standard_MRE_Pipeline ABORTED ###\n\n');
            results = [];
            return
        end
        
        % Copy working files to (temporary) processing directory...
        try
            copyfile(t1Filename,procdir);
            copyfile(aal35Filename,procdir);
            copyfile(aal120Filename,procdir);
            copyfile(wm22Filename,procdir);
            copyfile(wm26Filename,procdir);
            copyfile(gmFilename,procdir);
            copyfile(wmFilename,procdir);
            copyfile(csfFilename,procdir);
            copyfile(bmFilename,procdir);
            copyfile(ventFilename,procdir);
            copyfile(tivFilename,procdir);
        catch
            cprintf('*text', '    Required file(s) not found!');
            if ~exist(t1Filename, 'file')
                fprintf('\n    (Did Grinder select different series for T1?)');
            end
            fprintf('\n\n### Standard_MRE_Pipeline ABORTED ###\n\n');
            results = [];
            return
        end

        magn=abs(pdcimgs_interp);
        magn=mean(magn,5);
        magn=squeeze(magn(:,:,:,1));
        magnFn = fullfile(procdir,'magn.nii');

      % Affine transformation matrix (scaling and translation)...
      
      % V.mat = [ 3 0 0 -120;       These (former) hard-coded values were for 3-mm
      %           0 3 0 -120;       pixels (y, x, and z), translated to the center
      %           0 0 3  -72;       of the 240-mm(y,x) and 144-mm(z) 'magn' volume.
      %           0 0 0    1 ];     (That is, 240/2 = 120mm, and 144/2 = 72mm)
      
        half_xfov_mm = proc_xfov_meter * 500.0;
        half_yfov_mm = proc_yfov_meter * 500.0;
        half_zfov_mm = proc_zfov_meter * 500.0;
        
        V.mat = [ proc_yspacing_mm                 0                 0  -half_yfov_mm ; ...
                                 0  proc_xspacing_mm                 0  -half_xfov_mm ; ...
                                 0                 0  proc_zspacing_mm  -half_zfov_mm ; ... 
                                 0                 0                 0              1 ];
        V.fname = magnFn;
        V.dim = [iny inx nz];
        V.dt = [16 0];
        V.pinfo = [1; 0; 352];
        V.n = [1 1];
        V.descrip = '';
        V.private = [];

        spm_write_vol(V,magn);

        % register and reslice images to T1-FLAIR image
        other = {proc_aal35_fn,proc_aal120_fn,proc_wm22_fn,proc_wm26_fn,proc_gm_fn,proc_wm_fn,proc_csf_fn,proc_bm_fn,proc_vent_fn,proc_tiv_fn};

        mls_atlas_coreg(magnFn,proc_t1_fn,other,6,1);

        others_linInterp = char(magnFn,proc_t1_fn,proc_gm_fn,proc_wm_fn,proc_csf_fn);

        flags.mask   = 0;
        flags.mean   = 0;
        flags.interp = 1;
        flags.which  = 1;

        spm_reslice(others_linInterp,flags);

        others_nearNeigh = char(magnFn,proc_aal35_fn,proc_aal120_fn,proc_wm22_fn,proc_wm26_fn,proc_bm_fn,proc_vent_fn,proc_tiv_fn);
        flags.interp = 0;
        spm_reslice(others_nearNeigh,flags);

        user_wait = tic;
        
        % User MUST visually verify that registration is reasonable...
        system('mricron magn.nii -o rc1*.nii -o rc2*.nii -o rc3*.nii &');
        fprintf('\nVERIFY registration visually (MRIcron window)\n');
        fprintf('(You may then close the MRIcron window)\n\n');
        process_data = 'ask';
        if (reg_reslice_only)
            fprintf('Press any key to continue...\n');
            pause
        else
            while (~strcmpi(process_data,'yes') && ~strcmpi(process_data,'no'))
                process_data = input('Process regions of this data (yes/no)? : ','s');
            end
        end

        user_wait_t = toc(user_wait);
        
        if (reg_reslice_only)
            if ~exist(resultsdir,'dir')
                mkdir(resultsdir)
            end

            if length(dir(fullfile(resultsdir,'r*.nii'))) < 1
                copyfile(fullfile(procdir,'r*.nii'),resultsdir);
                copyfile(fullfile(procdir,'magn.nii'),resultsdir);
            end
        end
        
    else
        
        user_wait_t = 0;
        if (reg_reslice_only)
            fprintf('\nPreprocessed images already exist.\n');
            fprintf('(Will not re-register or reslice data)\n\n');
        else
            fprintf('\nPreprocessed (registered) images found;\n');
            fprintf('Copying to processing directory...\n');
            copyfile(fullfile(resultsdir,'r*.nii'),procdir);
            copyfile(fullfile(resultsdir,'magn.nii'),procdir);
        end
        
    end

    if (reg_reslice_only)
        mls_cleanup_dir(procdir);
        cd(currentDir);
        results = [];
        return
    end

    if (~strcmpi(process_data, 'yes'))
        mls_cleanup_dir(procdir);
        cd(currentDir);
        results = [];
        t = toc(elapsed);
        fprintf('Elapsed time %.1f minutes\n', t / 60);
        return
    end
        
    gm=spm_read_vols(spm_vol(prepend(proc_gm_fn,'r')));
    gm=cat(3,zeros(iny,inx),gm,zeros(iny,inx));   % load gm image and add slices of zeros to each side of the volume

    wm=spm_read_vols(spm_vol(prepend(proc_wm_fn,'r')));
    wm=cat(3,zeros(iny,inx),wm,zeros(iny,inx));   % same for wm

    csf=spm_read_vols(spm_vol(prepend(proc_csf_fn,'r')));
    csf=cat(3,zeros(iny,inx),csf,zeros(iny,inx)); % and csf

    if (special_csf_fraction <= 0)
        % Normal CSF threshold (about 50%)...
        brain_mask = csf < (gm+wm); % create the brain mask where gm+wm content is greater than csf content
    else
        % Special/custom CSF threshold fraction for testing...
        brain_mask = csf < (special_csf_fraction * (gm + wm + csf));
    end
    
    % load the 22 region lobar WM atlas
    atlas=spm_read_vols(spm_vol(prepend(proc_wm22_fn,'r')));
    atlas=cat(3,zeros(iny,inx),atlas,zeros(iny,inx)); % add slices of zeros

    % load 35 region AAL atlas
    aal35=spm_read_vols(spm_vol(prepend(proc_aal35_fn,'r')));
    aal35=cat(3,zeros(iny,inx),aal35,zeros(iny,inx)); % add slices of zeros

    % Save both atlas volumes (converted to unsigned, 8-bit integers) to files in motion directory...
    array_name = ['atlas_' patientid_nospaces '_scan' date9(t1_imgInfo.StudyDate)];
    eval([array_name ' = uint8(atlas);']);
    matFile = fullfile(directory, [array_name '.mat']);
    save(matFile, array_name);
    clear array_name;
    
    array_name = ['aal35_' patientid_nospaces '_scan' date9(t1_imgInfo.StudyDate)];
    eval([array_name ' = uint8(aal35);']);
    matFile = fullfile(directory, [array_name '.mat']);
    save(matFile, array_name);
    clear array_name;
    
    % Apply striping filter...
    vtk = 3; % number of samples to keep at the center of k-space
    vfc = (vtk-1)/2; % number of samples from dc
    lpfilter = zeros(iny,inx); % initialize filter as all zeros
    lpfilter((iny/2+1-vfc):(iny/2+1+vfc),(inx/2+1-vfc):(inx/2+1+vfc)) = 1; % set the central voxels to 1

    % take the 2D FFT
    ksp = ifftshift(ifftshift(ifft(ifft(pdcimgs_interp,[],1),[],2),1),2);
    % multiple k-space data by the lowpass filter made above and go
    % back to image space
    lpfiltered = fft(fft(fftshift(fftshift(ksp.*repmat(lpfilter,[1 1 nz nt 3]),1),2),[],1),[],2);
    % keep the magnitude from the original pdcimgs and set the
    % phase equal to the difference between the original phase of
    % pdcimgs_interp and the lowpass filtered images
    pdcimgs_interp = abs(pdcimgs_interp).*exp(1i*angle(pdcimgs_interp.*conj(lpfiltered)));

    % add a slice of zeros to the top and bottom of the volume
    pdcimgs_interp = cat(3,zeros(iny,inx,1,nt,3),pdcimgs_interp,zeros(iny,inx,1,nt,3));

    % calculate the magnitude image
    magn = mean(abs(pdcimgs_interp),5);

    % Save magnitude image to file in motion directory...
    array_name = ['disp_mag_' patientid_nospaces '_scan' date9(t1_imgInfo.StudyDate)];
    eval([array_name ' = single(magn);']);
    matFile = fullfile(directory, [array_name '.mat']);
    save(matFile, array_name);
    clear array_name;
    
    % calculate displacements in each the x, y and z directions
    ux = angle(pdcimgs_interp(:,:,:,:,1));
    uy = angle(pdcimgs_interp(:,:,:,:,2));
    uz = angle(pdcimgs_interp(:,:,:,:,3));

%### TESTING: Save displacement images to .mat files...
% try
%     magfilename = ['brain_' mc10(patientid_nospaces) '_' date9(motion_imgInfo.StudyDate) '_Mag.mat'];
%     xfilename   = ['brain_' mc10(patientid_nospaces) '_' date9(motion_imgInfo.StudyDate) '_Xdisp.mat'];
%     yfilename   = ['brain_' mc10(patientid_nospaces) '_' date9(motion_imgInfo.StudyDate) '_Ydisp.mat'];
%     zfilename   = ['brain_' mc10(patientid_nospaces) '_' date9(motion_imgInfo.StudyDate) '_Zdisp.mat'];
%     save(magfilename, 'magn');
%     save(xfilename, 'ux');
%     save(yfilename, 'uy');
%     save(zfilename, 'uz');
%     copyfile(fullfile(procdir, magfilename),resultsdir);
%     copyfile(fullfile(procdir, xfilename),resultsdir);
%     copyfile(fullfile(procdir, yfilename),resultsdir);
%     copyfile(fullfile(procdir, zfilename),resultsdir);
%     keyboard
% catch
%     disp('TESTING error');
%     return
% end


%### TESTING: Save displacement images to binary (.r4) files...
% try
%     magfilename = ['brain_' mc10(patientid_nospaces) '_' date9(motion_imgInfo.StudyDate) '_Mag.r4'];
%     xfilename   = ['brain_' mc10(patientid_nospaces) '_' date9(motion_imgInfo.StudyDate) '_Xdisp.r4'];
%     yfilename   = ['brain_' mc10(patientid_nospaces) '_' date9(motion_imgInfo.StudyDate) '_Ydisp.r4'];
%     zfilename   = ['brain_' mc10(patientid_nospaces) '_' date9(motion_imgInfo.StudyDate) '_Zdisp.r4'];
%     disp(magfilename);
%     disp(xfilename);
%     disp(yfilename);
%     disp(zfilename);
%     fwrite_single_binary_file(magn, magfilename);
%     fwrite_single_binary_file(ux, xfilename);
%     fwrite_single_binary_file(uy, yfilename);
%     fwrite_single_binary_file(uz, zfilename);
%     copyfile(fullfile(procdir, magfilename),resultsdir);
%     copyfile(fullfile(procdir, xfilename),resultsdir);
%     copyfile(fullfile(procdir, yfilename),resultsdir);
%     copyfile(fullfile(procdir, zfilename),resultsdir);
%     keyboard
% catch
%     disp('TESTING error');
%     return
% end
% keyboard

% Read and Save Total Intracranial Volume mask...
% tiv_mask = spm_read_vols(spm_vol(prepend(proc_tiv_fn,'r')));
% tiv_mask = cat(3,zeros(iny,inx),tiv_mask,zeros(iny,inx)); % (add slice of zeros to each side of the volume)
% fwrite_binary_file(tiv_mask,'tiv_mask.n1','char');            
% copyfile(fullfile(procdir,'tiv_mask.n1'),resultsdir);
% keyboard
%### (TESTING)

    % find some voxels that cover most of the brain to use as an ROI to
    % measure the divergence

    % create kernels for calculating derivatives
    dx = [-1 0 1];
    dy = permute(dx,[2 1 3]);
    dz = permute(dx,[3 1 2]);

    % determine coordinate system by minimizing the divergence
    ampOfDiv=zeros(1,8); % initialize variable that will hold the amplitude of the divergence for each of the 8 possible coordinate choices
    for ii=1:8 % for each of the candidate coordinate choices
        if ii==1     % +x, +y, +z
            px = exp(1i *  ux);
            py = exp(1i *  uy);
            pz = exp(1i *  uz);
        elseif ii==2 % -x, +y, +z
            px = exp(1i * -ux);
            py = exp(1i *  uy);
            pz = exp(1i *  uz);
        elseif ii==3 % +x, -y, +z
            px = exp(1i *  ux);
            py = exp(1i * -uy);
            pz = exp(1i *  uz);
        elseif ii==4 % +x, +y, -z
            px = exp(1i *  ux);
            py = exp(1i *  uy);
            pz = exp(1i * -uz);
        elseif ii==5 % +y, +x, +z
            px = exp(1i *  uy);
            py = exp(1i *  ux);
            pz = exp(1i *  uz);
        elseif ii==6 % -y, +x, +z
            px = exp(1i * -uy);
            py = exp(1i *  ux);
            pz = exp(1i *  uz);
        elseif ii==7 % +y, -x, +z
            px = exp(1i *  uy);
            py = exp(1i * -ux);
            pz = exp(1i *  uz);
        elseif ii==8 % +y, +x, -z
            px = exp(1i *  uy);
            py = exp(1i *  ux);
            pz = exp(1i * -uz);
        end

        % calculate the divergence
        divergence = real(-1i*conj(px).*(convn(px,dx,'same'))) + real(-1i*conj(py).*(convn(py,dy,'same'))) + real(-1i*conj(pz).*(convn(pz,dz,'same')));

        % take the first temporal harmonic
        divergence = 2*ifft(divergence,[],4);
        divergence = abs(divergence(:,:,:,2));

        % measure the median within the brain
        ampOfDiv(ii) = median(divergence(magn(:,:,:,1)>1000));

    end

    % determine which coordinate choice minimized the amplitude of the
    % divergence
    [~,coordinateChoice] = min(ampOfDiv);

    % reorganize ux, uy and uz accordingly
    if coordinateChoice==2
        ux=-ux;
    elseif coordinateChoice==3
        uy=-uy;
    elseif coordinateChoice==4
        uz=-uz;
    elseif coordinateChoice==5
        tmp=ux;
        ux=uy;
        uy=tmp;
        clear tmp
    elseif coordinateChoice==6
        uy=-uy;
        tmp=ux;
        ux=uy;
        uy=tmp;
        clear tmp
    elseif coordinateChoice==7
        ux=-ux;
        tmp=ux;
        ux=uy;
        uy=tmp;
        clear tmp
    elseif coordinateChoice==8
        uz=-uz;
        tmp=ux;
        ux=uy;
        uy=tmp;
        clear tmp
    end
    
    % create complex images with unit magnitude for use with SPUA
    px_motion = exp(1i*ux);
    py_motion = exp(1i*uy);
    pz_motion = exp(1i*uz);

    % apply adaptive curl calculation and filtering
    
    % initialize derivative kernels for edge voxels
    dx_left = [-1 1 0];
    dx_right = [0 -1 1];
    dy_top = [-1;1;0];
    dy_bottom = [0;-1;1];
    dz_top(1,1,1)=-1;
    dz_top(1,1,2)=1;
    dz_top(1,1,3)=0;
    dz_bottom(1,1,1)=0;
    dz_bottom(1,1,2)=-1;
    dz_bottom(1,1,3)=1;

    counter=1;
    % calculate stiffness for each region
    
    corrected_med_vs2_stiffness = zeros(length(rlist),1);
    med_vs2_stiffness           = zeros(length(rlist),1);
    mean_cmu_mag                = zeros(length(rlist),1);
    std_cmu_mag                 = zeros(length(rlist),1);
    med_cmu_mag                 = zeros(length(rlist),1);
    cmu_mag_iterative_cmedian   = zeros(length(rlist),1);
    cmu_mag_from_RImeans        = zeros(length(rlist),1);
    med_cmu_real                = zeros(length(rlist),1);
    med_cmu_imag                = zeros(length(rlist),1);
    cmu_mag_from_RImeds         = zeros(length(rlist),1);
    vs2_stiffness_from_RImeds = zeros(length(rlist),1);
    
    regionNames  = cell(length(rlist),1);
    noVoxels     = zeros(length(rlist),1);
    numberIters  = zeros(length(rlist),1);
    regionalSnrx = zeros(length(rlist),1);
    regionalSnry = zeros(length(rlist),1);
    regionalSnrz = zeros(length(rlist),1);
    regionalSnrModeX = zeros(length(rlist),1);
    regionalSnrModeY = zeros(length(rlist),1);
    regionalSnrModeZ = zeros(length(rlist),1);
    xMedianOfMostLikelySnrs = zeros(length(rlist),1);
    yMedianOfMostLikelySnrs = zeros(length(rlist),1);
    zMedianOfMostLikelySnrs = zeros(length(rlist),1);

    region_count = -0.5;
    region_waitbar = waitbar(0, ' ', 'Name',' Brain MRE Pipeline');
    newpos = get(region_waitbar,'Position');
    newpos(1) = 200;
    set(region_waitbar,'Position',newpos);

    disp(' ')

    % Find longest region name string...
    longest_rname = 0;
    for rr = rlist
        if (longest_rname < length(region_name(rr)))
            longest_rname = length(region_name(rr));
        end
    end
    
    % create a 6-connected neighborhood structural element (for erosions)...
    str_elem = zeros(3,3,3);
    str_elem(2,2,:) = 1;
    str_elem(:,2,2) = 1;
    str_elem(2,:,2) = 1;
    
    for rr = rlist

%         if (rr == rlist(1))
%             if (n_erode_roi == 0)
%                 cprintf('*text', '''Eroding'' all regions by ZERO layers...\n\n');
%             elseif (n_erode_roi < 0)
%                 cprintf('*text', 'DILATING all regions ');
%                 fprintf('(except no_mask) ');
%                 cprintf('*text', 'by %d layers...\n\n', abs(n_erode_roi));
%             elseif (n_erode_roi > 1)
%                 cprintf('*text', 'Eroding all regions ');
%                 fprintf('(except no_mask) ');
%                 cprintf('*text', 'by %d layers...\n\n', n_erode_roi);
%             end
%         end
        
        region_count = region_count + 1;
        waitbar(region_count/size(rlist, 2), region_waitbar, ['Processing region ' num2str(rr) ': ' region_name(rr) ' ...']);
        fprintf('Processing region %2d: %-*s', rr, longest_rname + 6, [region_name(rr) ' ... ']);
        region_name_str = region_name(rr);
        
        %%orig_region_voxels = 0;
        if (strcmpi(region_name_str,'no_mask'))
            % Read Total Intracranial Volume mask, (tiv_mask)...
            tiv_mask = spm_read_vols(spm_vol(prepend(proc_tiv_fn,'r')));
            tiv_mask = cat(3,zeros(iny,inx),tiv_mask,zeros(iny,inx));
                       % (add slices of zeros to each side of the volume)
            mask = ones(size(brain_mask));
            [region, rname] = define_region(rr, mask, atlas, aal35);
            regionNames{counter} = rname;
        elseif (strncmpi(region_name_str,'roi',3))
            % Build specific ROI region (roi_num)...
            roi_num = rr - 40;
            if (roi_num == 1)
                % ROI region 1 ("whole tumor") is defined as the union 
                % of region 1 (whole tumor outline/shell) and all tumor 
                % sub-volumes, regions 6, 7 (, 8, ...) ...
                region = (roi_array == roi_num) | (roi_array > 5);
            else
                region = (roi_array == roi_num);
            end
 
            rname = region_name(rr);
            %%orig_region_voxels = length(find(region));
            
            if (n_erode_roi < 1)
                % "Normal" ROI processing (no NET erosion)...
                % Dilate 'region' here, to cancel roi erosion later in processing...
                region = imdilate(region, str_elem);  % 6-connected-XYZ
                regionNames{counter} = rname;
            else %if (n_erode_roi == 1)
                % Leave 'region' as-is; net result (later) is a one-layer erosion of ROI...
                regionNames{counter} = [rname ' (E1)'];
%           else
%               % Erode ROI by (n_erode_roi - 1) layers; net result (later) is n_erode_roi layers...
%               for en = 1:(n_erode_roi - 1)
%                   region = imerode(region, str_elem);  % 6-connected-XYZ
%               end
%               regionNames{counter} = [rname ' (E' num2str(n_erode_roi) ')'];
            end
            
            mask = brain_mask | region;
        elseif (strncmpi(region_name_str,'cerebrum_MINUS_tumor',20))
            % Build special "cerebrum-minus-tumor" region...
            mask = brain_mask;
            cerebrum_region = define_region(1, mask, atlas, aal35);    % cerebrum
            
            % "whole tumor" is defined as the union of 
            % region 1 (whole tumor outline/shell) and all 
            % tumor sub-volumes, regions 6, 7 (, 8, ...) ...
            whole_tumor = (roi_array == 1) | (roi_array > 5);

            region = cerebrum_region - whole_tumor;
            rname = region_name_str;
            regionNames{counter} = rname;
        else
            mask = brain_mask;
            [region, rname] = define_region(rr, mask, atlas, aal35);
            regionNames{counter} = rname;
        end
       
        nt = int32(acq_nt);    % number of time offsets in the motion data
        px=px_motion;
        py=py_motion;
        pz=pz_motion; % use px, py and pz from the motion data

        maskt = repmat((mask & region),[1 1 1 nt]);

        % calculate duz/dy
        duzdy = zeros(iny,inx,inz,nt);
        curlMask_1 = zeros(iny,inx,inz,nt);
        
        ind = find(imerode(maskt,ones(3,1,1)));
        duzdy_temp = real(-1i*conj(pz).*convn(pz,dy,'same'));
        duzdy(ind) = duzdy_temp(ind);
        curlMask_1(ind) = 1;
        clear ind
        
        ind = find(~imerode(maskt,ones(3,1,1)) & imerode(maskt,abs(dy_top)));
        duzdy_temp = 2*real(-1i*conj(pz).*convn(pz,dy_bottom,'same'));
        duzdy(ind) = duzdy_temp(ind);
        curlMask_1(ind) = 1;
        clear ind
        
        ind = find(~imerode(maskt,ones(3,1,1)) & imerode(maskt,abs(dy_bottom)));
        duzdy_temp = 2*real(-1i*conj(pz).*convn(pz,dy_top,'same'));
        duzdy(ind) = duzdy_temp(ind);
        curlMask_1(ind) = 1;
        clear ind
        
        clear duzdy_temp

        % calculate duy/dz
        duydz = zeros(iny,inx,inz,nt);
        curlMask_2 = zeros(iny,inx,inz,nt);
        ind = find(imerode(maskt,ones(1,1,3))); duydz_temp = real(-1i*conj(py).*convn(py,dz,'same')); duydz(ind) = duydz_temp(ind); curlMask_2(ind) = 1; clear ind
        ind = find(~imerode(maskt,ones(1,1,3)) & imerode(maskt,abs(dz_top))); duydz_temp = 2*real(-1i*conj(py).*convn(py,dz_bottom,'same')); duydz(ind) = duydz_temp(ind); curlMask_2(ind) = 1; clear ind
        ind = find(~imerode(maskt,ones(1,1,3)) & imerode(maskt,abs(dz_bottom))); duydz_temp = 2*real(-1i*conj(py).*convn(py,dz_top,'same')); duydz(ind) = duydz_temp(ind); curlMask_2(ind) = 1; clear ind
        clear duydz_temp

        % calculate dux/dz
        duxdz = zeros(iny,inx,inz,nt);
        curlMask_3 = zeros(iny,inx,inz,nt);
        ind = find(imerode(maskt,ones(1,1,3))); duxdz_temp = real(-1i*conj(px).*convn(px,dz,'same')); duxdz(ind) = duxdz_temp(ind); curlMask_3(ind) = 1; clear ind
        ind = find(~imerode(maskt,ones(1,1,3)) & imerode(maskt,abs(dz_top))); duxdz_temp = 2*real(-1i*conj(px).*convn(px,dz_bottom,'same')); duxdz(ind) = duxdz_temp(ind); curlMask_3(ind) = 1; clear ind
        ind = find(~imerode(maskt,ones(1,1,3)) & imerode(maskt,abs(dz_bottom))); duxdz_temp = 2*real(-1i*conj(px).*convn(px,dz_top,'same')); duxdz(ind) = duxdz_temp(ind); curlMask_3(ind) = 1; clear ind
        clear duxdz_temp

        % calculate duz/dx
        duzdx = zeros(iny,inx,inz,nt);
        curlMask_4 = zeros(iny,inx,inz,nt);
        ind = find(imerode(maskt,ones(1,3,1))); duzdx_temp = real(-1i*conj(pz).*convn(pz,dx,'same')); duzdx(ind) = duzdx_temp(ind); curlMask_4(ind) = 1; clear ind
        ind = find(~imerode(maskt,ones(1,3,1)) & imerode(maskt,abs(dx_left))); duzdx_temp = 2*real(-1i*conj(pz).*convn(pz,dx_right,'same')); duzdx(ind) = duzdx_temp(ind); curlMask_4(ind) = 1; clear ind
        ind = find(~imerode(maskt,ones(1,3,1)) & imerode(maskt,abs(dx_right))); duzdx_temp = 2*real(-1i*conj(pz).*convn(pz,dx_left,'same')); duzdx(ind) = duzdx_temp(ind); curlMask_4(ind) = 1; clear ind
        clear duzdx_temp

        %calculate duy/dx
        duydx = zeros(iny,inx,inz,nt);
        curlMask_5 = zeros(iny,inx,inz,nt);
        ind = find(imerode(maskt,ones(1,3,1))); duydx_temp = real(-1i*conj(py).*convn(py,dx,'same')); duydx(ind) = duydx_temp(ind); curlMask_5(ind) = 1; clear ind
        ind = find(~imerode(maskt,ones(1,3,1)) & imerode(maskt,abs(dx_left))); duydx_temp = 2*real(-1i*conj(py).*convn(py,dx_right,'same')); duydx(ind) = duydx_temp(ind); curlMask_5(ind) = 1; clear ind
        ind = find(~imerode(maskt,ones(1,3,1)) & imerode(maskt,abs(dx_right))); duydx_temp = 2*real(-1i*conj(py).*convn(py,dx_left,'same')); duydx(ind) = duydx_temp(ind); curlMask_5(ind) = 1; clear ind
        clear duydx_temp

        % calculate dux/dy
        duxdy = zeros(iny,inx,inz,nt);
        curlMask_6 = zeros(iny,inx,inz,nt);
        ind = find(imerode(maskt,ones(3,1,1))); duxdy_temp = real(-1i*conj(px).*convn(px,dy,'same')); duxdy(ind) = duxdy_temp(ind); curlMask_6(ind) = 1; clear ind
        ind = find(~imerode(maskt,ones(3,1,1)) & imerode(maskt,abs(dy_top))); duxdy_temp = 2*real(-1i*conj(px).*convn(px,dy_bottom,'same')); duxdy(ind) = duxdy_temp(ind); curlMask_6(ind) = 1; clear ind
        ind = find(~imerode(maskt,ones(3,1,1)) & imerode(maskt,abs(dy_bottom))); duxdy_temp = 2*real(-1i*conj(px).*convn(px,dy_top,'same')); duxdy(ind) = duxdy_temp(ind); curlMask_6(ind) = 1; clear ind
        clear duxdy_temp

        curlMask = curlMask_1 & curlMask_2 & curlMask_3 & curlMask_4 & curlMask_5 & curlMask_6;
        clear curlMask_1 curlMask_2 curlMask_3 curlMask_4 curlMask_5 curlMask_6
        curlMask = squeeze(curlMask(:,:,:,1));

        curlx = duzdy - duydz;
        curly = duxdz - duzdx;
        curlz = duydx - duxdy;

        if (strcmpi(region_name(rr),'no_mask'))
            % Save curl data (for entire volume) to .mat file in motion directory...
            save([directory 'curl_data_no_mask.mat'], 'curlx', 'curly', 'curlz');
            % Save curl data (for entire TIV) to .mat file in motion directory...
            for i = 1:size(curlx, 4)
                curlx(:,:,:,i) = tiv_mask .* curlx(:,:,:,i); 
                curly(:,:,:,i) = tiv_mask .* curly(:,:,:,i); 
                curlz(:,:,:,i) = tiv_mask .* curlz(:,:,:,i); 
            end
            save([directory 'curl_data_tiv.mat'], 'curlx', 'curly', 'curlz');
            curlx = duzdy - duydz;
            curly = duxdz - duzdx;
            curlz = duydx - duxdy;
        elseif (strcmpi(region_name(rr),'cerebrum'))
            % Save curl data (for whole cerebrum) to .mat file in motion directory...
            save([directory 'curl_data_cerebrum.mat'], 'curlx', 'curly', 'curlz');
        elseif (strcmpi(region_name(rr),'brain'))
            % Save curl data (cerebrum plus cerebellum) to .mat file in motion directory...
            save([directory 'curl_data_brain.mat'], 'curlx', 'curly', 'curlz');
        end
        
        clear duzdy duydz duxdz duzdx duydx duxdy

        % process no-motion data

        clear pdcimgs pdcimgs_interp cimgs cimgs_unit_mag data
        nt = 2;

        % calculate first temporal harmonic
        curlxf = ifft(curlx,[],4)*2;
        curlyf = ifft(curly,[],4)*2;
        curlzf = ifft(curlz,[],4)*2;

        curlxf = squeeze(curlxf(:,:,:,2)).*(mask & region);
        curlyf = squeeze(curlyf(:,:,:,2)).*(mask & region);
        curlzf = squeeze(curlzf(:,:,:,2)).*(mask & region);

        if counter==1
            NOmotion_imgInfo = dicominfo(fullfile(directory_nm, 'i001.dcm'));
            nx = int32(NOmotion_imgInfo.Columns);
            ny = int32(NOmotion_imgInfo.Rows);
            noImgs = nz*nt*megDirns*2;
            if fopen([directory_nm 'cimgs.mat'])==-1
                h=waitbar(0,'loading no-motion images...', 'Name', 'Brain MRE Pipeline');
                data=zeros(ny,nx,noImgs);
                for ii=1:noImgs
                    data(:,:,ii) = dicomread([directory_nm 'i' sprintf('%03d',ii) '.dcm']);
                    waitbar(single(ii)/single(noImgs),h)
                end
                delete(h);  % Close the waitbar window
                pause(0.1);

                data=single(data);
                data=reshape(data,[ny,nx,2,nz*nt*megDirns]);
                cimgs = squeeze(data(:,:,1,:)+1i*data(:,:,2,:));
                cimgs = reshape(cimgs,[ny,nx,megDirns,nz,nt]);
                cimgs = permute(cimgs,[1 2 4 5 3]);

                info_slice1=dicominfo([directory_nm 'i' sprintf('%03d',1) '.dcm']);
                info_sliceEnd=dicominfo([directory_nm 'i' sprintf('%03d',noImgs) '.dcm']);

                if info_slice1.SliceLocation>info_sliceEnd.SliceLocation
                    cimgs=flipdim(cimgs,3);
                end

                cimgs = permute(cimgs,[2 1 3 4 5]);
                cimgs = flipdim(cimgs,2);
                cimgs = flipdim(cimgs,1);

                save([directory_nm 'cimgs.mat'],'cimgs')
            else
                load([directory_nm 'cimgs.mat'],'cimgs')
            end

            pdcimgs(:,:,:,:,1) = sqrt(abs(cimgs(:,:,:,:,1)).*abs(cimgs(:,:,:,:,2))).*exp(1i*angle(cimgs(:,:,:,:,1).*conj(cimgs(:,:,:,:,2))));
            pdcimgs(:,:,:,:,2) = sqrt(abs(cimgs(:,:,:,:,3)).*abs(cimgs(:,:,:,:,4))).*exp(1i*angle(cimgs(:,:,:,:,3).*conj(cimgs(:,:,:,:,4))));
            pdcimgs(:,:,:,:,3) = sqrt(abs(cimgs(:,:,:,:,5)).*abs(cimgs(:,:,:,:,6))).*exp(1i*angle(cimgs(:,:,:,:,5).*conj(cimgs(:,:,:,:,6))));

            if inx~=nx || iny~=ny
                ksp = ifftshift(ifftshift(ifft(ifft(pdcimgs,[],1),[],2),1),2);
                try
                    pdcimgs_interp = ksp(((ny-iny)/2+1):((ny-iny)/2+iny),((nx-inx)/2+1):((nx-inx)/2+inx),:,:,:);
                catch
                    fprintf('### Motion/No-motion array dimension mismatch ###\n\n');
                    mls_cleanup_dir(procdir);
                    cd(currentDir);
                    results = [];
                    return
                end
                pdcimgs_interp = fft(fft(fftshift(fftshift(pdcimgs_interp,1),2),[],1),[],2);
            else
                pdcimgs_interp=pdcimgs;
            end

            ksp = ifftshift(ifftshift(ifft(ifft(pdcimgs_interp,[],1),[],2),1),2);
            lpfiltered = fft(fft(fftshift(fftshift(ksp.*repmat(lpfilter,[1 1 nz nt 3]),1),2),[],1),[],2);
            pdcimgs_interp = abs(pdcimgs_interp).*exp(1i*angle(pdcimgs_interp.*conj(lpfiltered)));

            pdcimgs_interp = cat(3,zeros(iny,inx,1,nt,3),pdcimgs_interp,zeros(iny,inx,1,nt,3));

            ux = angle(pdcimgs_interp(:,:,:,:,1));
            uy = angle(pdcimgs_interp(:,:,:,:,2));
            uz = angle(pdcimgs_interp(:,:,:,:,3));

            if coordinateChoice==2
                ux=-ux;
            elseif coordinateChoice==3
                uy=-uy;
            elseif coordinateChoice==4
                uz=-uz;
            elseif coordinateChoice==5
                tmp=ux;
                ux=uy;
                uy=tmp;
                clear tmp
            elseif coordinateChoice==6
                uy=-uy;
                tmp=ux;
                ux=uy;
                uy=tmp;
                clear tmp
            elseif coordinateChoice==7
                ux=-ux;
                tmp=ux;
                ux=uy;
                uy=tmp;
                clear tmp
            elseif coordinateChoice==8
                uz=-uz;
                tmp=ux;
                ux=uy;
                uy=tmp;
                clear tmp
            end

            px_nomo = exp(1i*ux); py_nomo = exp(1i*uy); pz_nomo = exp(1i*uz);

        end

        maskt = repmat((mask & region),[1 1 1 nt]);
        px=px_nomo;
        py=py_nomo;
        pz=pz_nomo;

        % calculate duz/dy
        duzdy = zeros(iny,inx,inz,nt);
        ind = find(imerode(maskt,ones(3,1,1))); duzdy_temp = real(-1i*conj(pz).*convn(pz,dy,'same')); duzdy(ind) = duzdy_temp(ind); clear ind
        ind = find(~imerode(maskt,ones(3,1,1)) & imerode(maskt,abs(dy_top))); duzdy_temp = 2*real(-1i*conj(pz).*convn(pz,dy_bottom,'same')); duzdy(ind) = duzdy_temp(ind); clear ind
        ind = find(~imerode(maskt,ones(3,1,1)) & imerode(maskt,abs(dy_bottom))); duzdy_temp = 2*real(-1i*conj(pz).*convn(pz,dy_top,'same')); duzdy(ind) = duzdy_temp(ind); clear ind
        clear duzdy_temp

        % calculate duy/dz
        duydz = zeros(iny,inx,inz,nt);
        clear tmp
        ind = find(imerode(maskt,ones(1,1,3))); duydz_temp = real(-1i*conj(py).*convn(py,dz,'same')); duydz(ind) = duydz_temp(ind); clear ind
        ind = find(~imerode(maskt,ones(1,1,3)) & imerode(maskt,abs(dz_top))); duydz_temp = 2*real(-1i*conj(py).*convn(py,dz_bottom,'same')); duydz(ind) = duydz_temp(ind); clear ind
        ind = find(~imerode(maskt,ones(1,1,3)) & imerode(maskt,abs(dz_bottom))); duydz_temp = 2*real(-1i*conj(py).*convn(py,dz_top,'same')); duydz(ind) = duydz_temp(ind); clear ind
        clear duydz_temp

        % calculate dux/dz
        duxdz = zeros(iny,inx,inz,nt);
        ind = find(imerode(maskt,ones(1,1,3))); duxdz_temp = real(-1i*conj(px).*convn(px,dz,'same')); duxdz(ind) = duxdz_temp(ind); clear ind
        ind = find(~imerode(maskt,ones(1,1,3)) & imerode(maskt,abs(dz_top))); duxdz_temp = 2*real(-1i*conj(px).*convn(px,dz_bottom,'same')); duxdz(ind) = duxdz_temp(ind); clear ind
        ind = find(~imerode(maskt,ones(1,1,3)) & imerode(maskt,abs(dz_bottom))); duxdz_temp = 2*real(-1i*conj(px).*convn(px,dz_top,'same')); duxdz(ind) = duxdz_temp(ind); clear ind
        clear duxdz_temp

        % calculate duz/dx
        duzdx = zeros(iny,inx,inz,nt);
        ind = find(imerode(maskt,ones(1,3,1))); duzdx_temp = real(-1i*conj(pz).*convn(pz,dx,'same')); duzdx(ind) = duzdx_temp(ind); clear ind
        ind = find(~imerode(maskt,ones(1,3,1)) & imerode(maskt,abs(dx_left))); duzdx_temp = 2*real(-1i*conj(pz).*convn(pz,dx_right,'same')); duzdx(ind) = duzdx_temp(ind); clear ind
        ind = find(~imerode(maskt,ones(1,3,1)) & imerode(maskt,abs(dx_right))); duzdx_temp = 2*real(-1i*conj(pz).*convn(pz,dx_left,'same')); duzdx(ind) = duzdx_temp(ind); clear ind
        clear duzdx_temp

        %calculate duy/dx
        duydx = zeros(iny,inx,inz,nt);
        ind = find(imerode(maskt,ones(1,3,1))); duydx_temp = real(-1i*conj(py).*convn(py,dx,'same')); duydx(ind) = duydx_temp(ind); clear ind
        ind = find(~imerode(maskt,ones(1,3,1)) & imerode(maskt,abs(dx_left))); duydx_temp = 2*real(-1i*conj(py).*convn(py,dx_right,'same')); duydx(ind) = duydx_temp(ind); clear ind
        ind = find(~imerode(maskt,ones(1,3,1)) & imerode(maskt,abs(dx_right))); duydx_temp = 2*real(-1i*conj(py).*convn(py,dx_left,'same')); duydx(ind) = duydx_temp(ind); clear ind
        clear duydx_temp

        % calculate dux/dy
        duxdy = zeros(iny,inx,inz,nt);
        ind = find(imerode(maskt,ones(3,1,1))); duxdy_temp = real(-1i*conj(px).*convn(px,dy,'same')); duxdy(ind) = duxdy_temp(ind); clear ind
        ind = find(~imerode(maskt,ones(3,1,1)) & imerode(maskt,abs(dy_top))); duxdy_temp = 2*real(-1i*conj(px).*convn(px,dy_bottom,'same')); duxdy(ind) = duxdy_temp(ind); clear ind
        ind = find(~imerode(maskt,ones(3,1,1)) & imerode(maskt,abs(dy_bottom))); duxdy_temp = 2*real(-1i*conj(px).*convn(px,dy_top,'same')); duxdy(ind) = duxdy_temp(ind); clear ind
        clear duxdy_temp


        curlx_nm = duzdy - duydz;
        curly_nm = duxdz - duzdx;
        curlz_nm = duydx - duxdy;

        clear duzdy duydz duxdz duzdx duydx duxdy

        % create complex no motion imag
        curlx0 = curlx_nm(:,:,:,1) + 1i*curlx_nm(:,:,:,2);
        curly0 = curly_nm(:,:,:,1) + 1i*curly_nm(:,:,:,2);
        curlz0 = curlz_nm(:,:,:,1) + 1i*curlz_nm(:,:,:,2);

        brainmask = single(curlMask & region);
        % local noise is sqrt(2)*std(real and imaginary parts of complex no motion image)
        localNoise = zeros(iny,inx,inz);
        for ii = 2:(inx-1)
            for jj = 2:(inx-1)
                for kk = 2:(inz-1)
                    maskNeighborhood = brainmask(ii-1:ii+1,jj-1:jj+1,kk-1:kk+1);
                    maskNeighborhood=reshape(maskNeighborhood,27,1);
                    
                    datax=curlx0(ii-1:ii+1,jj-1:jj+1,kk-1:kk+1);
                    datax=reshape(datax,27,1);
                    datax(maskNeighborhood==0)=[];
                    
                    datay=curly0(ii-1:ii+1,jj-1:jj+1,kk-1:kk+1);
                    datay=reshape(datay,27,1);
                    datay(maskNeighborhood==0)=[];
                    
                    dataz=curlz0(ii-1:ii+1,jj-1:jj+1,kk-1:kk+1);
                    dataz=reshape(dataz,27,1);
                    dataz(maskNeighborhood==0)=[];
                    
                    localNoise(ii,jj,kk) = std([real(datax);imag(datax);real(datay);imag(datay);real(dataz);imag(dataz)]);
                end
            end
        end

        % snr is (amplitude of 1st harmonic of curl with motion)/(std of curl without motion)

        snrx = (convn(abs(curlxf),ones(3,3,3)/27,'same')./convn((curlMask & region),ones(3,3,3)/27,'same'))./localNoise;
        snrx(isnan(snrx))=0;
        snrx(isinf(snrx))=0;
        snry = (convn(abs(curlyf),ones(3,3,3)/27,'same')./convn((curlMask & region),ones(3,3,3)/27,'same'))./localNoise;
        snry(isnan(snry))=0;
        snry(isinf(snry))=0;
        snrz = (convn(abs(curlzf),ones(3,3,3)/27,'same')./convn((curlMask & region),ones(3,3,3)/27,'same'))./localNoise;
        snrz(isnan(snrz))=0;
        snrz(isinf(snrz))=0;

        % temporal ft of curl data...
        xf = 2 * ifft(curlx,[],4);
        xf = xf(:,:,:,2) .* brainmask;
        yf = 2 * ifft(curly,[],4);
        yf = yf(:,:,:,2) .* brainmask;
        zf = 2 * ifft(curlz,[],4);
        zf = zf(:,:,:,2) .* brainmask;

%         if (0)
%             %%% TESTING %%%
%             %%% Smooth curl data with ### NOTE: 2-D ### Bilateral filter (instead of 3-D Romano)...
%             radius = 3;
%             sigma = 0.77312;    % Equivalent to Romano 3...
%             filtered_xf = zeros(size(xf));
%             filtered_yf = zeros(size(yf));
%             filtered_zf = zeros(size(zf));
%             for zloop = 1:size(xf,3)
%                 filtered_xf(:,:,zloop) = BilateralFilter(xf(:,:,zloop), radius, sigma);
%             end
%             for zloop = 1:size(yf,3)
%                 filtered_yf(:,:,zloop) = BilateralFilter(yf(:,:,zloop), radius, sigma);
%             end
%             for zloop = 1:size(zf,3)
%                 filtered_zf(:,:,zloop) = BilateralFilter(zf(:,:,zloop), radius, sigma);
%             end
%         else
            % generate altered romano filter (entire 3x3x3 kernel is non-zero)
            filt = romano_filter(3,3,3);

            % Smooth curl data with edge-adaptive Romano filter...
            filtered_xf = (convn(xf,filt,'same') ./ convn(brainmask,filt,'same')) .* brainmask;
            filtered_yf = (convn(yf,filt,'same') ./ convn(brainmask,filt,'same')) .* brainmask;
            filtered_zf = (convn(zf,filt,'same') ./ convn(brainmask,filt,'same')) .* brainmask;
%         end
        
        filtered_xf(isnan(filtered_xf)) = 0;
        filtered_yf(isnan(filtered_yf)) = 0;
        filtered_zf(isnan(filtered_zf)) = 0;

        % Direct inversion of smoothed curl data...
        
        % Original DI (Now by calling di_JDT2() with DI_option = 1, below)...
        % [ vs2, cmu ] = di(fov, f, [1 1 1], 3, filtered_xf, filtered_yf, filtered_zf);
            % NOTE: [1 1 1] means NO data smoothing inside the di() function.
        
        JDT_MASK = ones(size(filtered_xf));
        [ vs2, cmu ] = di_JDT2(fov, f, [1 1 1], 3, JDT_MASK, DI_option, filtered_xf, filtered_yf, filtered_zf);
        
        % Final 1-layer erosion of masked regions...
        final_region = curlMask & region;
        if (rr ~= 0)
            final_region = imerode(final_region, str_elem);
        end
        
        ind = find(final_region);

        %-----------------------------------
        
        % Original pipeline median of squared-wavespeed values (NOT SNR-corrected)...
        med_vs2_stiffness(counter) = nanmedian(vs2(ind));

        % Median of magnitude-of-complex-modulus (shear-modulus) values...
        mag_cmu_array = abs(cmu(ind));
        med_cmu_mag(counter) = nanmedian(mag_cmu_array);

        if (CALC_CMEDIAN_MAG && ~strcmpi(region_name_str,'no_mask'))
            % Independent magnitude of iterative-median-of-complex-modulus
            % (Calculated with Josh Trzasko's cmedian.m)...
            cmu_mag_iterative_cmedian(counter) = cmedian(cmu(ind));
        end
        
        % MEAN and standard deviation of magnitude-of-complex-modulus (shear-modulus) values...
        mean_cmu_mag(counter) = nanmean(mag_cmu_array);
        std_cmu_mag(counter)  = nanstd(mag_cmu_array);
        
        % Independent median of storage-modulus values...
        real_cmu_array = real(cmu);
        med_cmu_real(counter) = nanmedian(real_cmu_array(ind));
        cmuR = med_cmu_real(counter);

        % Independent median of loss-modulus values...
        imag_cmu_array = imag(cmu);
        med_cmu_imag(counter) = nanmedian(imag_cmu_array(ind));
        cmuI = med_cmu_imag(counter);
        
        % Magnitude of complex-modulus calculated from independent ~R and ~I median values...
        cmu_mag_from_RImeds(counter) = (cmuR^2 + cmuI^2) ^ 0.5;

        if (CALC_MAG_FROM_RI_MEANS)
            % Magnitude of complex-modulus calculated from independent ~R and ~I MEAN values...
            cmuRmean = nanmean(real_cmu_array(ind));
            cmuImean = nanmean(imag_cmu_array(ind));
            cmu_mag_from_RImeans(counter) = (cmuRmean^2 + cmuImean^2) ^ 0.5;
        end
        
        % Squared-wavespeed calculated from independent ~R and ~I median values...
        vs2_stiffness_from_RImeds(counter) = (2*(cmuR.^2 + cmuI.^2))./(cmuR + (cmuR.^2 + cmuI.^2).^.5);

        noVoxels(counter) = length(ind);

        % Report median CMU magnitude and number of voxels in region...
        %%fprintf('[n=%d]  %.2f  (n=%d)\n', orig_region_voxels, med_cmu_mag(counter), noVoxels(counter));
        fprintf('%.2f  (%d)\n', med_cmu_mag(counter), noVoxels(counter));
        
        %-----------------------------------

        if (strcmpi(region_name_str,'cerebrum'))
            % Save ERODED cerebrum stiffness images to (one) .mat file...
            vs2_eroded = zeros(size(vs2));
            vs2_eroded(ind) = vs2(ind);
            vs2_eroded(isnan(vs2_eroded)) = 0;  %#ok<NASGU> % Change NaNs to zeros
            cmu_eroded = zeros(size(cmu));
            cmu_eroded(ind) = cmu(ind); %#ok<NASGU>
            if ~exist(resultsdir2,'dir')
                mkdir(resultsdir2)
            end
            
            try
                if (special_csf_fraction ~= 0)
                    save(fullfile(resultsdir2, [patientid_nospaces '_' date9(t1_imgInfo.StudyDate) '_s' num2str(motion_imgInfo.SeriesNumber) '_cerebrum_eroded_vs2_cmu_SPECIAL_CSF.mat']), 'vs2_eroded', 'cmu_eroded');
                else
                    save(fullfile(resultsdir2, [patientid_nospaces '_' date9(t1_imgInfo.StudyDate) '_s' num2str(motion_imgInfo.SeriesNumber) '_cerebrum_eroded_vs2_cmu.mat']), 'vs2_eroded', 'cmu_eroded');
                end
            catch
            end
        end

        if (strcmpi(region_name_str,'cerebellum'))
            % Save ERODED cerebellum stiffness images to (one) .mat file...
            vs2_eroded = zeros(size(vs2));
            vs2_eroded(ind) = vs2(ind);
            vs2_eroded(isnan(vs2_eroded)) = 0;  %#ok<NASGU> % Change NaNs to zeros
            cmu_eroded = zeros(size(cmu));
            cmu_eroded(ind) = cmu(ind); %#ok<NASGU>
            if ~exist(resultsdir2,'dir')
                mkdir(resultsdir2)
            end
            
            try
                if (special_csf_fraction ~= 0)
                    save(fullfile(resultsdir2, [patientid_nospaces '_' date9(t1_imgInfo.StudyDate) '_s' num2str(motion_imgInfo.SeriesNumber) '_cerebellum_eroded_vs2_cmu_SPECIAL_CSF.mat']), 'vs2_eroded', 'cmu_eroded');
                else
                    save(fullfile(resultsdir2, [patientid_nospaces '_' date9(t1_imgInfo.StudyDate) '_s' num2str(motion_imgInfo.SeriesNumber) '_cerebellum_eroded_vs2_cmu.mat']), 'vs2_eroded', 'cmu_eroded');
                end
            catch
            end
        end
        
        save_no_mask_images = 0;
        if (original_rlist(1) == 0)
            % User specified region '0' (no_mask)...
            save_no_mask_images = 1;
        end
        
        if (strcmpi(region_name_str,'no_mask'))
            if (save_no_mask_images)
                % Save full volume stiffness images to .mat file...
                if ~exist(resultsdir2,'dir')
                    mkdir(resultsdir2)
                end
                save(fullfile(resultsdir2, [patientid_nospaces '_' date9(t1_imgInfo.StudyDate) ...
                                           '_s' num2str(motion_imgInfo.SeriesNumber) '_no_mask_vs2.mat']), 'vs2');
                save(fullfile(resultsdir2, [patientid_nospaces '_' date9(t1_imgInfo.StudyDate) ...
                                           '_s' num2str(motion_imgInfo.SeriesNumber) '_no_mask_cmu.mat']), 'cmu');
            end

            % Always save TIV stiffness images to .mat file...
            vs2_tiv = tiv_mask .* vs2; %#ok<NASGU>
            cmu_tiv = tiv_mask .* cmu;  %#ok<NASGU>
            if ~exist(resultsdir2,'dir')
                mkdir(resultsdir2)
            end
            save(fullfile(resultsdir2, [patientid_nospaces '_' date9(t1_imgInfo.StudyDate) ...
                                        '_s' num2str(motion_imgInfo.SeriesNumber) '_tiv_vs2.mat']), 'vs2_tiv');
            save(fullfile(resultsdir2, [patientid_nospaces '_' date9(t1_imgInfo.StudyDate) ...
                                        '_s' num2str(motion_imgInfo.SeriesNumber) '_tiv_cmu.mat']), 'cmu_tiv');
        end
        
        if (save_cmu_files)
            % Save each region's cmu image to binary file...
            try
                JDT_DI = '';
                if (DI_option == 2)
                    JDT_DI = '_JDT_LS';
                elseif (DI_option == 3)
                    JDT_DI = '_JDT_MAD';
                end
                
                cmufilename = ['standard_' mc10(patientid_nospaces) '_' date9(t1_imgInfo.StudyDate) JDT_DI '_cmu_' rname '.cr4'];
                if (strncmpi(region_name_str,'roi',3))
                    % ROI region...
                    erode_str = '';
                    if (n_erode_roi ~= 0)
                        erode_str = ['_E' num2str(n_erode_roi)];
                    end
                    cmufilename = ['standard_' mc10(patientid_nospaces) '_' date9(t1_imgInfo.StudyDate) JDT_DI '_cmu_roi_' num2str(rr - 40) erode_str '.cr4'];
                end

                fprintf('--> saving cmu image to %s\n\n', cmufilename);
                cmu(isnan(cmu)) = 0;    % Change NaNs to zeros
                cmu_roi = zeros(size(cmu));
                cmu_roi(ind) = cmu(ind);
                write_binary_file(cmu_roi, cmufilename, 'single', 'n', 0, 1);
                if ~exist(resultsdir2,'dir')
                    mkdir(resultsdir2)
                end
                copyfile(fullfile(procdir, cmufilename),resultsdir2);
            catch
                disp('Error saving cmu image');
                return
            end
        end % (save_cmu_files)
            
        regionalSnrx(counter) = nanmedian(snrx(ind));
        regionalSnry(counter) = nanmedian(snry(ind));
        regionalSnrz(counter) = nanmedian(snrz(ind));

        binSize=.25;

        snrxArray=snrx(ind);
        [h b]=hist(snrxArray,0:binSize:200);
        regionalSnrModeX(counter) = mean(b(h==max(h)));
        [sortedH sortOrder]=sort(h,'descend');
        sortedB=b(sortOrder);
        sumOfVoxels=sortedH(1);
        cnt=1;
        while sumOfVoxels<(.5*sum(sortedH))
            sumOfVoxels=sumOfVoxels+sortedH(cnt+1);
            cnt=cnt+1;
        end
        xMedianOfMostLikelySnrs(counter) = median(snrxArray(snrxArray>(min(sortedB(1:cnt))-binSize/2) & snrxArray<(max(sortedB(1:cnt))+binSize/2)));
        clear h b sortedH sortedB sortOrder

        snryArray=snry(ind);
        [h b]=hist(snryArray,0:binSize:200);
        regionalSnrModeY(counter) = mean(b(h==max(h)));
        [sortedH sortOrder]=sort(h,'descend');
        sortedB=b(sortOrder);
        sumOfVoxels=sortedH(1);
        cnt=1;
        while sumOfVoxels<(.5*sum(sortedH))
            sumOfVoxels=sumOfVoxels+sortedH(cnt+1);
            cnt=cnt+1;
        end
        yMedianOfMostLikelySnrs(counter) = median(snryArray(snryArray>(min(sortedB(1:cnt))-binSize/2) & snryArray<(max(sortedB(1:cnt))+binSize/2)));
        clear h b sortedH sortedB sortOrder

        snrzArray=snrz(ind);
        [h b]=hist(snrzArray,0:binSize:200);
        regionalSnrModeZ(counter) = mean(b(h==max(h)));
        [sortedH sortOrder]=sort(h,'descend');
        sortedB=b(sortOrder);
        sumOfVoxels=sortedH(1);
        cnt=1;
        while sumOfVoxels<(.5*sum(sortedH))
            sumOfVoxels=sumOfVoxels+sortedH(cnt+1);
            cnt=cnt+1;
        end
        zMedianOfMostLikelySnrs(counter) = median(snrzArray(snrzArray>(min(sortedB(1:cnt))-binSize/2) & snrzArray<(max(sortedB(1:cnt))+binSize/2)));
        clear h b sortedH sortedB sortOrder

        tol=0.001;
        maxIters=50;
        iters=0;

        measuredVs2 = med_vs2_stiffness(counter);
        measuredSnrx = xMedianOfMostLikelySnrs(counter);
        measuredSnry = yMedianOfMostLikelySnrs(counter);
        measuredSnrz = zMedianOfMostLikelySnrs(counter);

        % initial guess
        [~, delta] = snr_correction_3med([.12 .12 .12], f, [40 40 40 8], [measuredSnrx measuredSnry measuredSnrz], measuredVs2, measuredVs2);
        corrVs2 = measuredVs2-delta;

        while abs(delta)>tol && iters<maxIters
            [~, delta] = snr_correction_3med([.12 .12 .12], f, [40 40 40 8], [measuredSnrx measuredSnry measuredSnrz], measuredVs2, corrVs2);
            corrVs2 = corrVs2-delta;
            iters=iters+1;
        end

        corrected_med_vs2_stiffness(counter)=corrVs2;
        numberIters(counter)=iters;

        clear pdcimgs cimgs ux uy uz px py pz
        counter=counter+1;
        %disp(['done with region ' num2str(rr)]);
    end
    
    delete(region_waitbar);  % Close the waitbar window
    pause(0.1);

    if ~exist(resultsdir,'dir')
        mkdir(resultsdir)
    end
    
    if length(dir(fullfile(resultsdir,'r*.nii'))) < 1
        copyfile(fullfile(procdir,'r*.nii'),resultsdir);
        copyfile(fullfile(procdir,'magn.nii'),resultsdir);
    end

    proctime=now;

    if (no_results == 0)
%         try
            % Create Comma-Separated-Values(CSV) results file...
            prefix = 'standard_mreResults_';
            if (DI_option == 2)
                prefix = [prefix 'JDT_LS_'];
            elseif (DI_option == 3)
                prefix = [prefix 'JDT_MAD_'];
            end
            
            if (special_csf_fraction ~= 0)
                csvFile = fullfile(procdir,[prefix patientid_nospaces '_scan' date9(t1_imgInfo.StudyDate) '_proc' datestr(proctime,'yyyy_mmdd_HHMM') '_SPECIAL_CSF' '.csv']);
            else
                csvFile = fullfile(procdir,[prefix patientid_nospaces '_scan' date9(t1_imgInfo.StudyDate) '_proc' datestr(proctime,'yyyy_mmdd_HHMM') '.csv']);
            end

            fid = fopen(csvFile,'w');

            % 
            [ ~, fname, ext] = fileparts(csvFile);
            final_csvFile = [fname ext];
            fprintf(fid,'%s\n', final_csvFile);
            fprintf(fid,'%s\n', directory);

            if (abs(proc_xspacing_mm - 3.0) > 0.01 || abs(proc_yspacing_mm - 3.0) > 0.01)
                fprintf(fid,'(%.3f mm scan)\n', proc_xspacing_mm);
            end
            
            if (special_csf_fraction ~= 0)
                fprintf(fid,'Special CSF fraction threshold = %.3f\n', special_csf_fraction);
            end
           
            if (DI_option == 2)
                fprintf(fid,'Using Josh''s least-squares(XYZ) DI\n');
            elseif (DI_option == 3)
                fprintf(fid,'Using Josh''s minimum-absolute-deviation (MAD) DI\n');
            end
            
            fprintf(fid,'%d Hz\n', f);
            fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s,', ...
                'Patient_ID','Scan_Date','Region_Number','Region_Name', ...
                'Corrected_Vs2_Stiff','Uncorrected_Vs2_Stiff','MEAN_CMU_Mag','STD_CMU_Mag','CMU_Mag_med');

            if (CALC_CMEDIAN_MAG)
                fprintf(fid,'%s,', 'CMU_Iter_cmed_Mag');
            end
            
            if (CALC_MAG_FROM_RI_MEANS)
                fprintf(fid,'%s,', 'CMU_Mag_from_RImeans');
            end

            fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,\n', ...
                'CMU_Real_med','CMU_Imag_med','CMU_Mag_from_RImeds','Vs2_stiffness_from_RImeds','Number_of_Voxels','Number_of_Iterations', ...
                'X_Median_SNR','Y_Median_SNR','Z_Median_SNR','X_Mode_SNR','Y_Mode_SNR','Z_Mode_SNR','X_MML_SNR','Y_MML_SNR','Z_MML_SNR');

            for ii=1:length(rlist)
                if (rlist(ii) == 0 && no_reg0)
                    % Do not include region '0' (no-mask) in spreadsheet report
                    % (region '0' was not specified in the original 'rlist')...
                else
                    % Write spreadsheet line for region rlist(ii)...
                    fprintf(fid,'%s,%s,%d,%s,', patientid_nospaces, t1_imgInfo.StudyDate, rlist(ii), regionNames{ii} );
                    fprintf(fid,'%7.14f,%7.14f,', corrected_med_vs2_stiffness(ii), med_vs2_stiffness(ii) );
                    fprintf(fid,'%7.14f,%7.14f,%7.14f,', mean_cmu_mag(ii), std_cmu_mag(ii), med_cmu_mag(ii) );

                    if (CALC_CMEDIAN_MAG)
                        fprintf(fid,'%7.14f,', cmu_mag_iterative_cmedian(ii) );
                    end

                    if (CALC_MAG_FROM_RI_MEANS)
                        fprintf(fid,'%7.14f,', cmu_mag_from_RImeans(ii) );
                    end

                    fprintf(fid,'%7.14f,%7.14f,', med_cmu_real(ii), med_cmu_imag(ii));
                    fprintf(fid,'%7.14f,%7.14f,', cmu_mag_from_RImeds(ii), vs2_stiffness_from_RImeds(ii));
                    fprintf(fid,'%d,%d,', noVoxels(ii), numberIters(ii));
                    fprintf(fid,'%7.14f,%7.14f,%7.14f,', regionalSnrx(ii), regionalSnry(ii), regionalSnrz(ii));
                    fprintf(fid,'%7.14f,%7.14f,%7.14f,', regionalSnrModeX(ii), regionalSnrModeY(ii), regionalSnrModeZ(ii));
                    fprintf(fid,'%7.14f,%7.14f,%7.14f,', xMedianOfMostLikelySnrs(ii), yMedianOfMostLikelySnrs(ii), zMedianOfMostLikelySnrs(ii));
                    fprintf(fid,'\n');
                end
            end

            if (~isempty(roi_filename))
                fprintf(fid,',\n,,roi filename:,%s\n', roi_filename);
            end

            fclose(fid);
            
            if ~exist(resultsdir2,'dir')
                mkdir(resultsdir2);
            end
            copyfile(csvFile,resultsdir2);
%         catch
%             try fclose(fid); catch ; end
%         end
    end
    
    % (Testing?)... save T1 image to .mat file...
    t1_filename = fullfile(resultsdir,['rM_' filename_idToSeries '_MT1_GW_N3m.nii']);
    if (exist(t1_filename, 'file'))
        a = load_nii(t1_filename);
        t1 = single(a.img); %#ok<NASGU>
        t1_mat_filename = fullfile(resultsdir,['rM_' filename_idToSeries '_MT1_GW_N3m.mat']);
        save(t1_mat_filename, 't1');
    end
    
    if (no_results)
        results = [];
    else
        % Feb 19, 2013: Added 'subject_id' and 'scan_date' fields to 'results' structure... 
        subject_id = cell(length(rlist), 1);
        subject_id{1} = mc10(patientid_nospaces);
        scan_date = cell(length(rlist), 1);
        scan_date{1} = date9(t1_imgInfo.StudyDate);
        for s = 2: length(rlist)
            subject_id{s} = subject_id{1};
            scan_date{s} = scan_date{1};
        end

        if (~isempty(roi_filename))
            results=struct( 'subject_id',{subject_id}, 'scan_date',{scan_date}, ...
                            'region_number',rlist', 'region_names',{regionNames}, ...
                            ...
                            'corrected_med_vs2_stiffness', corrected_med_vs2_stiffness, ...
                            'uncorrected_med_vs2_stiffness', med_vs2_stiffness, ...
                            'mean_cmu_mag', mean_cmu_mag, ...
                            'std_cmu_mag', std_cmu_mag, ...
                            'med_cmu_mag', med_cmu_mag, ...
                            'cmu_mag_iterative_cmedian', cmu_mag_iterative_cmedian, ...
                            'cmu_mag_from_RImeans', cmu_mag_from_RImeans, ...
                            'med_cmu_real', med_cmu_real, ...
                            'med_cmu_imag', med_cmu_imag, ...
                            'cmu_mag_from_RImeds', cmu_mag_from_RImeds, ...
                            'vs2_stiffness_from_RImeds', vs2_stiffness_from_RImeds, ...
                            ...
                            'number_voxels',noVoxels, 'number_iterations',numberIters,...
                            'x_median_snr',regionalSnrx, 'y_median_snr',regionalSnry, 'z_median_snr',regionalSnrz,...
                            'x_mode_snr',regionalSnrModeX, 'y_mode_snr',regionalSnrModeY, 'z_mode_snr',regionalSnrModeZ,...
                            'x_mml_snr',xMedianOfMostLikelySnrs, 'y_mml_snr',yMedianOfMostLikelySnrs, 'z_mml_snr',zMedianOfMostLikelySnrs, ...
                            'roi_filename',{roi_filename} ...
                          );
        else
            results=struct( 'subject_id',{subject_id}, 'scan_date',{scan_date}, ...
                            'region_number',rlist', 'region_names',{regionNames}, ...
                            ...
                            'corrected_med_vs2_stiffness', corrected_med_vs2_stiffness, ...
                            'uncorrected_med_vs2_stiffness', med_vs2_stiffness, ...
                            'mean_cmu_mag', mean_cmu_mag, ...
                            'std_cmu_mag', std_cmu_mag, ...
                            'med_cmu_mag', med_cmu_mag, ...
                            'cmu_mag_iterative_cmedian', cmu_mag_iterative_cmedian, ...
                            'cmu_mag_from_RImeans', cmu_mag_from_RImeans, ...
                            'med_cmu_real', med_cmu_real, ...
                            'med_cmu_imag',  med_cmu_imag, ...
                            'cmu_mag_from_RImeds', cmu_mag_from_RImeds, ...
                            'vs2_stiffness_from_RImeds', vs2_stiffness_from_RImeds, ...
                            ...
                            'number_voxels',noVoxels, 'number_iterations',numberIters,...
                            'x_median_snr',regionalSnrx, 'y_median_snr',regionalSnry, 'z_median_snr',regionalSnrz,...
                            'x_mode_snr',regionalSnrModeX, 'y_mode_snr',regionalSnrModeY, 'z_mode_snr',regionalSnrModeZ,...
                            'x_mml_snr',xMedianOfMostLikelySnrs, 'y_mml_snr',yMedianOfMostLikelySnrs, 'z_mml_snr',zMedianOfMostLikelySnrs ...
                          );
        end

        prefix = 'standard_mreResults_';
        if (DI_option == 2)
            prefix = [prefix 'JDT_LS_'];
        elseif (DI_option == 3)
            prefix = [prefix 'JDT_MAD_'];
        end
        
        if (special_csf_fraction ~= 0)
            matFile = fullfile(procdir,[prefix patientid_nospaces '_scan' date9(t1_imgInfo.StudyDate) '_proc' datestr(proctime,'yyyy_mmdd_HHMM') '_SPECIAL_CSF' '.mat']);
        else
            matFile = fullfile(procdir,[prefix patientid_nospaces '_scan' date9(t1_imgInfo.StudyDate) '_proc' datestr(proctime,'yyyy_mmdd_HHMM') '.mat']);
        end

        save(matFile,'-struct','results');
        if ~exist(resultsdir2,'dir')
            mkdir(resultsdir2);
        end
        copyfile(matFile,resultsdir2);
    end
    
    mls_cleanup_dir(procdir);
    
    try  cd(currentDir);  catch;  end  % Trap error, in case currentDir was moved/renamed/deleted.

    disp(' ');
    if (pause_flag == 1)
        disp('Standard_MRE_Pipeline completed');
        fprintf('(%s, %s, %s)\n\n', mc10(patientid_nospaces), date9(t1_imgInfo.StudyDate), directory);
    end
    
    t = toc(elapsed);
    
    % Set file attributes (all users) for all results files...
    set_results_files_rw(directory);
    set_results_files_rw(resultsdir);   % same as fullfile(directory,'pipeline_results')
    if ~strcmp(resultsdir,resultsdir2)
        % resultsdir and resultsdir2 are different directories...
        set_results_files_rw(resultsdir2);
    end
    
    if (user_wait_t == 0)
        fprintf('Processing time %.1f minutes\n\n', t / 60);
    else
        fprintf('Processing time %.1f minutes\n', (t - user_wait_t) / 60);
        fprintf('   Elapsed time %.1f minutes\n\n', t / 60);
    end
    
% catch
%     
%    % Error...
%    try delete(h);              catch; end
%    try delete(region_waitbar); catch; end
%    results = [];
%     
% end

end  % function Standard_MRE_Pipeline.m
