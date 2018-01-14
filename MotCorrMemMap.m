clear all classes;
addpath(genpath('../../../ca_source_extraction'));  % add packages to matlab path
addpath(genpath('../../../NoRMCorre'));
%% complete pipeline for calcium imaging data pre-processing
clear;
addpath(genpath('../NoRMCorre'));               % add the NoRMCorre motion correction package to MATLAB path
gcp;        % start a parallel engine
foldername = '/Volumes/Samsung USB/Hollis lab_videos 2min recording_7.11.17/#717_7.11.17/';% #717_7.10.17-pm_002.tif';   
        % folder where all the files are located. Currently supported .tif,
        % .hdf5, .raw, .avi, and .mat files
files = subdir(fullfile(foldername,'*.tif'));   % list of filenames (will search all subdirectories)
FOV = [512,512]; % Pixel width and height of video
numFiles = length(files);

%% motion correct and save as tif files, as well as saving the shifts and templates as .mat files
% register files one by one. use template obtained from file n to
% initialize template of file n + 1;
% Also memmaps the rigid files to .mat files for faster initialization
% later.

motion_correct = true;      % perform motion correction 
non_rigid = false;          % Perform non-rigid motion correction otherwise rigid motion correction
template = [];
for i = 1:numFiles
    fullname = files(i).name;
    [folder_name,file_name,ext] = fileparts(fullname);
    if motion_correct && isempty(strfind(file_name,'_rig')) && isempty(strfind(file_name,'_nr'))
        if non_rigid
            if ~exist(fullfile(folder_name,[file_name,'_nr.tif']),'file')
                options_nonrigid = NoRMCorreSetParms('d1',FOV(1),'d2',FOV(2),'grid_size',[128,128],...
                    'overlap_pre',64,'mot_uf',4,'bin_width',300,'max_shift',24,'max_dev',8,'us_fac',50,...
                    'output_type','tif','tiff_filename',fullfile(folder_name,[file_name,'_nr.tif']));
                [M,shifts,template] = normcorre_batch(fullname,options_nonrigid,template);
                save(fullfile(folder_name,[file_name,'_shifts_nr.mat']),'shifts','-v7.3');
                save(fullfile(folder_name,[file_name,'_template_nr.mat']),'template','-v7.3');
            end
            if ~exist(fullfile(folder_name,[file_name,'_nr.mat']),'file')
                data = memmap_file(options_nonrigid.tiff_filename);
            end
        else  
            if ~exist(fullfile(folder_name,[file_name,'_rig.tif']),'file')
                % perform rigid motion correction (faster, could be less accurate)
                options_rigid = NoRMCorreSetParms('d1',FOV(1),'d2',FOV(2),'bin_width',300,'max_shift',32,...
                'output_type','tif','tiff_filename',fullfile(folder_name,[file_name,'_rig.tif']));
                [M,shifts,template] = normcorre_batch(fullname,options_rigid,template);
                save(fullfile(folder_name,[file_name,'_shifts_rig.mat']),'shifts','-v7.3');
                save(fullfile(folder_name,[file_name,'_template_rig.mat']),'template','-v7.3');
            end
            if ~exist(fullfile(folder_name,[file_name,'_rig.mat']),'file')
                data = memmap_file(options_rigid.tiff_filename);
            end
        end      
    end
end  