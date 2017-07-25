% demo script for splitting the field of view in patches and processing in parallel
% through memory mapping. See also run_pipeline.m for the complete  
% pre-processing pipeline of large datasets

clear;
%% load file

path_to_package = '../ca_source_extraction';   % path to the folder that contains the package
addpath(genpath(path_to_package));
             
%filename = '/Users/Brandon/Desktop/TestMotionCorrection/C1-018_rig.tif';      % path to stack tiff file
filename = '/Volumes/Samsung USB/Hollis lab_videos 2min recording_7.11.17/#717_7.11.17/#717_7.10.17-pm_002_rig.tif';
%foldername = '/Volumes/Samsung USB/Hollis lab_videos 2min recording_7.11.17/#717_7.11.17/#717_7.10.17-pm_002_rig';    % path to folder that contains a sequence of tiff files
%%
if exist([filename(1:end-3),'mat'],'file')
    data = matfile([filename(1:end-3),'mat'],'Writable',true);
else
    sframe=1;						% user input: first frame to read (optional, default 1)
    num2read=[];					% user input: how many frames to read   (optional, default until the end)
    chunksize=5000;                 % user input: read and map input in chunks (optional, default read all at once)
    data = memmap_file(filename,sframe,num2read,chunksize);
    %data = memmap_file_sequence(foldername);
end

%% Set parameters
sizY = size(data.Y);                  % size of data matrix
patch_size = [32,32];                   % size of each patch along each dimension (optional, default: [32,32])
overlap = [4,4];                        % amount of overlap in each dimension (optional, default: [4,4])

patches = construct_patches(sizY(1:end-1),patch_size,overlap);
K = 3;                                            % number of components to be found
tau = 6;                                          % std of gaussian kernel (size of neuron) 
p = 0;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.70;                                 % merging threshold

options = CNMFSetParms(...
    'd1',sizY(1),'d2',sizY(2),...
    'search_method','ellipse','dist',3,...      % search locations when updating spatial components
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps 
    'init_method','HALS',...
    'ssub',2,...
    'tsub',10,...
    'fudge_factor',0.98,...                     % bias correction for AR coefficients
    'merge_thr',merge_thr,...                   % merging threshold
    'gSig',tau... 
    );

%% Run on patches

[A,b,C,f,S,P,RESULTS,YrA] = run_CNMF_patches(data,K,patches,tau,p,options);

%% classify components
[ROIvars.rval_space,ROIvars.rval_time,ROIvars.max_pr,ROIvars.sizeA,keep] = classify_components(data.Y,A,C,b,f,YrA,options);

%% run GUI for modifying component selection (optional, close twice to save values)
Cn = reshape(P.sn,sizY(1),sizY(2));  % background image for plotting
run_GUI = false;
if run_GUI 
    Coor = plot_contours(A,Cn,options,1);
    GUIout = ROI_GUI(A,options,Cn,Coor,keep,ROIvars);   
    options = GUIout{2};
    keep = GUIout{3};    
end

%% Merge Accepted Components (optional, shouldn't have to do again)
%[Am,Cm,K_m,merged_ROIs,Pm,Sm] = merge_components(Yr,A(:,keep),b,C(keep,:),f,P,S,options);


%% view contour plots of selected and rejected components (optional)
if run_GUI
    throw = ~keep;
    figure;
        ax1 = subplot(121); plot_contours(A(:,keep),Cn,options,0,[],Coor,1,find(keep)); title('Selected components','fontweight','bold','fontsize',14);
        ax2 = subplot(122); plot_contours(A(:,throw),Cn,options,0,[],Coor,1,find(throw));title('Rejected components','fontweight','bold','fontsize',14);
        linkaxes([ax1,ax2],'xy')
end 
%% keep only the active components  

A_keep = A(:,keep);
C_keep = C(keep,:);

%% re-estimate temporal components
P.p = 2;      % perform deconvolution
[C2,f2,P2,S2,YrA2] = update_temporal_components_fast(data,A_keep,b,C_keep,f,P,options);


%% plot results
options.sx = 64;
plot_components_GUI(double(data.Y),A_keep,C2,b,f2,Cn,options);

[C_df,Df] = extract_DF_F_new(A_keep,C2,b,f2,P2,options);

%% Save Results
[dir_nm,file_nm,ext] = fileparts(filename);

[Coor, json_file] = plot_contours(A_keep,Cn,options,1);
saveas(gcf,strcat(dir_nm, '/' ,file_nm,'.png')); close;

nneurons=length(json_file);
c_df = C_df';


for i=1:nneurons
   [json_file(i).dff]=deal(full(c_df(:,i)));
end

file_save_loc = strcat(dir_nm, '/' ,file_nm,'.json');
%savejson('jmesh',json_file,file_save_loc);

%save(strcat(dir_nm, '/' ,file_nm,'_workspace','.mat'))

%% Video
%make_patch_video(A_keep,C2,b,f2,data.Yr,Coor,options)