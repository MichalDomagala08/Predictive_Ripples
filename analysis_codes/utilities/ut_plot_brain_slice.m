
% This scritp plots a brain slices images with selected coloring: 


%%% --- SETUP AND PARAMETERS

params.region = "hippocampus"; % Select a region

%%% Paths
iELVis_root     = 'D:\Documents_Dell\Predictive_Ripples_2025\toolboxes\iELVis-master'; % iELVIS package path
subjects_parent = 'D:\Documents_Dell\Predictive_Ripples_2025\data\freeview_photos_freesurfer'; % Freesurfer image data Path
dataFolder      = 'D:\Documents_Dell\Predictive_Ripples_2025\data\reref'; %Folder where preprocessed Electrophysiological data are bieng housed
image_dir       = 'D:\Documents_Dell\Predictive_Ripples_2025\data\electrodes_viz';
if exist(iELVis_root,'dir'); addpath(genpath(iELVis_root)); end
setenv('SUBJECTS_DIR', subjects_parent);
setenv('FREESURFER_HOME', 'C:\Program Files\FreeSurfer');

%%% Get files
files = {dir(fullfile(dataFolder, '*.mat')).name};
subjNames = cellfun(@(x) erase(extractAfter(x, 'sub'), '.mat'), files, 'UniformOutput', false);
sbn = {subjNames{13:13}};


set(0, 'DefaultFigureVisible', 'off');

for i = 1:length(sbn)

    % Local Paths
    currSubjName = sbn{i};
    subjDir = fullfile(subjects_parent, currSubjName);
    elec_recon_dir = fullfile(subjDir,'elec_recon');
    subj_img_dir =  fullfile(image_dir,currSubjName);
    
    % quick checks of important Freesurfer files and paths
    assert(isfolder(subjDir),'Subject dir not found: %s', subjDir);
    assert(isfile(fullfile(subjDir,'mri','orig.mgz')),'orig.mgz missing');
    assert(isfile(fullfile(subjDir,'mri','aparc+aseg.mgz')),'aparc+aseg.mgz missing');
    assert(isfolder(elec_recon_dir),'elec_recon folder missing');
    disp('Basic files present. Proceed to iELVis plotting commands.');
    
    % slice plotting for depth electrodes (with anat overlay option)
    cfg = [];
    cfg.printFigs = 1;            % saves JPGs
    cfg.anatOverlay = 1;          % overlay aparc+aseg coloration on MRI (hipokamp will appear)
    cfg.fsurfSubDir = subjects_parent;
    cfg.region = params.region;
    cfg.printFigs = subj_img_dir;
    
    plot_all_depth_on_slices_mod(currSubjName, 'mgrid', cfg);
    % Output: elec_recon/PICS/*.jpg (slices per depth electrode)

end



