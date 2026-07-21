% Przygotowanie (EDYTUJ ścieżkę do iELVis jeśli trzeba)
iELVis_root = 'D:\Documents_Dell\Predictive_Ripples_2025\toolboxes\iELVis-master'; % <- ustaw jeśli iELVis nie jest na MATLAB path
if exist(iELVis_root,'dir')
    addpath(genpath(iELVis_root));
end

subjects_parent = 'D:\Documents_Dell\Predictive_Ripples_2025\data\freeview_photos_freesurfer';



% set FreeSurfer SUBJECTS_DIR environment variable so iELVis finds subjects non-interactively
setenv('SUBJECTS_DIR', subjects_parent);
setenv('FREESURFER_HOME', 'C:\Program Files\FreeSurfer');
dataFolder =     'D:\Documents_Dell\Predictive_Ripples_2025\data\reref'; %Folder where preprocessed Electrophysiological data are bieng housed


files = {dir(fullfile(dataFolder, '*.mat')).name};
subjNames = cellfun(@(x) erase(extractAfter(x, 'sub'), '.mat'), files, 'UniformOutput', false);


set(0, 'DefaultFigureVisible', 'off');

sbn = {subjNames{13:13}};
for i = 1:length(sbn)
    currSubjName = sbn{i};
    


    subjDir = fullfile(subjects_parent, currSubjName);
    elec_recon_dir = fullfile(subjDir,'elec_recon');
    subj_img_dir =  fullfile('D:\Documents_Dell\Predictive_Ripples_2025\data\electrodes_viz',currSubjName);
    
    % quick checks
    assert(isfolder(subjDir),'Subject dir not found: %s', subjDir);
    assert(isfile(fullfile(subjDir,'mri','orig.mgz')),'orig.mgz missing');
    assert(isfile(fullfile(subjDir,'mri','aparc+aseg.mgz')),'aparc+aseg.mgz missing');
    assert(isfolder(elec_recon_dir),'elec_recon folder missing');
    disp('Basic files present. Proceed to iELVis plotting commands.');
    
    % slice plotting for depth electrodes (with anat overlay option)
    cfg = [];
    cfg.printFigs = 1;            % saves JPGs
    cfg.anatOverlay = 1;         % overlay aparc+aseg coloration on MRI (hipokamp will appear)
    cfg.fsurfSubDir = subjects_parent;
    cfg.region = 'hippocampus';
    cfg.printFigs = subj_img_dir;
    
    plotAllDepthsOnSlices_mod(currSubjName, 'mgrid', cfg);
    % Output: elec_recon/PICS/*.jpg (slices per depth electrode)

end



