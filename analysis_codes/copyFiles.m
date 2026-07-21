srcRoot = 'D:\freeview_photos_freesurfer\';
dstRoot = 'C:\Users\barak\Documents\Predictive_Ripples_2025\data\electrodes_viz';

for i = 1:numel(subjNames)

    subj = subjNames{i};

    srcFolder = fullfile(srcRoot, subj, 'elec_recon', 'PICS');
    dstFolder = fullfile(dstRoot, subj);

    if ~isfolder(srcFolder)
        fprintf('Brak folderu: %s\n', srcFolder);
        continue
    end

    if ~isfolder(dstFolder)
        mkdir(dstFolder);
    end

    copyfile(fullfile(srcFolder, '*'), dstFolder);

    fprintf('%s: skopiowano zawartość PICS\n', subj);

end


% Columns: data_name | folder_name | pic_name | elec_locs_name
subject_map = {
    'subNS128_02',    'NS128_02',    'NS128_02',    'NS128_02';
    'subNS140_02',    'NS140_02',    'NS140_02',    'NS140_02';
    'subNS167',       'NS167',       'NS167',       'NS167';
    'subNS173',       'NS173',       'NS173',       'NS173';
    'subNS174',       'NS174_02',    'NS174_02',    'NS174_02';
    'subNS174_03',    'NS174_03',    'NS174_03',    'NS174_03';
    'sublij135',      'NS135',       'NS135',       'NS135';
    'sublij136',      'NS136',       'NS136',       'NS136';
    'sublij137',      'NS137',       'NS137',       'NS137';
    'sublij138',      'NS138',       'NS138',       'NS138';
    'sublij140',      'NS140',       'NS140',       'NS140';
    'sublij2',        'NS134_impl2', 'NS134_impl2', 'NS134_impl2';
    'subns127_02',    'NS127_02',    'NS127_02',    'NS127_02';
    'subns142',       'NS142',       'NS142',       'NS142';
    'subns144_02',    'NS144_02',    'NS144_02',    'NS144_02';
    'subns148_02',    'NS148_02',    '',            'NS148_02';
    'subns149_02',    'NS149',       'NS149',       'NS149';
    'subns151',       'NS151',       'NS151',       'NS151';
    'subns153',       'NS153',       'NS153',       'NS153';
    'subns155',       'NS155',       'NS155',       'NS155';
    'subns166',       'NS166',       'NS166',       'NS166';
};