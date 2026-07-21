ft_defaults


currSubjName = 'NS128_02';
talDir = 'C:\Users\barak\Documents\Predictive_Ripples_2025\data\elec_coordinates';

% Chosen Electrode for now:
coord = [15.5786200000000,	-6.08913100000000,	-22.8181510000000];
%coord = [12.3023160000000,	-6.65930900000000,	54.9207670000000];

d = dir(fullfile(talDir, ['*' currSubjName '*xfm']));
if isempty(d)
    error('No talairach xfm found for %s in %s', currSubjName, talDir);
end
found = fullfile(d(1).folder, d(1).name);
fprintf('Using transform file: %s\n', found);

% read_talxfm must return 3x4 row-wise (you already have this function)
M3x4 = read_talxfm(found);
T_fs2mni = [ M3x4; 0 0 0 1 ]; % make homogeneous 4x4


candA = (T_fs2mni * [coord(:);1])';
candB = (inv(T_fs2mni) * [coord(:);1])';

disp('candA ='); disp(candA);
disp('candB ='); disp(candB);


coord_mni =  candB(1:3); % ft_warp_apply(T_fs2mni, coord, 'homogeneous');

%%% Get MNE Template 

ftdir = fileparts(which('ft_defaults'));

mri = ft_read_mri(fullfile(ftdir,'template','anatomy','single_subj_T1_1mm.nii'));
aal = ft_read_atlas(fullfile(ftdir,'template','atlas','aal','ROI_MNI_V4.nii'));


% Create A  hippocampus Mask:
cfg = [];
cfg.atlas = aal;
cfg.roi = {'Hippocampus_L','Hippocampus_R'};
cfg.inputcoord = 'mni';

hippo = ft_volumelookup(cfg,aal); % Looks for specific volume

hippoVol = aal;

hippoVol.tissue = hippo;
hippoVol.dim = size(hippo);


% Resmaple to MRI
cfg = [];
cfg.parameter = 'tissue';
cfg.interpmethod = 'nearest';

hippoMRI = ft_sourceinterpolate(cfg,hippoVol,mri);

size(hippoMRI.tissue)
size(mri.anatomy)



vox = ft_warp_apply(inv(mri.transform),coord_mni,'homogeneous'); % Transfrom all to MNI coordinates

vox = round(vox);

%%% Extract Slices:

% Axial
MRI_ax = squeeze(mri.anatomy(:,:,vox(3)));
HIP_ax = squeeze(hippoMRI.tissue(:,:,vox(3)));

%Coronal
MRI_cor = squeeze(mri.anatomy(:,vox(2),:));
HIP_cor = squeeze(hippoMRI.tissue(:,vox(2),:));

% Sagital
MRI_sag = squeeze(mri.anatomy(vox(1),:,:));
HIP_sag = squeeze(hippoMRI.tissue(vox(1),:,:));

figure('Color','w','Position',[100 100 1400 450])


%% AXIAL
subplot(1,3,1)

imagesc(rot90(MRI_ax))
colormap gray
axis image off
hold on

contour(rot90(HIP_ax),[1 1],'r','LineWidth',2)

plot(size(MRI_ax,2)-vox(2),vox(1),...
    'bo','MarkerFaceColor','b','MarkerSize',10)

title('Axial')


%% CORONAL
subplot(1,3,2)

imagesc(rot90(MRI_cor))
colormap gray
axis image off
hold on

contour(rot90(HIP_cor),[1 1],'r','LineWidth',2)

plot(size(MRI_cor,2)-vox(3),vox(1),...
    'bo','MarkerFaceColor','b','MarkerSize',10)

title('Coronal')


%% SAGITTAL
subplot(1,3,3)

imagesc(rot90(MRI_sag))
colormap gray
axis image off
hold on

contour(rot90(HIP_sag),[1 1],'r','LineWidth',2)

plot(size(MRI_sag,2)-vox(3),vox(2),...
    'bo','MarkerFaceColor','b','MarkerSize',10)

title('Sagittal')