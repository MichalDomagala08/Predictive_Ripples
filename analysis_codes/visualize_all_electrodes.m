% Visualize all subjects' electrodes on the fsaverage brain surface using iELVis
%
% Requires iELVis toolbox (D:\toolboxes\iELVis-master) and FreeSurfer
% subjects directory (E:\freesurfer) with fsaverage surface.

%% Paths
iELVis_path   = 'D:\toolboxes\iELVis-master';
% fs_dir        = 'E:\freesurfer';
fs_dir        = 'D:\freeview_photos_freesurfer'
elec_coord_dir = 'D:\UJ_more_data\projects\intracPhotos\data\elec_coordinates';

addpath(genpath(iELVis_path));

%% Collect electrodes across all subjects

subjects = {'NS128_02'};
        
        
% 
% subjects = {'NS127_02', 'NS138'};


all_coords = [];   % Nx3 RAS in fsaverage space
all_isLeft = [];   % Nx1 binary (1=left, 0=right)
all_names  = {};   % cell array of electrode labels (SubjectName_ElecName)

for si = 1:length(subjects)
    sub = subjects{si};

    % --- Load RAS coordinates (skip 2-line header) ---
    coord_file = fullfile(elec_coord_dir, [sub '.FSAVERAGE']);
    fid = fopen(coord_file, 'r');
    fgetl(fid); % timestamp line
    fgetl(fid); % "R A S" header
    coords_raw = textscan(fid, '%f %f %f');
    fclose(fid);
    coords = [coords_raw{1}, coords_raw{2}, coords_raw{3}];

    % --- Load electrode names and hemisphere labels ---
    names_file = fullfile(fs_dir, sub, 'elec_recon', [sub '.electrodeNames']);
    fid = fopen(names_file, 'r');
    fgetl(fid); % timestamp line
    fgetl(fid); % "Name, Depth/Strip/Grid, Hem" header
    names_raw = textscan(fid, '%s %s %s');
    fclose(fid);
    elec_names = names_raw{1};   % cell array of electrode name strings
    hemispheres = names_raw{3};  % cell array of 'L' or 'R'

    n_elec = size(coords, 1);
    if length(elec_names) ~= n_elec
        warning('%s: coordinate count (%d) != name count (%d), skipping.', sub, n_elec, length(elec_names));
        continue;
    end

    is_left = coords(:,1) < 0;   % negative x = left hemisphere in RAS

    % Prefix subject name to disambiguate across subjects
    prefixed_names = cellfun(@(n) [sub '_' n], elec_names, 'UniformOutput', false);

    all_coords = [all_coords; coords];
    all_isLeft = [all_isLeft; is_left];
    all_names  = [all_names; prefixed_names];
end

fprintf('Total electrodes loaded: %d from %d subjects\n', size(all_coords,1), length(subjects));

%% Build the Nx4 matrix expected by plotPialSurf (x, y, z, isLeft)
elec_matrix = [all_coords, double(all_isLeft)];

%% Assign a distinct color per subject
n_subjects = length(subjects);
cmap = get_subject_colors(subjects);

elec_colors = zeros(size(all_coords,1), 3);
idx = 1;
for si = 1:length(subjects)
    sub = subjects{si};
    % Count electrodes for this subject
    coord_file = fullfile(elec_coord_dir, [sub '.FSAVERAGE']);
    fid = fopen(coord_file, 'iELVis_pathr');
    fgetl(fid); fgetl(fid);
    raw = textscan(fid, '%f %f %f');
    fclose(fid);
    n = length(raw{1});
    elec_colors(idx:idx+n-1, :) = repmat(cmap(si,:), n, 1);
    idx = idx + n;
end

%% Plot — 2×2 subplots: dorsal (top row), medial (bottom row)
hfig = figure(1); clf;
set(hfig, 'Color', 'w', 'Position', [50 50 1400 900]);

panel_views  = {'l',             'r',              'lm',            'rm'};
panel_titles = {'Left lateral',  'Right lateral',  'Left medial', 'Right medial'};


cfg = [];
cfg.fsurfSubDir     = fs_dir;
cfg.elecCoord       = elec_matrix;
cfg.elecNames       = all_names;
cfg.elecColors      = elec_colors;
cfg.elecCbar        = 'n';   % required when elecColors is an RGB matrix
cfg.elecSize        = 2;
cfg.elecShape       = 'sphere';
cfg.ignoreDepthElec = 'n';
cfg.showLabels      = 'n';
cfg.opaqueness      = 0.2;
cfg.figId           = hfig;
cfg.clearFig        = 0;     % don't wipe the figure between subplots

ax = gobjects(4,1);
for vi = 1:4
    ax(vi) = subplot(2, 2, vi);
end

for vi = 1:4
    cfg.axis  = ax(vi);
    cfg.view  = panel_views{vi};
    cfg.title = panel_titles{vi};
    plotPialSurf('fsaverage', cfg);
end

% Tighten spacing between subplots
gap_h         = -0.2;
gap_v         =  0.1;
margin_side   =  0.04;
margin_top    =  0.12;
margin_bottom =  0.08;
w = (1 - 2*margin_side - gap_h) / 2;
h = (1 - margin_top - margin_bottom - gap_v) / 2;
positions = [
    margin_side,          margin_bottom+h+gap_v,  w, h;
    margin_side+w+gap_h,  margin_bottom+h+gap_v,  w, h;
    margin_side,          margin_bottom,           w, h;
    margin_side+w+gap_h,  margin_bottom,           w, h;
];
for vi = 1:4
    set(ax(vi), 'Position', positions(vi,:));
end

%% Subject-color legend (separate figure)
hleg = figure(2); clf;
set(hleg, 'Color', 'w', 'Position', [1500 50 300 40*n_subjects+60]);
ax_leg = axes(hleg);
hold(ax_leg, 'on');
h_leg = gobjects(n_subjects, 1);
for si = 1:n_subjects
    h_leg(si) = scatter(nan, nan, 100, cmap(si,:), 'filled', ...
                        'DisplayName', strrep(subjects{si}, '_', '\_'));
end
legend(h_leg, 'Location', 'west', 'FontSize', 12);
axis(ax_leg, 'off');

% %% Save figures
% out_dir = 'D:\UJ\projects\intracPhotos\figs';
% savefig(hfig, fullfile(out_dir, 'electrodes_brain.fig'));
% exportgraphics(hfig, fullfile(out_dir, 'electrodes_brain.png'), 'Resolution', 300);
% savefig(hleg, fullfile(out_dir, 'electrodes_legend.fig'));
% exportgraphics(hleg, fullfile(out_dir, 'electrodes_legend.png'), 'Resolution', 300);
