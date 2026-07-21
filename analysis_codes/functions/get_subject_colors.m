function cmap = get_subject_colors(subj_names)
% GET_SUBJECT_COLORS  Return consistent per-subject RGB colors.
%   subj_names — cell array of NS-style subject IDs (e.g. 'NS127_02')
%   Returns an n x 3 RGB matrix, one row per entry in subj_names.
%
% Colors are fixed by the master list order so the same subject always
% gets the same color regardless of how many others are present.
% Unknown subjects receive grey [0.5 0.5 0.5].

master_list = { ...
    'NS127_02','NS128_02','NS134_impl2','NS135','NS136','NS137', ...
    'NS138','NS140','NS140_02','NS142','NS144_02','NS148_02', ...
    'NS149','NS151','NS153','NS155','NS166','NS167','NS173', ...
    'NS174_02','NS174_03'};

% n_master = length(master_list);
% hues     = linspace(0, 1, n_master + 1);  hues = hues(1:end-1);
% sats     = repmat([0.90, 0.55], 1, ceil(n_master/2));  sats = sats(1:n_master);
% vals     = repmat([0.92, 0.72], 1, ceil(n_master/2));  vals = vals(1:n_master);
% palette  = hsv2rgb([hues(:), sats(:), vals(:)]);

palette = [
  0.894, 0.102, 0.110;  % red               NS127_02
  0.216, 0.494, 0.722;  % blue              NS128_02
  0.302, 0.686, 0.290;  % green             NS134_impl2
  0.596, 0.306, 0.639;  % purple            NS135
  1.000, 0.498, 0.000;  % orange            NS136
  1.000, 1.000, 0.200;  % yellow            NS137
  0.651, 0.337, 0.157;  % brown             NS138
  0.969, 0.506, 0.749;  % pink              NS140
  0.969, 0.506, 0.749;  % pink              NS140_02
  0.400, 0.761, 0.647;  % teal              NS142
  0.988, 0.553, 0.384;  % salmon            NS144_02
  0.553, 0.627, 0.796;  % lavender blue     NS148_02
  0.702, 0.886, 0.502;  % light green       NS149  
  0.800, 0.800, 0.200;  % olive             NS151
  0.992, 0.706, 0.384;  % peach             NS153
  0.200, 0.200, 0.600;  % dark blue         NS155
  0.600, 0.200, 0.200;  % dark red          NS166
  0.200, 0.600, 0.200;  % dark green        NS167
  0.745, 0.729, 0.855;  % lilac             NS173
  0.000, 0.749, 0.749;  % cyan              NS174_02
  0.749, 0.000, 0.749;  % magenta           NS174_03
];




n    = length(subj_names);
cmap = zeros(n, 3);
for i = 1:n
    idx = find(strcmp(master_list, subj_names{i}), 1);
    if isempty(idx)
        cmap(i, :) = [0.5 0.5 0.5];
    else
        cmap(i, :) = palette(idx, :);
    end
end
end
