%PR_3_RippleAnalysis
%
% This script plots ripple densities and computes aggregated measures and
% table for later testings in case of RIPPLES:

% This scirpt in particular creates several tangible tables
% (1) Ripple_Channelwise Table: An old tabl having aggregated rippples
% channelwise
% (2) Ripplewise table: A new ripple containing in each row a ripple
% measurment
% (3) Saccadewise table: a supplemented original table of saccade counts

% (4) Computes Ripple densities
% (5) Computes saccadic densities




%%% --- PARAMETERS and SETUP --- %%%

%%% This script Finaly analyses the Ripples that are in a specific thing:
params.ripple_detection_scheme_name = "RippleDetection_franz_et_al_RipplePeak_manualCorr";
params.analysisName  = "RippleDensity_Franz_improved_saccadic";
params.artif_chan_rej     = true; % If we reject artefactual channels
params.blink_saccades_rej = true; % if we are rejecting blink saccades
params.sameRipples        = false; % removing ripples that coincide in time across channels
params.trialNum = 80; % How many trials are there
params.viz = 0;

addpath('D:\Documents_Dell\Predictive_Ripples_2025\analysis_codes\functions'); addpath("utilities\")
dataFolder =  fullfile("D:\Documents_Dell\Predictive_Ripples_2025\preprocessed",params.ripple_detection_scheme_name);
saveFolder =  fullfile("D:\Documents_Dell\Predictive_Ripples_2025\Results",params.analysisName); mkdir(saveFolder); 
saveDataFolder = fullfile("D:\Documents_Dell\Predictive_Ripples_2025\preprocessed",params.analysisName); mkdir(saveDataFolder)
imageFolder = "D:\Documents_Dell\Predictive_Ripples_2025\data\stimuli";
saccadeFile = 'D:\Documents_Dell\Predictive_Ripples_2025\data\saccadic_data\content_all_subjects_improved.mat';

%%% Save parameters in a log file:
T_params = struct2table(params, 'AsArray', true);
logFile = fullfile(saveFolder, 'processing_parameters.txt');
writetable(stack(T_params, 1:width(T_params)), logFile, ...
    'WriteVariableNames', false, 'Delimiter', '\t');

%%% Load data
ripple_data = load(fullfile(dataFolder,"rippleData.mat"));

%%% Load image dirs::
images  = {dir(fullfile(imageFolder,'*.png')).name};

%%% Load saccade file:
saccades_table = load(saccadeFile).content;

%%% Prepare saccadic info: 

%% 
%%% --- PREPROCESSING --- %%%

%%% Saccade Information
info_columns = [20:23, 27:38];
info_column_names = saccades_table(:,[20:23, 27:38]).Properties.VariableNames;

saccade_timestamps_all = cell(length([20:23, 27:38]),19,params.trialNum);

%%
% Change Subject Names to new ones (for proper synergy)

ut_subject_map
[found, idx] = ismember(string(saccades_table.subject), string(subject_map(:,1)));
saccades_table.subject_new = string(saccades_table.subject);             % domyślnie stara nazwa
saccades_table.subject_new(found) = string(subject_map(idx(found), 3));  % podmiana wg mapy

%%
%%% Saccade Kernel Densities
saccade_density_sbj = zeros(length(info_columns),length(ripple_data.all_subjNames),7001);
saccade_density_tri = zeros(length(info_columns),80,7001);
s_xi = zeros(length(info_columns),80,7001);

%%% Ripple Density Preprocessing

if params.artif_chan_rej; ut_bad_channels; else; artfi_subj_chan = cell(1,19);end % get bad channels

ripple_timestamps_all = cell(1,params.trialNum);
ripple_densities_per_subj = cell(19,params.trialNum);

%%% Ripple Kernel Densirties
ripple_density_tri = zeros(length(ripple_timestamps_all), 100);
ripple_density_timecourse = zeros(length(ripple_timestamps_all),100);
ripple_density_tri_x_sb = zeros(length(ripple_timestamps_all),length(ripple_data.all_subjNames),7001);
ripple_density_sbj = zeros(length(ripple_data.all_subjNames), 100);

%%% (1) Prelocating ripple channel table data:
subject = {}; subject_n= {};trial = {};channel = {};channel_name = {};
ripple_number = {};ripple_length  = {};sacc_number = {}; blink_trial  = {};
mean_sacc_velocity = {};mean_sacc_amplitude = {};mean_sacc_latency= {};
mean_sacc_duration = {};X_mean  = {};Y_mean = {};X_std = {};Y_std = {};
aws50px_mean = {};aws50px_std  = {};lum50px_mean = {};lum50px_std  = {};
ripple_length_cl = {}; ripple_number_cl={};

%%% (2)  Adding ripple features to saccades_table data 
saccades_table.closestRipple = zeros(size(saccades_table,1),1); saccades_table.closestRipple_TimingDiff = zeros(size(saccades_table,1),1);
saccades_table.closestRipple_Length = zeros(size(saccades_table,1),1); saccades_table.rippleNumIn200ms = zeros(size(saccades_table,1),1); 
saccades_table.rippleNumIn200ms_after = zeros(size(saccades_table,1),1); saccades_table.rippleNumIn200ms_before = zeros(size(saccades_table,1),1); 
saccades_table.closestRipple_idx = zeros(size(saccades_table,1),1);  saccades_table.closestRipple_channel = strings(height(saccades_table), 1);  

saccades_table.fixation_beg = zeros(size(saccades_table,1),1); saccades_table.fixation_end = zeros(size(saccades_table,1),1);
saccades_table.fixation_len= zeros(size(saccades_table,1),1); saccades_table.ripples_in_fixation = zeros(size(saccades_table,1),1); 
saccades_table.closest_rip_to_fix_beg_diff = zeros(size(saccades_table,1),1); saccades_table.closest_rip_to_fix_beg = zeros(size(saccades_table,1),1);
saccades_table.closest_rip_to_fix_half_diff = zeros(size(saccades_table,1),1); saccades_table.closest_rip_to_fix_half = zeros(size(saccades_table,1),1);

saccades_table.closestRipple_random_Length  = zeros(size(saccades_table,1),1); saccades_table.closestRipple_random  = zeros(size(saccades_table,1),1);
saccades_table.rippleNumBefore_saccade  = zeros(size(saccades_table,1),1); saccades_table.rippleNumIn_saccade  = zeros(size(saccades_table,1),1);
saccades_table.rippleNumAfter_saccade  = zeros(size(saccades_table,1),1); saccades_table.rippleNumAfter2_saccade  = zeros(size(saccades_table,1),1);
saccades_table.rippleNumAfter3_saccade  = zeros(size(saccades_table,1),1); saccades_table.rippleNumAfter4_saccade  = zeros(size(saccades_table,1),1);

saccades_table.rippleNumBefore_saccade_M  = zeros(size(saccades_table,1),1); saccades_table.rippleNumIn_saccade_M  = zeros(size(saccades_table,1),1);
saccades_table.rippleNumAfter_saccade_M  = zeros(size(saccades_table,1),1); saccades_table.rippleNumAfter2_saccade_M  = zeros(size(saccades_table,1),1);
saccades_table.rippleNumAfter3_saccade_M  = zeros(size(saccades_table,1),1); saccades_table.rippleNumAfter4_saccade_M  = zeros(size(saccades_table,1),1);

saccades_table.rippleAroundSaccade = zeros(size(saccades_table,1),1); %If Closest Ripple has been encountered, around 100ms around saccade onset/offset
saccades_table.rippleAfterSaccade = zeros(size(saccades_table,1),1); saccades_table.rippleBeforeSaccade = zeros(size(saccades_table,1),1);


% random index used in a ripple precedence computation and as a light permutation
rng(0);
saccades_table.random_sac_timing_index = 6*rand(height(saccades_table),1);

%%% Ripple in Saccades distribution vol 2:
% Get bin boundaries:
binBoundaries = [sort(-(13:round(nanmean( saccades_table.duration/2)):501)) 13:round(nanmean( saccades_table.duration/2)):501];

saccRippTiming = zeros(height(saccades_table),length(binBoundaries));


%%% (3) Get All Ripple distribution per channel as a building block for ripple table 
allRipples_timecourse_across_channels_all = cell(1,19);
allRipples_across_channels_channames_all = cell(1,19);
allRipples_across_channels_all = cell(1,19);


unfolded_subject = []; unfolded_channels = [];unfolded_beg = []; unfolded_repeat = [];  % Creating ripple wise table
unfolded_trials = []; unfolded_median=[]; unfolded_rip_num = []; unfolded_end = []; unfolded_times = [];

for sb = 1:length(ripple_data.all_subjNames )
    currentRippleSubj = ripple_data.all_hipothetical_ripples_cluster{sb};


    %%% manually omti subject 10 (no channels) - im tired of this shit,no
    %%% failsafe is working here XD
    if sb == 10; continue; end

    %%%% (1) (3) ripple per- channel computations
    allRipples_across_channels = cell(params.trialNum,1);
    allRipples_across_channels_channames = cell(params.trialNum,1);
    allRipples_timecourse_across_channels =  cell(1,params.trialNum);


    for chan = 1:length(ripple_data.all_chanNames{sb})

        %%%% (1) (3) ripple per- channel computations
        % Remove bad Ripple channels
        currentRippleSubjChan = currentRippleSubj{chan};
        if isempty(currentRippleSubjChan);  continue; end
        if ~isempty(artfi_subj_chan{sb})
            if sum(ismember(ripple_data.all_chanNames{sb}{chan},artfi_subj_chan{sb})); continue;  end
        end

        % Get channel name nfo
        currentRippleSubjChanname = currentRippleSubjChan;
        idx = find(cellfun(@(c) iscell(c) && ~isempty(c), currentRippleSubjChanname));
        for k = 1:length(idx), currentRippleSubjChanname{idx(k)} = repmat({ripple_data.all_chanNames{sb}{chan}}, size(currentRippleSubjChanname{idx(k)})); end

        allRipples_across_channels = cellfun(@(a,b) [a b], allRipples_across_channels, currentRippleSubjChan, 'UniformOutput', false);
        allRipples_across_channels_channames = cellfun(@(a,b) [a b], allRipples_across_channels_channames, currentRippleSubjChanname, 'UniformOutput', false);
        allRipples_timecourse_across_channels = cellfun(@(a,b)[a; b], allRipples_timecourse_across_channels, ripple_data.all_timecourseSWR_trial{sb}{chan}, 'UniformOutput', false);


        
        for it = 1:params.trialNum

            %%% (4) Get ripple densities
            sel = strcmp(saccades_table.subject_new, ripple_data.all_subjNames{sb}) & saccades_table.trial_number == it;
            currentSaccades = saccades_table(sel, :);
            ripple_timestamps_all{it} = [ripple_timestamps_all{it} cellfun(@median,currentRippleSubjChan{it})];
            ripple_densities_per_subj{sb,it}  = [ripple_densities_per_subj{sb,it}  cellfun(@median,currentRippleSubjChan{it})];

            %%% Ripple Kernel Desnity per Subject x Trial
            if ~isempty(ripple_densities_per_subj{sb,it})
                [ripple_density_tri_x_sb(it,sb,:), ~] = ksdensity(ripple_densities_per_subj{sb,it},1:7001,  'Bandwidth', 200); % density with kernel estimataion
            else;ripple_density_tri_x_sb(it,sb,:) = zeros(1,7001);end
            % 
            % % 
            % %%%% (1) Save Ripple Features to a table:
            % % Ripple Feaat and Headers
            % subject{end+1} =  ripple_data.all_subjNames{sb};  subject_n{end+1} = sb; trial{end+1} =  it; 
            % channel{end+1} = chan; channel_name{end+1} = ripple_data.all_chanNames{sb}{chan};
            % ripple_number{end+1} = length(currentRippleSubj{chan}{it}); 
            % ripple_length{end+1} = mean(cellfun(@length, currentRippleSubj{chan}{it}));
            % 
            % % cleaned ripple features (not in the first 500ms after trial start
            % ripple_number_cl{end+1} = length({currentRippleSubj{chan}{it}{find(cellfun(@median,currentRippleSubj{chan}{it}) >2250)'}}); 
            % ripple_length_cl{end+1} = mean(cellfun(@length, {currentRippleSubj{chan}{it}{find(cellfun(@median,currentRippleSubj{chan}{it}) >2250)'}}));
            % 
            % % Saccading features
            % sacc_number{end+1} = size(currentSaccades,1); blink_trial{end+1} = currentSaccades.blinkTrial(1);
            % mean_sacc_velocity{end+1} = mean(currentSaccades.peakVelocity,"omitnan");
            % mean_sacc_amplitude{end+1} = mean(currentSaccades.amplitude,"omitnan");
            % mean_sacc_latency{end+1} = mean(currentSaccades.latency,"omitnan");
            % mean_sacc_duration{end+1} = mean(currentSaccades.duration,"omitnan");
            % X_mean{end+1} =  mean(currentSaccades.Xpx_original,"omitnan");  Y_mean{end+1} =mean(currentSaccades.Ypx_original,"omitnan"); 
            % X_std{end+1} =  std(currentSaccades.Xpx_original,"omitnan");  Y_std{end+1} =std(currentSaccades.Ypx_original,"omitnan"); 
            % 
            % % Image Features 
            % aws50px_mean{end+1} = mean(currentSaccades.aws_r50px,"omitnan"); aws50px_std{end+1} = std(currentSaccades.aws_r50px,"omitnan"); 
            % lum50px_mean{end+1} = mean(currentSaccades.lum_r50px,"omitnan");  lum50px_std{end+1} = std(currentSaccades.lum_r50px,"omitnan"); 
        end
    end

    % (3) Failsafe when there are no ripples for a given subject (happens sometimes) 
    allRipples_across_channels(cellfun(@(x) ~iscell(x) || isempty(x), allRipples_across_channels)) = {{}}; 
    allRipples_across_channels_channames(cellfun(@(x) ~iscell(x) || isempty(x), allRipples_across_channels_channames)) = {{}};
    allRipples_timecourse_across_channels(cellfun(@(x) ~iscell(x) || isempty(x), allRipples_timecourse_across_channels)) = {{}}; 
   
     validIDs = []; %cell(1,80);
     for it = 1:params.trialNum
        %%% (5) Get saccade density 
        %%% make sure that we do not re-count the same ripples from different channels
      %  if ~isempty(vertcat(allRipples_across_channels_channames{:}))

        validIDs= [ validIDs cellfun(@(x) ~any(cellfun(@(y) ~isempty(intersect(x,y)), ...
            allRipples_across_channels{it}(1:(find(cellfun(@(z)isequal(z,x),allRipples_across_channels{it}),1)-1)))),...
            allRipples_across_channels{it})];    
     %   end
        sel = strcmp(saccades_table.subject_new, ripple_data.all_subjNames{sb}) & saccades_table.trial_number == it;
        currentSaccades = saccades_table(sel, :);

        parfor i = 1:length(info_columns)
            saccade_timestamps_all{i,sb,it} = table2array(currentSaccades(:,info_columns(i)))';
        end


        %%% (2) Save ripple features to Saccade Table
        curr_saccadeTiming = currentSaccades.saccade_onset_time*500 + 2000;
        curr_fixationTiming_beg = currentSaccades.saccade_offset_time*500 + 2000;
        curr_fixationTiming_end = ([currentSaccades.saccade_onset_time(2:end); 6])*500 + 2000;

        saccades_table(sel,:).fixation_beg = currentSaccades.saccade_offset_time;
        saccades_table(sel,:).fixation_end = ([currentSaccades.saccade_onset_time(2:end); 6]);
        saccades_table(sel,:).fixation_len = saccades_table(sel,:).fixation_end - saccades_table(sel,:).fixation_beg;
        % jeśli brak rippli w tym trialu -> zapisz bezpieczne wartości
        if isempty(allRipples_across_channels{it}) || all(cellfun(@isempty, allRipples_across_channels{it}))
            saccades_table(sel,:).closestRipple = nan(size(curr_saccadeTiming));  saccades_table(sel,:).closestRipple_TimingDiff = nan(size(curr_saccadeTiming));
            saccades_table(sel,:).closestRipple_Length = zeros(size(curr_saccadeTiming)); saccades_table(sel,:).rippleNumIn200ms = zeros(size(curr_saccadeTiming));
            saccades_table(sel,:).rippleNumIn200ms_after = zeros(size(curr_saccadeTiming)); saccades_table(sel,:).rippleNumIn200ms_before = zeros(size(curr_saccadeTiming));
        else
            %curr_rippleTiming = cellfun(@(c) (isempty(c)||~isnumeric(c))*NaN + (~(isempty(c)||~isnumeric(c)))*median(c), ripCells);
            curr_rippleTiming = cellfun(@median, allRipples_across_channels{it});
            curr_channelNames = string(allRipples_across_channels_channames{it});
            [~, idxNearest] = min(abs(curr_saccadeTiming(:) - curr_rippleTiming(:)'), [], 2);

            saccades_table(sel,:).closestRipple_idx = idxNearest;
            if length(curr_rippleTiming) == 1; idxNearest = idxNearest'; end

            % długości i konwersje (Fs = 500)
            saccades_table(sel,:).closestRipple_Length = cellfun(@numel, allRipples_across_channels{it}(idxNearest))'/500;
            saccades_table(sel,:).closestRipple = (curr_rippleTiming(idxNearest)'-2000)/500;
            saccades_table(sel,:).closestRipple_channel = curr_channelNames(idxNearest)';

            saccades_table(sel,:).closestRipple_TimingDiff = (curr_saccadeTiming - curr_rippleTiming(idxNearest)')/500;

            % Odległość Ripple od losowego przedziału
            curr_random_Timing=saccades_table(sel,:).random_sac_timing_index*500 + 2000;
            [~, idxNearest_random] = min(abs(curr_random_Timing(:) - curr_rippleTiming(:)'), [], 2);
            if length(curr_rippleTiming) == 1; idxNearest_random = idxNearest_random'; end
            saccades_table(sel,:).closestRipple_random_Length = cellfun(@numel, allRipples_across_channels{it}(idxNearest_random))'/500;
            saccades_table(sel,:).closestRipple_random = (curr_rippleTiming(idxNearest_random)'-2000)/500;


            % liczenie ripple w +/-200 ms (thr = 100 próbek)
            delta = curr_rippleTiming(:)' - curr_saccadeTiming(:);
            saccades_table(sel,:).rippleNumIn200ms = sum(abs(delta) <= 100, 2);
            saccades_table(sel,:).rippleNumIn200ms_after = sum(delta >= 0 & delta <= 100, 2);
            saccades_table(sel,:).rippleNumIn200ms_before = sum(delta <= 0 & delta >= -100, 2);

            % Ripples in a closest proximity ( a duration of saccades in 4 consecutive saccade like intervals. With disclaimer on fixation length
            saccades_table(sel,:).rippleNumBefore_saccade = sum(delta <= 0 & delta >= -saccades_table(sel,:).duration/2, 2);
            saccades_table(sel,:).rippleNumIn_saccade     = sum(delta >= 0 & delta <=  saccades_table(sel,:).duration/2, 2);
            saccades_table(sel,:).rippleNumAfter_saccade  = sum(delta <= 2*saccades_table(sel,:).duration/2 & delta >=   saccades_table(sel,:).duration/2 & 2*saccades_table(sel,:).duration/2 < saccades_table(sel,:).fixation_len*500, 2);
            saccades_table(sel,:).rippleNumAfter2_saccade = sum(delta <= 3*saccades_table(sel,:).duration/2 & delta >= 2*saccades_table(sel,:).duration/2 & 3*saccades_table(sel,:).duration/2 < saccades_table(sel,:).fixation_len*500, 2);
            saccades_table(sel,:).rippleNumAfter3_saccade = sum(delta <= 4*saccades_table(sel,:).duration/2 & delta >= 3*saccades_table(sel,:).duration/2 & 4*saccades_table(sel,:).duration/2 < saccades_table(sel,:).fixation_len*500, 2);
            saccades_table(sel,:).rippleNumAfter4_saccade = sum(delta <= 5*saccades_table(sel,:).duration/2 & delta >= 4*saccades_table(sel,:).duration/2 & 5*saccades_table(sel,:).duration/2 < saccades_table(sel,:).fixation_len*500, 2);

            % Ripples in a closest proximity ( a mean duration of saccades in 4 consecutive saccade like intervals. With disclaimer on fixation length
            saccades_table(sel,:).rippleNumBefore_saccade_M = sum(delta <= 0 & delta >= -nanmean( saccades_table.duration/2), 2);
            saccades_table(sel,:).rippleNumIn_saccade_M     = sum(delta >= 0 & delta <=  nanmean( saccades_table.duration/2), 2);
            saccades_table(sel,:).rippleNumAfter_saccade_M  = sum(delta <= 2*nanmean( saccades_table.duration/2) & delta >=   nanmean( saccades_table.duration/2) & 5*nanmean( saccades_table.duration/2) < saccades_table(sel,:).fixation_len*500, 2);
            saccades_table(sel,:).rippleNumAfter2_saccade_M = sum(delta <= 3*nanmean( saccades_table.duration/2) & delta >= 2*nanmean( saccades_table.duration/2) & 5*nanmean( saccades_table.duration/2) < saccades_table(sel,:).fixation_len*500, 2);
            saccades_table(sel,:).rippleNumAfter3_saccade_M = sum(delta <= 4*nanmean( saccades_table.duration/2) & delta >= 3*nanmean( saccades_table.duration/2) & 5*nanmean( saccades_table.duration/2) < saccades_table(sel,:).fixation_len*500, 2);
            saccades_table(sel,:).rippleNumAfter4_saccade_M = sum(delta <= 5*nanmean( saccades_table.duration/2) & delta >= 4*nanmean( saccades_table.duration/2) & 5*nanmean( saccades_table.duration/2) < saccades_table(sel,:).fixation_len*500, 2);

            delta2 =  curr_rippleTiming(:)' - curr_fixationTiming_beg(:);
            saccades_table(sel,:).rippleAroundSaccade = sum(delta > -50 & delta2 < 50,2);
            saccades_table(sel,:).rippleBeforeSaccade = sum(delta < -50 & delta > -100,2);
            saccades_table(sel,:).rippleAfterSaccade = sum(delta2 < 100 & delta2 > 50,2);



            % Get Ripples in a closest proximit
            for i = 1:length(sort(-(13:round(nanmean( saccades_table.duration/2)):501)))
                saccRippTiming(sel,38-(i-1)) = sum(delta <= -round(nanmean( saccades_table.duration/2))*(i-1) & delta >= -round(nanmean( saccades_table.duration/2))*(i),2);
            end
            for i = 1:length(13:round(nanmean( saccades_table.duration/2)):501)
                saccRippTiming(sel,38+i)= sum(delta >= round(nanmean(saccades_table.duration/2))*(i-1) & delta <= round(nanmean(saccades_table.duration/2))*(i), 2);        
            end

            % liczenie ripple podczas fiksacji:
            saccades_table(sel,:).ripples_in_fixation = sum(curr_rippleTiming <= curr_fixationTiming_end & curr_fixationTiming_beg <= curr_rippleTiming,2);
            [~, idxNearest_fix] = min(abs(curr_fixationTiming_beg(:) - curr_rippleTiming(:)'), [], 2);
            if length(curr_rippleTiming) == 1; idxNearest_fix = idxNearest_fix'; end

            saccades_table(sel,:).closest_rip_to_fix_beg = (curr_rippleTiming(idxNearest_fix)'-2000)/500;
            saccades_table(sel,:).closest_rip_to_fix_beg_diff = (curr_fixationTiming_beg - curr_rippleTiming(idxNearest_fix)')/500;


            curr_fixationTiming_half = curr_fixationTiming_beg + saccades_table(sel,:).fixation_len*500;
            [~, idxNearest_fixh] = min(abs(curr_fixationTiming_half(:) - curr_rippleTiming(:)'), [], 2);
            if length(curr_rippleTiming) == 1; idxNearest_fixh = idxNearest_fixh'; end

            saccades_table(sel,:).closest_rip_to_fix_half_diff = (curr_fixationTiming_half - curr_rippleTiming(idxNearest_fixh)')/500;
            saccades_table(sel,:).closest_rip_to_fix_half = (curr_rippleTiming(idxNearest_fixh)'-2000)/500;

            % saccades_table(sel,:).closest_rip_to_fix_end_diff =
            % saccades_table(sel,:).closest_rip_to_fix_end =
        end
     end

     %%% (3) Create a quick, ripple-wise table (for better searches)

    allRipples_timecourse_across_channels_all{sb} = allRipples_timecourse_across_channels;
    allRipples_across_channels_channames_all{sb} = allRipples_across_channels_channames;
    allRipples_across_channels_all{sb} = allRipples_across_channels;

  %  if ~isempty(vertcat(allRipples_across_channels_channames{:})) % Only add when there are even ripples to begin with for a current subject
        aaaa = repelem(1:numel(allRipples_across_channels_channames),  cellfun(@numel, allRipples_across_channels_channames))';
        bbbb = cellfun(@(x) x(end), [allRipples_across_channels{:}]');
        unfolded_channels = [unfolded_channels; string([allRipples_across_channels_channames{:}])']; %Ripple channel unfolding
        unfolded_beg      = [unfolded_beg; cellfun(@(x) x(1), [allRipples_across_channels{:}]')];    %Ripple beginning timestamp unfolding
        unfolded_end      = [unfolded_end; cellfun(@(x) x(end), [allRipples_across_channels{:}]')];  %Ripple Ending  timestamp unfolding
        unfolded_median   = [unfolded_median; cellfun(@median,  [allRipples_across_channels{:}]')];  %Ripple medinan timestamp unfolding % PEAK TBD! 
        unfolded_times    = [unfolded_times; vertcat(allRipples_timecourse_across_channels{:})];     %Ripple Timecourse 
        unfolded_trials   = [unfolded_trials; repelem(1:numel(allRipples_across_channels_channames),  ... % Ripple trial number unfolding
                                            cellfun(@numel, allRipples_across_channels_channames))'];
        unfolded_repeat   = [unfolded_repeat; validIDs'];                                 % Whether ripple is coinciding with previous ripple for easy controll
        tmp = cellfun(@(c) (1:numel(c))', allRipples_across_channels, 'UniformOutput', false); 
        unfolded_rip_num  = [unfolded_rip_num; vertcat(tmp{:}) ];                                  % Current ripple number unfolded
        unfolded_subject  = [string(unfolded_subject); repmat(ripple_data.all_subjNames{sb},length([allRipples_across_channels_channames{:}]),1)]; % Subject unfolding 
        
 %   end

    %%% (5) get saccade Kernel Density per Subject
    for i = 1:length(info_columns)
        [saccade_density_sbj(i,sb,:), ~] = ksdensity(cell2mat(squeeze(saccade_timestamps_all(i,sb,:))')*500 +2000,1:7001, 'Bandwidth', 200); % density with kernel estimataion
    end

    %%% (4) get current Ripple Kernel Density per Subject 
    [ripple_density_sbj(sb,:), ~] = ksdensity(cell2mat(ripple_densities_per_subj(sb,:)), 'Bandwidth', 200); % density with kernel estimataion
end


%%% (4) (5) Densties per tral for rpple and saccade
kVec = 1:numel(info_columns);  
parfor it = 1:length(ripple_timestamps_all)
    %%% get current Ripple Kernel Density per Trial 
    [ripple_density_tri(it,:), ripple_density_timecourse(it,:)] = ksdensity(ripple_timestamps_all{it}, 'Bandwidth', 200); % density with kernel estimataion

    %%% get current Saccade Kernel Density per Trial 
    for k = kVec
        [saccade_density_tri(k,it,:), s_xi(k,it,:)] = ksdensity(cell2mat(saccade_timestamps_all(k,:,it))*500 +2000,1:7001, 'Bandwidth', 200); % density with kernel estimataion
    end
end

x_common_ripple = linspace(min(min(ripple_density_timecourse)), max(max(ripple_density_timecourse)), 100);
%x_common_saccade = linspace(min(min(s_xi(3,:,:))), max(max(s_xi(3,:,:))), 100);




%%% (3) Create a Ripple - wise table, for better computations 
T_text = table(unfolded_subject, unfolded_channels, 'VariableNames', {'subject', 'channel'});
T_numeric = array2table(unfolded_times, 'VariableNames', compose("sample_%d", 1:size(unfolded_times,2)));
T_other = table(unfolded_trials,unfolded_rip_num,unfolded_beg, unfolded_end,unfolded_median,unfolded_repeat,'VariableNames', { 'trial','rip_num','beg', 'end','median','is_repeated'});
ripplewise_table = [T_text, T_other, T_numeric];
save(fullfile(saveDataFolder,"ripple_table.mat"),"ripplewise_table");
 
%%% (1) Create aggregated channel-wise ripple table:
% ripple_table = table(...
%     string(subject(:)), cell2mat(subject_n(:)), cell2mat(trial(:)), ...
%     cell2mat(channel(:)), string(channel_name(:)), cell2mat(ripple_number(:)), ...
%     cell2mat(ripple_length(:)), cell2mat(ripple_number_cl(:)),cell2mat(ripple_length_cl(:)) ,...
%     cell2mat(sacc_number(:)), cell2mat(blink_trial(:)), ...
%     cell2mat(mean_sacc_velocity(:)), cell2mat(mean_sacc_amplitude(:)), ...
%     cell2mat(mean_sacc_latency(:)), cell2mat(mean_sacc_duration(:)), ...
%     cell2mat(X_mean(:)), cell2mat(Y_mean(:)), cell2mat(X_std(:)), cell2mat(Y_std(:)), ...
%     cell2mat(aws50px_mean(:)), cell2mat(aws50px_std(:)), ...
%     cell2mat(lum50px_mean(:)), cell2mat(lum50px_std(:)), ...
%     'VariableNames', {'Subject','SubjectIdx','Trial','ChannelIdx','ChannelName',...
%                       'RippleNumber','RippleLength','RippleNumber_cleaned','RippleLength_cleaned','SaccNumber','BlinkTrial',...
%                       'MeanSaccVelocity','MeanSaccAmplitude','MeanSaccLatency','MeanSaccDuration',...
%                       'X_mean','Y_mean','X_std','Y_std','aws50px_mean','aws50px_std',...
%                       'lum50px_mean','lum50px_std'});

save(fullfile(saveDataFolder,"rippleDensity.mat"),'saccades_table','saccade_timestamps_all','ripple_timestamps_all','x_common_ripple','allRipples_across_channels_all',...
                                                    'ripple_density_tri_x_sb','saccade_density_sbj','saccade_density_tri','ripple_density_tri','ripple_density_sbj','saccRippTiming');




%%% --- VISUALISATION --- %%%


%%% 0 Ripple peri-saccadic timing:


%%
if params.viz
    %%% (1) Densities of Probability of Fixations,Saccades and Ripples in Time
    all_ripples = cell2mat(ripple_timestamps_all);
    
    saccadic = squeeze(saccade_timestamps_all(3,:,:));
    all_saccades = cell2mat(reshape(saccadic,1,size(saccadic,2)*size(saccadic,1)));
    all_saccades = all_saccades(~isnan(all_saccades))*500 +2000;
    
    fixations = squeeze(saccade_timestamps_all(4,:,:));
    all_fixations = cell2mat(reshape(fixations,1,size(fixations,2)*size(fixations,1)));
    all_fixations = all_fixations(~isnan(all_fixations))*500 +2000;
    
    figure("Visible","On","PaperOrientation","portrait","Units","normalized","Position",[0 0 1 1]);
    subplot(3,1,1)
    plot_rug_density(all_fixations,'Overall Ripple Density KDE with Rug Plot',0.05,100);xlim([1,7001])
    
    subplot(3,1,2)
    plot_rug_density(all_saccades,'Overall Saccade Density KDE with Rug Plot',0.05,100);xlim([1,7001])
    
    subplot(3,1,3)
    plot_rug_density(all_ripples,'Overall Fixation Density KDE with Rug Plot',0.05);xlim([1,7001])
    
    exportgraphics(gcf,fullfile(saveFolder,"rippleDensity_Trial.png"))
    
    %%
    %%% (2) Trial SEM Density:
    
    
    figure("Visible","On","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
    
    subplot(3,1,1)
    stdshade( ripple_density_tri, 0.5,[0.8500 0.3250 0.0980],x_common_ripple); % Niebieski
    xticks(1:500:7001); xticklabels([-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10]);xlim([1,7001])
    xlabel('Time (s)'); ylabel('Density');title("Ripples")
    
    subplot(3,1,2)
    stdshade( squeeze(saccade_density_tri(3,:,:)), 0.5,[0.8500 0.3250 0.0980],1:7001); % Niebieski
    xticks(1:500:7001); xticklabels([-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10]);xlim([1,7001])
    xlabel('Time (s)'); ylabel('Density');title("Saccades")
    
    subplot(3,1,3)
    stdshade( squeeze(saccade_density_tri(4,:,:)), 0.5,[0.8500 0.3250 0.0980],1:7001); % Niebieski
    xticks(1:500:7001); xticklabels([-1,0,1,2,3,4,5,6,7]);xlim([1,7001])
    xlabel('Time (s)'); ylabel('Density');title("Fixations")
    
    sgtitle("Mean Density across Trials (n = 19)")
    
    exportgraphics(gcf,fullfile(saveFolder,"rippleDensity_perTrial.png"))
    
    
    %%
    %%% (3) Subject SEM Density:
    
    
    figure("Visible","On","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
    
    subplot(3,1,1)
    stdshade( ripple_density_sbj, 0.5,[0.8500 0.3250 0.0980],x_common_ripple); % Niebieski
    xticks(1:500:7001); xticklabels([-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10]);xlim([1,7001])
    xlabel('Time (s)'); ylabel('Density');title("Ripples")
    
    subplot(3,1,2)
    stdshade( squeeze(saccade_density_sbj(3,:,:)), 0.5,[0.8500 0.3250 0.0980],1:7001); % Niebieski
    xticks(1:500:7001); xticklabels([-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10]);xlim([1,7001])
    xlabel('Time (s)'); ylabel('Density');title("Saccades")
    
    subplot(3,1,3)
    stdshade( squeeze(saccade_density_sbj(4,:,:)), 0.5,[0.8500 0.3250 0.0980],1:7001); % Niebieski
    xticks(1:500:7001); xticklabels([-1,0,1,2,3,4,5,6,7]);xlim([1,7001])
    xlabel('Time (s)'); ylabel('Density');title("Fixations")
    
    sgtitle("Mean Density across Subjects (N = 19)")
    exportgraphics(gcf,fullfile(saveFolder,"rippleDensity_perSubject.png"))
    
    
    
    
    %%
    
    %%% (4) PER  Trial density:
    parfor it = 1:length(ripple_timestamps_all)
        trial_density = squeeze(ripple_density_tri_x_sb(it,:,:));
        
        figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
    
        subplot(3,1,1)
        stdshade( trial_density, 0.5,[0.8500 0.3250 0.0980],1:7001); % Niebieski
        xticks(1:500:7001); xticklabels([-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10])
        xlabel('Time (s)'); ylabel('Density'); title("Ripple")
    
        subplot(3,1,2)
        plot(1:7001,squeeze(saccade_density_tri(3,it,:))'); % Niebieski
        xticks(1:500:7001); xticklabels([-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10])
        xlabel('Time (s)'); ylabel('Density'); title("Saccades")
          
        subplot(3,1,3)
        plot(1:7001,squeeze(saccade_density_tri(4,it,:))'); % Niebieski
        xticks(1:500:7001); xticklabels([-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10])
        xlabel('Time (s)'); ylabel('Density'); title("Fixation")
    
        sgtitle(sprintf("Distribution per Trial: %d (Ripple n. = %d)",it,length(cell2mat(ripple_densities_per_subj(:,it)'))))
    
    
        print(gcf,fullfile(saveFolder,sprintf('tmp_rippleDensity_sbj_%04d.ps', it)),'-dpsc','-append','-fillpage')
        close(gcf)  
    end
    
    % Scal wszystkie tmp .ps w jeden plik
    tmpPaths = fullfile(saveFolder, { dir(fullfile(saveFolder, 'tmp_rippleDensity_sbj_*.ps')).name});
    gsCmd = sprintf('gswin64c -q -dBATCH -dNOPAUSE -sDEVICE=ps2write -sOutputFile="%s" %s', fullfile(saveFolder, 'rippleDensity_Trial_sbj.ps'), strjoin(tmpPaths, ' '));
    status = system(gsCmd); if status ~= 0; error('Ghostscript ps2write failed with status %d', status); end
    
    system(sprintf('gswin64c -sDEVICE=pdfwrite -o "%s" "%s"', strrep(fullfile(saveFolder, 'rippleDensity_Trial_sbj.ps'), '.ps', '.pdf'), fullfile(saveFolder, 'rippleDensity_Trial_sbj.ps')));
    delete(fullfile(saveFolder, 'tmp_rippleDensity_sbj_*.ps')); delete(fullfile(saveFolder, 'rippleDensity_Trial_sbj.ps'));
    
    
    
    %%
    %%%  (5) Plot Per-trial rpple Density 
    % parfor it = 1:length(ripple_timestamps_all)
    %     trial_density = ripple_timestamps_all{it};
    %     [ripple_density_tri, ripple_density_timecourse] = ksdensity(trial_density, 'Bandwidth', 200);
    % 
    %     figure("Visible","Off","PaperOrientation","portrait","Units","normalized","Position",[0 0 1 1]);
    %     ax1 =subplot(3,1,[1 2]);
    %     imshow(imread(fullfile(imageFolder,images{it})),'Parent',ax1);
    % 
    % 
    %     subplot(3,1,3)
    %     plot_rug_density(trial_density,sprintf("Ripple Density KDE & Rug Plot: trial %d/%d", it, length(ripple_timestamps_all)),0)    
    % 
    %     print(gcf,fullfile(saveFolder,sprintf('tmp_rippleDensity_%04d.ps', it)),'-dpsc','-append','-fillpage')
    %     close(gcf)  
    % end
    % 
    % % Scal wszystkie tmp .ps w jeden plik
    % tmpPaths = fullfile(saveFolder, { dir(fullfile(saveFolder, 'tmp_rippleDensity_*.ps')).name});
    % gsCmd = sprintf('gswin64c -q -dBATCH -dNOPAUSE -sDEVICE=ps2write -sOutputFile="%s" %s', fullfile(saveFolder, 'rippleDensity_Trial.ps'), strjoin(tmpPaths, ' '));
    % status = system(gsCmd); %if status ~= 0; %error('Ghostscript ps2write failed with status %d', status); end
    % 
    % system(sprintf('gswin64c -sDEVICE=pdfwrite -o "%s" "%s"', strrep(fullfile(saveFolder, 'rippleDensity_Trial.ps'), '.ps', '.pdf'), fullfile(saveFolder, 'rippleDensity_Trial.ps')));
    % delete(fullfile(saveFolder, 'tmp_rippleDensity_*.ps')); delete(fullfile(saveFolder, 'rippleDensity_Trial.ps'));
    
    
    %%
    
    %%%  (6) Plot Per-subject:
    parfor sb = 1:length(ripple_data.all_subjNames )
        trial_density = squeeze(ripple_density_tri_x_sb(:,sb,:));
        
    
        subplot(3,1,1)
        stdshade( trial_density, 0.5,[0.8500 0.3250 0.0980],1:7001); % Niebieski
        xticks(1:500:7001); xticklabels([-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10])
        xlabel('Time (s)'); ylabel('Density'); title("Ripple")
    
        subplot(3,1,2)
        plot(1:7001,squeeze(saccade_density_sbj(3,sb,:))'); % Niebieski
        xticks(1:500:7001); xticklabels([-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10])
        xlabel('Time (s)'); ylabel('Density'); title("Saccades")
          
        subplot(3,1,3)
        plot(1:7001,squeeze(saccade_density_sbj(4,sb,:))'); % Niebieski
        xticks(1:500:7001); xticklabels([-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10])
        xlabel('Time (s)'); ylabel('Density'); title("Fixation")
    
    
        sgtitle(sprintf("Distribution per Subject: %s (Ripple n. = %d",...
            ripple_data.all_subjNames{sb},length(cell2mat(ripple_densities_per_subj(sb,:)))))
    
        print(gcf,fullfile(saveFolder,sprintf('tmp_rippleDensity_tri_%04d.ps', sb)),'-dpsc','-fillpage')
        close(gcf)  
    end
    
    % Scal wszystkie tmp .ps w jeden plik
    tmpPaths = fullfile(saveFolder, { dir(fullfile(saveFolder, 'tmp_rippleDensity_tri_*.ps')).name});
    gsCmd = sprintf('gswin64c -q -dBATCH -dNOPAUSE -sDEVICE=ps2write -sOutputFile="%s" %s', fullfile(saveFolder, 'rippleDensity_Trial_tri.ps'), strjoin(tmpPaths, ' '));
    status = system(gsCmd); if status ~= 0; error('Ghostscript ps2write failed with status %d', status); end
    
    system(sprintf('gswin64c -sDEVICE=pdfwrite -o "%s" "%s"', strrep(fullfile(saveFolder, 'rippleDensity_Trial_tri.ps'), '.ps', '.pdf'), fullfile(saveFolder, 'rippleDensity_Trial_tri.ps')));
    delete(fullfile(saveFolder, 'tmp_rippleDensity_tri_*.ps')); delete(fullfile(saveFolder, 'rippleDensity_Trial_tri.ps'));
    
    
    %%
    %%%  (7) Plot Per-trial rpple Density 
    parfor sb = 1:length(ripple_data.all_subjNames )
        trial_density = cell2mat(ripple_densities_per_subj(sb,:));
        [ripple_density_tri, ripple_density_timecourse] = ksdensity(trial_density, 'Bandwidth', 200);
    
        figure("Visible","Off","PaperOrientation","portrait","Units","normalized","Position",[0 0 1 1]);
        subplot(2,1,1);
        stdshade( squeeze(ripple_density_tri_x_sb(:,sb,:)), 0.5,[0.8500 0.3250 0.0980],1:7001); % Niebieski
        xticks(1:500:7001); xticklabels([-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10])
        xlabel('Time (s)'); ylabel('Density')
    
        title(sprintf("Distribution per Subject: %s (Ripple n. = %d",...
            ripple_data.all_subjNames{sb},length(cell2mat(ripple_densities_per_subj(sb,:)))))
    
    
        subplot(2,1,2)
        plot_rug_density(trial_density,sprintf("Ripple Density KDE & Rug Plot: trial %d/%d", it, length(ripple_timestamps_all)),0)    
    
        print(gcf,fullfile(saveFolder,sprintf('tmp_rippleDensity_sbjwise_%04d.ps', it)),'-dpsc','-append','-fillpage')
        close(gcf)  
    end
    
    % Scal wszystkie tmp .ps w jeden plik
    tmpPaths = fullfile(saveFolder, { dir(fullfile(saveFolder, 'tmp_rippleDensity_sbjwise_*.ps')).name});
    gsCmd = sprintf('gswin64c -q -dBATCH -dNOPAUSE -sDEVICE=ps2write -sOutputFile="%s" %s', fullfile(saveFolder, 'rippleDensity_Trial_sbjwise.ps'), strjoin(tmpPaths, ' '));
    status = system(gsCmd); if status ~= 0; error('Ghostscript ps2write failed with status %d', status); end
    
    system(sprintf('gswin64c -sDEVICE=pdfwrite -o "%s" "%s"', strrep(fullfile(saveFolder, 'rippleDensity_Trial_sbjwise.ps'), '.ps', '.pdf'), fullfile(saveFolder, 'rippleDensity_Trial_sbjwise.ps')));
    delete(fullfile(saveFolder, 'tmp_rippleDensity_sbjwise_*.ps')); delete(fullfile(saveFolder, 'rippleDensity_Trial_sbjwise.ps'));
end
