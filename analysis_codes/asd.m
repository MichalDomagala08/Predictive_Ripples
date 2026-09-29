
%%% This script computes representational vectors of a whole brain and
%%% Saccade represenation



%%% PARAMETERS AND SETUP

% General parametes
params.ripple_detection_scheme_name = "RippleDetection_franz_et_al_RipplePeak_manualCorr";
params.analysisName                 = "RippleDensity_Franz_reocurring_ripple";
params.artif_chan_rej               = true; % If we reject artefactual channels
params.blink_saccades_rej           = true; % if we are rejecting blink saccades
params.trialNum = 80; % How many trials are there
params.viz = 0;

% Parameters of Representation
params.box_length = 13;             %length of ms of a boxcar for getting samples: Default 26 ms
params.make_repr = 1;               % to make representaitons actually 
params.nR = 75;


%%% Analysis Parameters
params.ripple_no_coincidence = true;   % Remove from counting ripples that are coinciding with oneanother
params.before_cond =  5;          % N. samples before Sacc Onset that are not coinciding with saccades
params.after_cond  =  5;          % N. samples after Sacc Offset!!!  that are nto coinciding with saccade
params.max_time    = 50;          % maximul n of samples in which we will count ripple (for selected analyses)
params.rip_sac_coincidence = 20;  % Making sure that there is no Saccade - ripple coincidence (as it will drive the Correlation for sure) 
                                          % - If numeric, then use to count  ONLY when distance between saccade and ripple is greater than N samples:
                                          % - If "exact", get exactly when Beg and End of saccades and ripple are divergent

params.ripple_beginning = 500; % if we limit ourselves only to ripples that happened during trials AND NOT DURING first N ms! 

%%% Spike Detection Parameters %%%
params.spikeCtsThresh  = 2;      % Zawsze 1 na czas debugowania!
params.spikePeakWin    = 0.3;   % Zwiększ okno dopasowania (hp) do 150ms
params.spikeZThresh    = 4;    % Obniż próg (z-score na dużych danych rzadko dobija do 4 dla rozlazłych fal)
params.spikeAmpScale   = 2.5;    % Obniż sumaryczny wymóg Peak-to-Trough
params.spikeMNegPeakW  = 150;    % KLUCZ: Pozwól negatywnej fazie trwać do 300ms
params.spikeTrackPeaks = false;   % Szukaj od pozytywnego (tak jak na Twoim obrazku)
params.spikeWindow     = 250;    % Zwiększ margines wycinania artefaktu wokół IED

% Artifact Rejection parameters:
params.artPadding = 70; % how many samples should we pad the artifacts with
params.cluster_tolerance = 20;  % How many samples is it "close" 
params.range_threshold = 6; % How many stds from mean is considered an artifact 
params.iqr_w  = 3; %How many IQR do we have to surpass to have an artefact
params.artif_mad = 8;
analysisName = 'ArtifactRejection_alternative2';



%%% Paths
addpath('D:\Documents_Dell\Predictive_Ripples_2025\analysis_codes\functions'); addpath("utilities\")
dataFolder =  fullfile("D:\Documents_Dell\Predictive_Ripples_2025\data\reref");
saveFolder =  fullfile("D:\Documents_Dell\Predictive_Ripples_2025\Results",params.analysisName); mkdir(saveFolder); 
saveDataFolder = fullfile("D:\Documents_Dell\Predictive_Ripples_2025\preprocessed",params.analysisName); mkdir(saveDataFolder)
imageFolder = "D:\Documents_Dell\Predictive_Ripples_2025\data\stimuli";
saccadeFile = 'D:\Documents_Dell\Predictive_Ripples_2025\data\sacextr_ekm\content_all_subjects.mat';



data = load(fullfile(saveDataFolder,"rippleDensity.mat"));
saccades_table = data.saccades_table;

%%% Load ripple data:
ripple_data = load(fullfile("D:\Documents_Dell\Predictive_Ripples_2025\preprocessed",params.ripple_detection_scheme_name,"rippleData.mat"));

%
ut_bad_channels




%%% --- PREPROCESSING ---


% Concatenate Ripples across channels withholding coincident ripples (for
% both Peaks, Timing and Timecourse

ripple_table = load(fullfile(saveDataFolder,"ripple_table.mat")).ripplewise_table;
ripple_table = ripple_table(logical(ripple_table.is_repeated),:);


subject_level_repr_rp = cell(1,19);
subject_level_repr_sac = cell(1,19);
subject_level_repr_fix = cell(1,19);
files = {dir(fullfile(dataFolder, '*.mat')).name};


%%% CREATE REPRESENTATIONAL DATA
if params.make_repr 
    for sb = 1:length(ripple_data.all_subjNames )
    
        
        fprintf("\n SUBJECT: %s\n\n",ripple_data.all_subjNames{sb})
    
        %%% load and unpack the data structu res
        raw_data = load(fullfile(dataFolder,files{sb}),'data').data.data_eeg;
        currSubjName = ripple_data.all_subjNames{sb};
    
    
        subject_level_repr_rp{sb} = cell(1,params.trialNum);
        subject_level_repr_sac{sb} = cell(1,params.trialNum);
        subject_level_repr_fix{sb} = cell(1,params.trialNum);
        %%% Get Artifact Estimation across all channels - 
        [artifacts_byChan,ieds_byChan,iqr_byChan,range_byChan,data_artifact,ied_timestamps] =  func_artifact_rejection_alternative(raw_data,params);
        %data_artifact = raw_data;
    
    
        for it = 1:length(raw_data.trial)
            
    
            currTrialData = data_artifact.trial{it} - mean(data_artifact.trial{it},2,"omitnan");
            %%% (1) Ripple Locked Representation (Average) 
            currentRipples = table2array(ripple_table(ripple_table.subject ==ripple_data.all_subjNames{sb} & ripple_table.trial == it,[5,6] ));
            subject_level_repr_rp{sb}{it} = zeros(size(data_artifact.trial{it},1),height(currentRipples));
            for rip = 1:size(currentRipples,1)
                 subject_level_repr_rp{sb}{it}(:,rip) = mean(currTrialData(:,currentRipples(rip,1):currentRipples(rip,2)),2,"omitnan");
            end
    
    
            %%% (2) Saccade Loced Representation (Average) 
            curr_sac = saccades_table(strcmp(saccades_table.subject_new,ripple_data.all_subjNames{sb}) & saccades_table.trial_number == it,[22,23]);
            curr_sac  = curr_sac(~isnan(curr_sac.saccade_onset_time) & ~isnan(curr_sac.saccade_offset_time),:);
            subject_level_repr_sac{sb}{it} = zeros(size(data_artifact.trial{it},1),length(height(curr_sac)));
    
            for sac = 1:height(curr_sac)
                currsactime = round((curr_sac(sac,:).saccade_onset_time*500+2000):(curr_sac(sac,:).saccade_offset_time*500+2000));
                subject_level_repr_sac{sb}{it}(:,sac) = mean(currTrialData(:,currsactime),2,"omitnan");
            end


            %%% (3) Fixations Locked Representation (Overlapping windows of 25 ms (or N) 
            curr_fix = saccades_table(strcmp(saccades_table.subject_new,ripple_data.all_subjNames{sb}) & saccades_table.trial_number == it,[48:50]);
            curr_fix  = curr_fix(~isnan(curr_fix.fixation_beg) & ~isnan(curr_fix.fixation_end),:);
            curr_fix  = curr_fix(1:end-1,:);

            subject_level_repr_fix{sb}{it} = zeros(size(data_artifact.trial{it},1),length(height(curr_fix)));

            nWin = ceil(max(curr_fix.fixation_len)*500/params.box_length);
%%
            subject_level_repr_fix{sb}{it} = nan(size(currTrialData,1), height(curr_fix),max(nWin));
            for fix = 1:height(curr_fix)
                fixationSamples =  round((curr_fix(fix,:).fixation_beg*500+2000):(curr_fix(fix,:).fixation_end*500+2000));
                disp(ceil(length(fixationSamples)/params.box_length))

                for win = 1:ceil(length(fixationSamples)/params.box_length)
                    idx = (win-1)*params.box_length+1 : min(win*params.box_length, numel(fixationSamples));
                    subject_level_repr_fix{sb}{it}(:,fix,win) = mean(currTrialData(:,fixationSamples(idx)),2,"omitnan");
                end
            end
%%

    
        end
    
    end

    save(fullfile(saveFolder,"representational_data.mat"),"subject_level_repr_sac","subject_level_repr_rp","subject_level_repr_fix");
else
    repr = load(fullfile(saveFolder,"representational_data.mat"));


end
subject_level_repr_rp = repr.subject_level_repr_rp;
subject_level_repr_sac = repr.subject_level_repr_sac;
%%% Saccade x Ripple Representation Similarity (All against all! for each
%%% trial x Subject:
%%
%%% REPRESENTATIONAL SIMILATRITY: (1) Saccase X Ripple; (2) Ripple x Ripple; (3) Saccade x Saccade
T_all = nan; T_rip = nan; T_sac = nan; T_fix = nan;

%%
for sb = 1:length(ripple_data.all_subjNames )
    for it = 1:params.trialNum

        current_sac_rep = subject_level_repr_sac{sb}{it};
        current_rip_rep = subject_level_repr_rp{sb}{it};
        current_fix_rep = subject_level_repr_fix{sb}{it};

        curr_sac = saccades_table(strcmp(saccades_table.subject_new,ripple_data.all_subjNames{sb}) & saccades_table.trial_number == it,:);
        curr_sac  = table2array(round(curr_sac(~isnan(curr_sac.saccade_onset_time) & ~isnan(curr_sac.saccade_offset_time),[22,23]).*500 + 2000));
        curr_fix = saccades_table(strcmp(saccades_table.subject_new,ripple_data.all_subjNames{sb}) & saccades_table.trial_number == it,:);
        curr_fix  = table2array(round(curr_fix(~isnan(curr_fix.fixation_beg) & ~isnan(curr_fix.fixation_end),[48:50]).*500 + 2000));
        curr_fix = curr_fix(1:end-1,:);
        curr_rip =  table2array(ripple_table(strcmp(ripple_table.subject,ripple_data.all_subjNames{sb}) & ripple_table.trial == it,1:7));

        %%% (1) Saccade vs Ripple Representational Similarity 
        if ~isempty(curr_sac) & ~isempty(curr_rip)

            R = nan(size(current_sac_rep,2), size(current_rip_rep,2));
            for i = 1:size(current_sac_rep,2)
                for j = 1:size(current_rip_rep,2)
                    ok = ~isnan(current_sac_rep(:, i)) & ~isnan(current_rip_rep(:, j));
    
                    if sum(ok) >= 5  % minimalna liczba punktów, ustaw wg. siebie
                        R(i,j) = corr(current_sac_rep(ok, i), current_rip_rep(ok, j), 'type', 'Pearson');
                    else; R(i,j) = NaN; end
                end
            end
    
            [I,J] = ndgrid(1:size(current_sac_rep,2), 1:size(current_rip_rep,2));
            SacIdx = I(:);   RipIdx = J(:);
            SacOn = curr_sac(SacIdx,1); SacOff = curr_sac(SacIdx,2);
            Subject = curr_rip(RipIdx,1);  Trial = curr_rip(RipIdx,3);  Channel = curr_rip(RipIdx,2); 
            
            RipBeg = curr_rip(RipIdx,5); RipMid = curr_rip(RipIdx,7); RipEnd = curr_rip(RipIdx,6);
            Rvec = R(:);
            T = table(Subject, Trial, Channel,SacIdx,RipIdx, SacOn, SacOff, RipBeg, RipEnd,RipMid, Rvec, ...
                'VariableNames',{'subject','trial','channel','sac_idx','rip_idx','sac_onset','sac_offset','rip_beg','rip_end','rip_middle','R'});
            T = T(~isnan(T.R), :); T.rip_beg = str2double( T.rip_beg); T.rip_end = str2double( T.rip_end); T.rip_middle = str2double( T.rip_middle); 
            T.R_fisher = atanh(min(max(T.R, -0.9999), 0.9999));  T.trial = str2double( T.trial);
            T.rip_sac_timing = T.sac_onset - T.rip_beg;
    
            if ~istable(T_all); T_all = T([],:); end
            T_all = [T_all; T];
        end
        
        %%% (2) Ripple vs Ripple Representational Similarity 
        if ~isempty(curr_rip)

            R_rip  = nan(size(current_rip_rep,2), size(current_rip_rep,2));
            for i = 1:size(current_rip_rep,2)-1
                for j = i+1:size(current_rip_rep,2)
                    ok = ~isnan(current_rip_rep(:, i)) & ~isnan(current_rip_rep(:, j));
                    if sum(ok) >= 5  % minimalna liczba punktów, ustaw wg. siebie (Nie posiadających NANów
                        R_rip(i,j) = corr(current_rip_rep(ok, i), current_rip_rep(ok, j), 'type', 'Pearson');
                    else; R_rip(i,j) = NaN; end
    
                end
            end
    
    
            [I,J] = ndgrid(1:size(current_rip_rep,2), 1:size(current_rip_rep,2));
            Rip2Idx = I(:); RipIdx = J(:);  Rvec = R_rip(:);
            Subject = curr_rip(RipIdx,1);  Trial = curr_rip(RipIdx,3);  Channel = curr_rip(RipIdx,2); 
            Rip1Beg =str2double(curr_rip(RipIdx,5)); Rip1Mid = str2double(curr_rip(RipIdx,7)); Rip1End = str2double(curr_rip(RipIdx,6));
            Rip2Beg = str2double(curr_rip(Rip2Idx,5)); Rip2Mid = str2double(curr_rip(Rip2Idx,7)); Rip2End =str2double( curr_rip(Rip2Idx,6));
            
            T = table(Subject, Trial, Channel,RipIdx,Rip2Idx, Rip1Beg, Rip1End, Rip1Mid, Rip2Beg, Rip2End, Rip2Mid, Rvec, ...
                'VariableNames',{'subject','trial','channel','rip_idx_1','rip_idx_2','rip_beg_1','rip_end_1','rip_middle_1','rip_beg_2','rip_end_2','rip_middle_2','R'});
                    T = T(~isnan(T.R),:);
    
            T.R_fisher = atanh(min(max(T.R, -0.9999), 0.9999));  T.trial = str2double( T.trial);
            T.rip_sac_timing = T.rip_beg_1 - T.rip_beg_2;
    
            if ~istable(T_rip); T_rip = T([],:); end
            T_rip = [T_rip; T];
        end

        %%% (3) Saccade vs Saccade Representational Similarity 
        if ~isempty(curr_sac) 
    
            curr_sac_text = saccades_table(strcmp(saccades_table.subject_new,ripple_data.all_subjNames{sb}) & saccades_table.trial_number == it,:);
            curr_sac_text  = table2array(curr_sac_text(~isnan(curr_sac_text.saccade_onset_time) & ~isnan(curr_sac_text.saccade_offset_time),[39,4]));
    
            R_sac = nan(size(current_sac_rep,2), size(current_sac_rep,2));
            for i = 1:size(current_sac_rep,2)-1
                for j = i+1:size(current_sac_rep,2)
                    ok = ~isnan(current_sac_rep(:, i)) & ~isnan(current_sac_rep(:, j));
                    if sum(ok) >= 5  % minimalna liczba punktów, ustaw wg. siebie
                        R_sac(i,j) = corr(current_sac_rep(ok, i), current_sac_rep(ok, j), 'type', 'Pearson');
                    else; R_rip(i,j) = NaN; end
                end
            end
    
    
            [I,J] = ndgrid(1:size(current_sac_rep,2), 1:size(current_sac_rep,2));
            SacIdx1 = I(:); SacIdx2 = J(:);  Rvec = R_sac(:);
            Subject = curr_sac_text(SacIdx1,1);  Trial = curr_sac_text(SacIdx1,2);
            SacOn1 = curr_sac(SacIdx1,1); SacOff1 = curr_sac(SacIdx1,2);
            SacOn2 = curr_sac(SacIdx2,1); SacOff2 = curr_sac(SacIdx2,2);    
    
            T = table(Subject, Trial, SacIdx1,SacIdx2, SacOn1,SacOff1, SacOn2, SacOff2, Rvec, ...
                'VariableNames',{'subject','trial','sac_idx_1','sac_idx_2','sac_onset_1','sac_offset_1','sac_onset_2','sac_offset_2','R'});
            T = T(~isnan(T.R),:);
    
            T.R_fisher = atanh(min(max(T.R, -0.9999), 0.9999));  T.trial = str2double( T.trial);
            T.rip_sac_timing = T.sac_onset_1 - T.sac_onset_2;
    
            if ~istable(T_sac); T_sac = T([],:); end
            T_sac = [T_sac; T];
        end



        %%% (4) RIpple Fixation Repr Similarity
%%
        R = nan(size(current_fix_rep,2), size(current_rip_rep,2),size(current_fix_rep,3));
        for i = 1:size(current_fix_rep,2)
            for j = 1:size(current_rip_rep,2)
                for m = 1:size(current_fix_rep,3)

                    ok = ~isnan(current_fix_rep(:, i,m)) & ~isnan(current_rip_rep(:, j));
    
                    if sum(ok) >= 5  % minimalna liczba punktów, ustaw wg. siebie
                        R(i,j,m) = corr(current_fix_rep(ok, i,m), current_rip_rep(ok, j), 'type', 'Pearson');
                    else; R(i,j,m) = NaN; end
                end
            end
        end

        Rmat = reshape(R, [], size(R,3));     % (fix*rip) rows x windows cols
        Rmat75 = nan(size(Rmat,1), params.nR);
        Rmat75(:,1:min(size(Rmat,2), params.nR)) = Rmat(:,1:min(size(Rmat,2), params.nR));
                
        [M,J] = ndgrid(1:size(current_fix_rep,2), 1:size(current_rip_rep,2));
        FixIdx = M(:); RipIdx = J(:);
        
        FixBeg = curr_fix(FixIdx,1); FixEnd = curr_fix(FixIdx,2);
        Subject = curr_rip(RipIdx,1); Trial = str2double(curr_rip(RipIdx,3)); Channel = curr_rip(RipIdx,2);
        RipBeg = str2double(curr_rip(RipIdx,5)); RipEnd = str2double(curr_rip(RipIdx,6));
        
        Rnames = "R_" + string(1:params.nR);
        T = table(Subject, Trial, Channel, FixIdx, RipIdx, FixBeg, FixEnd, RipBeg, RipEnd, ...
            'VariableNames', {'subject','trial','channel','fix_idx','rip_idx','fix_beg','fix_end','rip_beg','rip_end'});
        T.is_excessive = sum(~isnan(Rmat), 2) > params.nR;
        T.rip_fix_timing_onsets = T.fix_beg - T.rip_end;
        T.rip_fix_timing_middle = (T.fix_beg+ (T.fix_end - T.fix_beg)/2) - (T.rip_end+ (T.rip_end - T.rip_beg)/2);

        T = [T, array2table(Rmat75, 'VariableNames', Rnames)];
        T = T(~any(~isnan(Rmat),2)==0, :);   % drop rows all-NaN across windows
        if ~istable(T_fix); T_fix = T([],:); end
        T_fix = [T_fix; T];


    end
end
%% 
%%


%%% GROUP VARIABLES AND COMPUTE ADDITIONAL INFO

% (1) Mask: Get mask for filtering Coincidental Ripples and fixations: 
segStart = T_fix.fix_beg + (0:params.nR-1)* params.box_length;           % N x 75, start of each segment
segEnd   = min(segStart +  params.box_length - 1, T_fix.fix_end); % N x 75, clipped end
no_overlap_mask = segEnd < T_fix.rip_beg | segStart > T_fix.rip_end;
mask_no_overlap_segments = [true(height(T_fix),width(T_fix) - params.nR), no_overlap_mask]; % N x 86

% (2) Mask:  Whether Ripple happened during that trial
mask_no_overlap_fixations = ~(T_fix.rip_beg <= T_fix.fix_end) & (T_fix.rip_end >= T_fix.fix_beg);

% (3) Mask:  for too short fixations
mask_short_fixations = T_fix.fix_end-T_fix.fix_beg > 5*params.box_length; % at least 5 segments of length)

%%% Cleaned table:
T_fix_cleaned = T_fix(mask_no_overlap_segments);
T_fix_cleaned = T_fix_cleaned(mask_short_fixations & ~T_fix_cleaned.is_excessive,:);

%%% Compute average R fix
T_fix_average = T_fix_cleaned(:,1:12);
T_fix_average.rip_fix_timing_middle_abs = abs(T_fix_average.rip_fix_timing_middle);
T_fix_average.R = table2array(mean(T_fix_cleaned(:,13:86),2,"omitnan"));
T_fix_average.R_abs = abs(T_fix_average.R);
T_fix_average.R_fisher = atanh(min(max(T_fix_average.R, -0.9999), 0.9999)); 
T_fix_average.R_fisher_abs = abs(T_fix_average.R_fisher);

% Get markers of Before/After:
T_fix_average.is_after  = T_fix_average.rip_fix_timing_middle  >= T_fix_average.fix_end - T_fix_average.fix_beg+params.after_cond; % Ripple is After saccade if timing is Positive AND it happened AFTER saccade +/- 5 samples 
T_fix_average.is_before = T_fix_average.rip_fix_timing_middle  <= -params.before_cond;                                  % if the Rip to Sac timing is negative and less than a boudnary of 5 samples
T_fix_average.is_after_restricted = T_fix_average.rip_fix_timing_middle  >= T_fix_average.fix_end - T_fix_average.fix_beg+params.after_cond & ...
                            T_fix_average.rip_fix_timing_middle  <= params.max_time+params.after_cond; % Restricting vounting up to 50ms of ripplesa fter (Seems Fucking reasnable) 
T_fix_average.is_before_restricted= T_fix_average.rip_fix_timing_middle  <= -params.before_cond & ...
                            T_fix_average.rip_fix_timing_middle  >= -params.before_cond-params.max_time;




%%% Make Long Segments: (Gospodi Pomilui) 
Rvars = T_fix_cleaned.Properties.VariableNames(startsWith(T_fix_cleaned.Properties.VariableNames,'R_'));
metaVars = T_fix_cleaned.Properties.VariableNames( ~startsWith(T_fix_cleaned.Properties.VariableNames,'R_'));

rowIdx = repelem((1:height(T_fix_cleaned))',numel(Rvars));
T_fix_tall = T_fix_cleaned(rowIdx,metaVars);
Rmat = T_fix_cleaned{:,Rvars};

T_fix_tall.Segment = repmat((1:numel(Rvars))',height(T_fix_cleaned),1);
T_fix_tall.R       = reshape(Rmat.',[],1);
T_fix_tall.R_abs = abs(T_fix_tall.R);
T_fix_tall.R_fisher = atanh(min(max(T_fix_tall.R, -0.9999), 0.9999)); 
T_fix_tall.R_fisher_abs = abs(T_fix_tall.R_fisher);
T_fix_tall = T_fix_tall(~isnan(T_fix_tall.R), :);

T_fix_tall.seg_beg = T_fix_tall.fix_beg + (T_fix_tall.Segment-1)*params.box_length;
T_fix_tall.seg_end = min(T_fix_tall.seg_beg + params.box_length - 1, T_fix_tall.fix_end);
seg_mid = (T_fix_tall.seg_beg + T_fix_tall.seg_end)/2;
rip_mid = (T_fix_tall.rip_beg + T_fix_tall.rip_end)/2;

T_fix_tall.seg_to_ripple_distance_beg = T_fix_tall.seg_beg - T_fix_tall.rip_beg;
T_fix_tall.seg_to_ripple_distance_mid = seg_mid - rip_mid;

T_fix_tall.is_before = T_fix_tall.seg_end < T_fix_tall.rip_beg - 5;
T_fix_tall.is_after  = T_fix_tall.seg_beg > T_fix_tall.rip_end + 5;
T_fix_tall.is_during = ~T_fix_tall.is_before & ~T_fix_tall.is_after;

T_fix_tall.is_before_restricted = T_fix_tall.is_before & T_fix_tall.rip_beg -T_fix_tall.seg_end <= params.max_time;
T_fix_tall.is_after_restricted  = T_fix_tall.is_after  & T_fix_tall.seg_beg -T_fix_tall.rip_end <= params.max_time;


%%% Ascertain whether the Simialrity is to Early, Middle or Late part of
%%% the ripple (by % of Not Nans

%%% TO DO: 

Rvars = startsWith(T_fix_cleaned.Properties.VariableNames, 'R_');
Rmat  = T_fix_cleaned{:, Rvars};

phaseMat = zeros(size(Rmat)); 
for r = 1:size(Rmat, 1)
    validSeg = find(~isnan(Rmat(r, :)));
    nValid   = numel(validSeg);
    if nValid >= 5
        phaseMat(r, validSeg) = ceil((1:nValid) * 5 / nValid);
    end
end

% Dodanie do T_tall
T_tall.fix_phase = reshape(phaseMat.', [], 1);
T_tall.fix_phase = categorical(T_tall.fix_phase, ...
    0:5, {'missing', 'Early2', 'Early1', 'Middle', 'Late1', 'Late2'});

% Usunięcie braków / fiksacji < 5 segmentów
T_tall = T_tall(~isnan(T_tall.R) & T_tall.fix_phase ~= 'missing', :);

% remove NAN- Rs - they are useless innit? 
T_fix_tall_cleaned = T_fix_tall_cleaned(~isnan(T_fix_tall.R),:);


%% 
%%% ANALYSIS 1 - MEAN Rs:


%%% (1) Are close ripple After vs Before Fixations Exhibiting different degree of Similarity 

% Analysis done on both restricted as well as Far away:

curr_est = []; curr_tStat = []; curr_p = []; nam = []; pred = [];anname = [];

for i = 1:6

    %%% For all Before vs After 
    lme = fitlme(T_fix_average((T_fix_average.is_before | T_fix_average.is_after) & ~overlap,:),sprintf("%s ~ is_after +  rip_sac_timing_abs +(1|subject)", string(gr_names{i})));
    curr_est = [curr_est;lme.Coefficients.Estimate(2:3)]; curr_tStat = [curr_tStat;lme.Coefficients.tStat(2:3)]; curr_p = [curr_p;lme.Coefficients.pValue(2:3)];
    pred = [pred; string(lme.Coefficients.Name(2:3))]; nam = [nam, string(gr_names{i}),string(gr_names{i})]; anname = [anname; "RipProx_Restr"; "RipProx_Restr"];


    figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
    sgtitle(sprintf("All: %s: P: %2.3f T: %2.3f",string(strrep(gr_names{i},'_','-')),lme.Coefficients.pValue(3),lme.Coefficients.tStat(3)))
    subplot(1,2,1)
    boxplot(T_fix_average(T_fix_average.is_before | T_fix_average.is_after & ~overlap,:).(string(gr_names{i})), T_fix_average(T_fix_average.is_before | T_fix_average.is_after  & ~overlap,:).is_after, 'Labels', {'after','before'}); % upewnij się kolejność etykiet
    ylabel(string(strrep(gr_names{i},'_','-'))); title(sprintf("Restricted Before vs After %s",string(strrep(gr_names{i},'_','-'))));
    title("After - Before Boxplot");

    subplot(1,2,2)
    plot_linear_model_bin2(T_fix_average( (~overlap)&(T_fix_average.is_after | T_fix_average.is_before),:),lme, string(gr_names{i}),'is_after',{'rip_sac_timing_abs',"is_after","subject"},1,1,250)
    xticks ([sort([xticks 25])])
    xticklabels(string(str2double(xticklabels) * 2));
    xlabel("Absolute Saccade to Ripple distance (ms)");    
    title("Closeness linear Model: ");

    print(gcf,fullfile(saveFolder,'RSA_ripples_fix_aft_vs_bef.ps'),'-dpsc','-append','-fillpage');


    % For restricted (-/+ params.max ms around Onset)
    lme = fitlme(T_fix_average((T_fix_average.is_before_restricted | T_fix_average.is_after_restricted) & (~overlap),:),sprintf("%s ~ is_after_restricted +  rip_sac_timing_abs +(1|subject) ", string(gr_names{i})));
    curr_est = [curr_est;lme.Coefficients.Estimate(2:3)]; curr_tStat = [curr_tStat;lme.Coefficients.tStat(2:3)]; curr_p = [curr_p;lme.Coefficients.pValue(2:3)];
    pred = [pred; string(lme.Coefficients.Name(2:3))]; nam = [nam, string(gr_names{i}),string(gr_names{i})]; anname = [anname; "RipProx_Restr"; "RipProx_Restr"];


    figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
    sgtitle(sprintf("Restricted (+/- %d ms): %s: P: %2.3f T: %2.3f",params.max_time/2,string(strrep(gr_names{i},'_','-')),lme.Coefficients.pValue(3),lme.Coefficients.tStat(3)))
    subplot(1,2,1)
    boxplot(T_fix_average(T_fix_average.is_before_restricted | T_fix_average.is_after_restricted  & ~overlap,:).(string(gr_names{i})), T_fix_average(T_fix_average.is_before_restricted | T_fix_average.is_after_restricted  & ~overlap,:).is_after_restricted, 'Labels', {'after','before'}); % upewnij się kolejność etykiet
    ylabel(string(strrep(gr_names{i},'_','-'))); title(sprintf("Restricted Before vs After %s",string(strrep(gr_names{i},'_','-'))));
    title("After - Before Boxplot");

    subplot(1,2,2)
    plot_linear_model_bin2(T_fix_average( (~overlap)&(T_fix_average.is_after_restricted | T_fix_average.is_before_restricted),:),lme, string(gr_names{i}),'is_after_restricted',{'rip_sac_timing_abs',"is_after_restricted","subject"},1,1,10)
    xticks ([sort([xticks 25])])
    xticklabels(string(str2double(xticklabels) * 2));
    xlabel("Absolute Saccade to Ripple distance (ms)");    
    title("Closeness linear Model: ");

    print(gcf,fullfile(saveFolder,'RSA_ripples_fix_aft_vs_bef.ps'),'-dpsc','-append','-fillpage');
end


psDir =fullfile(saveFolder,'RSA_ripples_fix_aft_vs_bef.ps');
system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 


T = table(anname,nam',pred,curr_est,curr_tStat,curr_p, ...
'VariableNames', {'anName','Dependent','Predictor','Coef','tStat','pValue'});


fn = fullfile(saveFolder,char(sprintf('LME_RSA_Restr_Fixations_%s_Coinc_%s.xlsx',string(params.max_time),string(params.rip_sac_coincidence)))); sheet = matlab.lang.makeValidName('saccades_aft_vs_bef');
writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
wb.Save(); wb.Close(false); ex.Quit(); delete(ex);




%%% (2) Is Ripple-Fixation similarity particular in specific part of
%%% saccades? 

