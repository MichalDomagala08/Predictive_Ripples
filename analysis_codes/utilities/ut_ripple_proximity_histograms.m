%%% this script plots ripple precednece and mean ripples in a close
%%% proximity   - WARNING THE HISTOGRAM PROCESSES ARE FAULTY AD WILL ALWAYS
%%% PRODUCE NORMAL-LIKE DISTRIBUTIONS!


figure("Visible","On","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);

% Ripple Saccade precednece histogram (FAULTY AS FUCK XDD
subplot(2,2,1)
x = saccades_table_cleaned.closestRipple_TimingDiff(saccades_table_cleaned.closestRipple_TimingDiff >= -params.bounds & ...
                                                    saccades_table_cleaned.closestRipple_TimingDiff <= params.bounds);         
binWidth = 2/params.nbins; edges = -params.bounds:binWidth:params.bounds;      
histogram(x, 'BinEdges', edges, 'FaceColor', [0.2 0.6 1], 'EdgeColor','none');
xlim([-params.bounds params.bounds]);
xlabel('Timing difference (s)'); ylabel('Count');
title("Saccade- Ripple Timing Difference Distribution")

% 2.2 Plot Ripples closest (100 ms AFTER saccade onset) 
current_saccade_table = saccades_table_cleaned(saccades_table_cleaned.closestRipple_TimingDiff > -0.01 & saccades_table_cleaned.closestRipple_TimingDiff <= 0.05,:);
 
all_sacc_ripples = {};
for sb = 1:length(ripple_data.all_subjNames  )
    current_saccade_table_sb = current_saccade_table(current_saccade_table.subject_new == ripple_data.all_subjNames{sb},:);
    allRipples_across_channels =  cell(1,params.trialNum);
    
    for ch = 1:numel( ripple_data.all_timecourseSWR_trial{sb})
        cur =  ripple_data.all_timecourseSWR_trial{sb}{ch}; if isempty(cur); continue; end
        allRipples_across_channels = cellfun(@(a,b)[a; b], allRipples_across_channels, cur, 'UniformOutput', false);
    end 

    if isempty(allRipples_across_channels);continue; end
    goodTrialNames = unique(current_saccade_table_sb.trial_number);
    
    for tr = 1:length(goodTrialNames)
        curr_final = current_saccade_table_sb(current_saccade_table_sb.trial_number == goodTrialNames(tr),:);
        if isempty(curr_final); continue; end
        all_sacc_ripples{end+1} = allRipples_across_channels{goodTrialNames(tr)}(curr_final.closestRipple_idx,:);
    end
end

subplot(2,2,3)
stdshade( cell2mat(all_sacc_ripples'), 0.5,[0.8500 0.3250 0.0980]); % Niebieski
title("Mean Ripple in first 50ms of Saccades")

% 2.3 SACCADE RIPPLE PRECEDENCE HISTOGRAM 



subplot(2,2,2)
x = saccades_table_cleaned.closest_rip_to_fix_beg_diff(saccades_table_cleaned.closest_rip_to_fix_beg_diff >= -params.bounds & ...
                                                        saccades_table_cleaned.closest_rip_to_fix_beg_diff <= params.bounds);     
binWidth = 2/params.nbins; edges = -params.bounds:binWidth:params.bounds;          
histogram(x, 'BinEdges', edges, 'FaceColor', [0.2 0.6 1], 'EdgeColor','none');
xlim([-params.bounds params.bounds]);
xlabel('Timing difference (s)'); ylabel('Count');
title("Fixation- Ripple Timing Difference Distribution")





current_fixation_table = saccades_table_cleaned(saccades_table_cleaned.closest_rip_to_fix_beg_diff > -0.01 & saccades_table_cleaned.closest_rip_to_fix_beg_diff <= 0.05,:);
 
all_fix_ripples = {};
for sb = 1:length(ripple_data.all_subjNames  )
    current_fixation_table_sb = current_fixation_table(current_fixation_table.subject_new == ripple_data.all_subjNames{sb},:);
    currentRippleSubj = ripple_data.all_timecourseSWR_trial{sb};
    allRipples_across_channels =  cell(1,params.trialNum);
    
    for ch = 1:numel(currentRippleSubj)
        cur = currentRippleSubj{ch}; if isempty(cur); continue; end
        allRipples_across_channels = cellfun(@(a,b)[a; b], allRipples_across_channels, cur, 'UniformOutput', false);
    end 

    if isempty(allRipples_across_channels);continue; end
    goodTrialNames = unique(current_fixation_table_sb.trial_number);
    
    for tr = 1:length(goodTrialNames)
        curr_final = current_fixation_table_sb(current_fixation_table_sb.trial_number == goodTrialNames(tr),:);
        if isempty(curr_final); continue; end
        all_fix_ripples{end+1} = allRipples_across_channels{goodTrialNames(tr)}(curr_final.closestRipple_idx,:);
    end
end

subplot(2,2,4)
stdshade( cell2mat(all_fix_ripples'), 0.5,[0.8500 0.3250 0.0980]); % Niebieski
title("Mean Ripple in first 50ms of Fixation")

if params.checks
    exportgraphics(gcf,fullfile(saveFolder,"occular_vs_ripple_time_difference.png"))
end
