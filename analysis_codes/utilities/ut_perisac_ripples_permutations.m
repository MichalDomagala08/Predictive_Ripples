
%%% This is a helper code that generates permutations 

allRipples_across_channels_all = data.allRipples_across_channels_all;


%%% --- PERMUTATIONS COMPUTATION ---


permNum = 1000;
window_len = 7001;

meanRipple = round(nanmean( saccades_table_cleaned.duration/2));
binBoundaries = [sort(-(13:round(nanmean( saccades_table_cleaned.duration/2)):501)) 13:round(nanmean( saccades_table_cleaned.duration/2)):501];

saccRippTiming_perms = zeros(permNum,length(binBoundaries));
saccRippTiming_perms_subwise = zeros(permNum,length(ripple_data.all_subjNames),length(binBoundaries)); % Remember! One subject is outed because it has no observations! 
saccRippTiming_perms_trialwise = zeros(permNum,81,length(binBoundaries));
parfor perm = 1:permNum
    saccRippTiming_temp = zeros(height(saccades_table_cleaned),length(binBoundaries));
    temp_sub = zeros(length(ripple_data.all_subjNames), length(binBoundaries));

    for sb = 1:length(ripple_data.all_subjNames )
        if isempty(allRipples_across_channels_all{sb}); continue; end
        for it = 1:params.trialNum
            if isempty(allRipples_across_channels_all{sb}{it}) || all(cellfun(@isempty, allRipples_across_channels_all{sb}{it})); continue; end
    
            sel = strcmp(saccades_table_cleaned.subject_new, ripple_data.all_subjNames{sb}) & saccades_table_cleaned.trial_number == it;
            curr_saccadeTiming =  saccades_table_cleaned(sel, :).saccade_onset_time*500 + 2000;
    
            curr_rippleTiming_perm = cellfun(@median,  allRipples_across_channels_all{sb}{it});
            rip_shifted = mod(curr_rippleTiming_perm - 1 + randi(window_len-1), 7001) + 1;
            
            delta = rip_shifted(:)' - curr_saccadeTiming(:);
    
            %%% Szybki Przed
            aa = sum(delta <= -round(nanmean( saccades_table_cleaned.duration/2))*(1-1) & delta >= -round(nanmean( saccades_table_cleaned.duration/2))*(1),2);
    
            %%% Szybki Po 
            bb = sum(delta <= -round(nanmean( saccades_table_cleaned.duration/2))*(1-1) & delta >= -round(nanmean( saccades_table_cleaned.duration/2))*(1),2);
    
            % Get Ripples in a closest proximit
            for i = 1:length(sort(-(13:round(nanmean( saccades_table_cleaned.duration/2)):501)))
                saccRippTiming_temp(sel,38-(i-1)) = sum(delta <= -round(nanmean( saccades_table_cleaned.duration/2))*(i-1) & delta >= -round(nanmean( saccades_table_cleaned.duration/2))*(i),2);
            end
            for i = 1:length(13:round(nanmean( saccades_table_cleaned.duration/2)):501)
                saccRippTiming_temp(sel,38+i)= sum(delta >= round(nanmean(saccades_table_cleaned.duration/2))*(i-1) & delta <= round(nanmean(saccades_table_cleaned.duration/2))*(i), 2);        
            end
        end
        temp_sub(sb,:) = sum(saccRippTiming_temp(strcmp(saccades_table_cleaned.subject_new, ripple_data.all_subjNames{sb}), :), 1);   
    end
    saccRippTiming_perms_subwise(perm,:,:) = temp_sub;
    saccRippTiming_perms_trialwise(perm,:,:) = splitapply(@sum,saccRippTiming_temp,findgroups(saccades_table_cleaned.trial_number));
    saccRippTiming_perms(perm,:)= sum(saccRippTiming_temp);
end
saccRippTiming_perms_trialwise = saccRippTiming_perms_trialwise(:,1:end-1,:);


save(fullfile(saveFolder,"Permutation_cl.mat"),"saccRippTiming_perms")

%% 

%%% --- TRUE SACC RIPPLE CLOSENESS ---

saccRippTiming = zeros(height(saccades_table_cleaned),length(binBoundaries));

for sb = 1:length(ripple_data.all_subjNames )
    if isempty(allRipples_across_channels_all{sb}); continue; end

    for it = 1:params.trialNum
        if isempty(allRipples_across_channels_all{sb}{it}) || all(cellfun(@isempty, allRipples_across_channels_all{sb}{it})); continue; end

        sel = strcmp(saccades_table_cleaned.subject_new, ripple_data.all_subjNames{sb}) & saccades_table_cleaned.trial_number == it;
        curr_saccadeTiming =  saccades_table_cleaned(sel, :).saccade_onset_time*500 + 2000;

        curr_rippleTiming_perm = cellfun(@median,  allRipples_across_channels_all{sb}{it});
        rip_shifted = mod(curr_rippleTiming_perm - 1 + randi(window_len-1), 7001) + 1;
        
        delta = curr_rippleTiming_perm(:)' - curr_saccadeTiming(:);

        %%% Szybki Przed
        aa = sum(delta <= -round(nanmean( saccades_table_cleaned.duration/2))*(1-1) & delta >= -round(nanmean( saccades_table_cleaned.duration/2))*(1),2);

        %%% Szybki Po 
        bb = sum(delta <= -round(nanmean( saccades_table_cleaned.duration/2))*(1-1) & delta >= -round(nanmean( saccades_table_cleaned.duration/2))*(1),2);

        % Get Ripples in a closest proximit
        for i = 1:length(sort(-(13:round(nanmean( saccades_table_cleaned.duration/2)):501)))
            saccRippTiming(sel,38-(i-1)) = sum(delta <= -round(nanmean( saccades_table_cleaned.duration/2))*(i-1) & delta >= -round(nanmean( saccades_table_cleaned.duration/2))*(i),2);
        end
        for i = 1:length(13:round(nanmean( saccades_table_cleaned.duration/2)):501)
            saccRippTiming(sel,38+i)= sum(delta >= round(nanmean(saccades_table_cleaned.duration/2))*(i-1) & delta <= round(nanmean(saccades_table_cleaned.duration/2))*(i), 2);        
        end
    end
end


%%
saccRippTiming_perms = load(fullfile(saveFolder,"Permutation_cl.mat")).saccRippTiming_perms;

%%% --- TESTING --- %%%

% (1) A general All file permutations

%%% Compute permutation
true_RipTiming = sum(saccRippTiming);
% deviation of True and Permuted from the mean permutation value (zero distribution)
obs_diff = true_RipTiming - mean(saccRippTiming_perms,1);
perm_diff = bsxfun(@minus, saccRippTiming_perms, mean(saccRippTiming_perms,1));

% Maxsimum deviation from the mean across all bins  
perm_max = max(abs(perm_diff),[],2);

% Famili-wise Error P: Correction checking in how manypermutation urrent
% observation is gretarrt than maximum across bins! 
p_fwer = (1 + sum(perm_max >= abs(obs_diff), 1)) ./ (1 + size(saccRippTiming_perms,1));    
p_trail = double(p_fwer<0.05); p_trail(p_trail ==0) = nan;


%%% Visualise permutation results
figure("Visible","On","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
subplot(2,1,1)
plot (p_fwer); hold on; yline(0.05); title("Permutations P-value")
xticks([1,39,76]); xticklabels(["-1s",'Sac',"1s"])
subplot(2,1,2); hold on; 
stdshade_std(saccRippTiming_perms, 0.5,[0.8500 0.3250 0.0980])
plot(true_RipTiming); xline(38); xline(40);
scatter(1:length(true_RipTiming),p_trail.*true_RipTiming*1.05,[],[0 0 0],'LineWidth',1);
xticks([1,39,76]); xticklabels(["-1s",'Sac',"1s"]); title("Permutations vs True timecourse")
exportgraphics(gcf,fullfile(saveFolder,"perisaccadic_rippleTiming.png"))


% (2) subject-wise mean permutations
    

true_RipTiming_subjet  = splitapply(@sum,saccRippTiming,findgroups(saccades_table_cleaned.subject_new));
obs_diff = mean(true_RipTiming_subjet,1,'omitnan') - squeeze(mean(mean(saccRippTiming_perms_subwise,2,'omitnan'),1,'omitnan'))';
perm_diff = bsxfun(@minus, squeeze(mean(saccRippTiming_perms_subwise,2,'omitnan')),squeeze(mean(mean(saccRippTiming_perms_subwise,2,'omitnan'),1))');
perm_max = max(abs(perm_diff),[],2);
p_fwer = (1 + sum(perm_max >= abs(obs_diff), 1)) ./ (1 + size(saccRippTiming_perms,1));    
p_trail = double(p_fwer<0.05); p_trail(p_trail ==0) = nan;

figure("Visible","On","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
subplot(2,1,1)
plot (p_fwer); hold on; yline(0.05); title("Permutations P-value")
xticks([1,39,76]); xticklabels(["-1s",'Sac',"1s"])
subplot(2,1,2); hold on; 
stdshade_std(squeeze(mean(saccRippTiming_perms_subwise,2,'omitnan')), 0.5,[0.8500 0.3250 0.0980])
%plot(true_RipTiming_subjet); xline(38); xline(40);
stdshade(true_RipTiming_subjet, 0.5,[0.0000, 0.4470, 0.7410]); xline(38); xline(40);

scatter(1:length(mean(true_RipTiming_subjet,1,'omitnan')),p_trail.*mean(true_RipTiming_subjet,1,'omitnan')*1.05,[],[0 0 0],'LineWidth',1);
xticks([1,39,76]); xticklabels(["-1s",'Sac',"1s"]); title("Permutations vs True timecourse")
exportgraphics(gcf,fullfile(saveFolder,"perisaccadic_rippleTiming_subject.png"))

% (3) trial-wise mean permutations


true_RipTiming_trial = mean(splitapply(@sum,saccRippTiming,findgroups(saccades_table_cleaned.trial_number)));
obs_diff = true_RipTiming_trial - squeeze(mean(mean(saccRippTiming_perms_trialwise,2,'omitnan'),1,'omitnan'))';
perm_diff = bsxfun(@minus, squeeze(mean(saccRippTiming_perms_trialwise,2,'omitnan')), squeeze(mean(mean(saccRippTiming_perms_trialwise,2,'omitnan'),1))');
perm_max = max(abs(perm_diff),[],2);
p_fwer = (1 + sum(perm_max >= abs(obs_diff), 1)) ./ (1 + size(saccRippTiming_perms,1));    
p_trail = double(p_fwer<0.05); p_trail(p_trail ==0) = nan;

figure("Visible","On","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
subplot(2,1,1)
plot (p_fwer); hold on; yline(0.05); title("Permutations P-value")
xticks([1,39,76]); xticklabels(["-1s",'Sac',"1s"])
subplot(2,1,2); hold on; 
stdshade_std( squeeze(mean(saccRippTiming_perms_trialwise,2,'omitnan')), 0.5,[0.8500 0.3250 0.0980])
stdshade(splitapply(@sum,saccRippTiming,findgroups(saccades_table_cleaned.trial_number)), 0.5,[0.0000, 0.4470, 0.7410]); xline(38); xline(40);
scatter(1:length(true_RipTiming_trial),p_trail.*true_RipTiming_trial*1.05,[],[0 0 0],'LineWidth',1);
xticks([1,39,76]); xticklabels(["-1s",'Sac',"1s"]); title("Permutations vs True timecourse")
exportgraphics(gcf,fullfile(saveFolder,"perisaccadic_rippleTiming_trial.png"))

