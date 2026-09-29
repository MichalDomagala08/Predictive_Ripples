
%%% This is a helper code that generates permutations 

params.perm_no_beg = 1; % consider beggining ripples
params.perm_no_repeated = 1; % consider repeated ripples

saccades_table = saccades_table_cleaned;

allRipples_across_channels_all = data.allRipples_across_channels_all;
ripplewise_table_curr = ripplewise_table;
if params.perm_no_beg;      ripplewise_table_curr = ripplewise_table_curr(ripplewise_table_curr.is_beg,:); end
if params.perm_no_repeated; ripplewise_table_curr = ripplewise_table_curr(logical(ripplewise_table_curr.is_repeated),:); end


permNum = 1000;
window_len = 7001;

meanRipple = round(nanmean( saccades_table.duration/2));
binBoundaries = [sort(-(13:round(nanmean( saccades_table.duration/2)):501)) 13:round(nanmean( saccades_table.duration/2)):501];

saccRippTiming_perms = zeros(permNum,length(binBoundaries));
parfor perm = 1:permNum
    saccRippTiming_temp = zeros(height(saccades_table),length(binBoundaries));

    for sb = 1:length(ripple_data.all_subjNames )
        for it = 1:params.trialNum
            current_Ripples = ripplewise_table_curr(ripplewise_table_curr.subject == ripple_data.all_subjNames{sb} &ripplewise_table_curr.trial ==it,:);
            if ~height(current_Ripples); continue; end % || all(cellfun(@isempty, allRipples_across_channels_all{sb}{it})); continue; end
    
            sel = strcmp(saccades_table.subject_new, ripple_data.all_subjNames{sb}) & saccades_table.trial_number == it;
            curr_saccadeTiming =  saccades_table(sel, :).saccade_onset_time*500 + 2000;
    
            curr_rippleTiming_perm = current_Ripples.median;
            rip_shifted = mod(curr_rippleTiming_perm - 1 + randi(window_len-1), 7001) + 1;
            
            delta = rip_shifted(:)' - curr_saccadeTiming(:);
    
            %%% Szybki Przed
            aa = sum(delta <= -round(nanmean( saccades_table.duration/2))*(1-1) & delta >= -round(nanmean( saccades_table.duration/2))*(1),2);
    
            %%% Szybki Po 
            bb = sum(delta <= -round(nanmean( saccades_table.duration/2))*(1-1) & delta >= -round(nanmean( saccades_table.duration/2))*(1),2);
    
            % Get Ripples in a closest proximit
            for i = 1:length(sort(-(13:round(nanmean( saccades_table.duration/2)):501)))
                saccRippTiming_temp(sel,38-(i-1)) = sum(delta <= -round(nanmean( saccades_table.duration/2))*(i-1) & delta >= -round(nanmean( saccades_table.duration/2))*(i),2);
            end
            for i = 1:length(13:round(nanmean( saccades_table.duration/2)):501)
                saccRippTiming_temp(sel,38+i)= sum(delta >= round(nanmean(saccades_table.duration/2))*(i-1) & delta <= round(nanmean(saccades_table.duration/2))*(i), 2);        
            end
        end
    end
    saccRippTiming_perms(perm,:)= sum(saccRippTiming_temp);
end



%% % True:

saccRippTiming = zeros(height(saccades_table),length(binBoundaries));

for sb = 1:length(ripple_data.all_subjNames )
    for it = 1:params.trialNum
        current_Ripples = ripplewise_table_curr(ripplewise_table_curr.subject == ripple_data.all_subjNames{sb} &ripplewise_table_curr.trial ==it,:);
        if ~height(current_Ripples); continue; end % || all(cellfun(@isempty, allRipples_across_channels_all{sb}{it})); continue; end

        %if isempty(allRipples_across_channels_all{sb}{it}) || all(cellfun(@isempty, allRipples_across_channels_all{sb}{it})); continue; end

        sel = strcmp(saccades_table.subject_new, ripple_data.all_subjNames{sb}) & saccades_table.trial_number == it;
        curr_saccadeTiming =  saccades_table(sel, :).saccade_onset_time*500 + 2000;

        curr_rippleTiming_perm = current_Ripples.median;%= cellfun(@median,  allRipples_across_channels_all{sb}{it});
        
        delta = curr_rippleTiming_perm(:)' - curr_saccadeTiming(:);

        %%% Szybki Przed
        aa = sum(delta <= -round(nanmean( saccades_table.duration/2))*(1-1) & delta >= -round(nanmean( saccades_table.duration/2))*(1),2);

        %%% Szybki Po 
        bb = sum(delta <= -round(nanmean( saccades_table.duration/2))*(1-1) & delta >= -round(nanmean( saccades_table.duration/2))*(1),2);

        % Get Ripples in a closest proximit
        for i = 1:length(sort(-(13:round(nanmean( saccades_table.duration/2)):501)))
            saccRippTiming(sel,38-(i-1)) = sum(delta <= -round(nanmean( saccades_table.duration/2))*(i-1) & delta >= -round(nanmean( saccades_table.duration/2))*(i),2);
        end
        for i = 1:length(13:round(nanmean( saccades_table.duration/2)):501)
            saccRippTiming(sel,38+i)= sum(delta >= round(nanmean(saccades_table.duration/2))*(i-1) & delta <= round(nanmean(saccades_table.duration/2))*(i), 2);        
        end
    end
end
%%
save(fullfile(saveFolder,"Permutation_cl.mat"),"saccRippTiming_perms")


saccRippTiming_perms = load(fullfile(saveFolder,"Permutation_cl.mat"),"saccRippTiming_perms").saccRippTiming_perms;

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

%%
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
