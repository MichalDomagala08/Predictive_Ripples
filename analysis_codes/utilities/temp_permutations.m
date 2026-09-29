
%%% This is a helper code that generates permutations 
permNum = 1000;
window_len = 7001;

meanRipple = round(nanmean( saccades_table.duration/2));
binBoundaries = [sort(-(13:round(nanmean( saccades_table.duration/2)):501)) 13:round(nanmean( saccades_table.duration/2)):501];

saccRippTiming_perms = zeros(permNum,length(binBoundaries));
parfor perm = 1:permNum
    saccRippTiming_temp = zeros(height(saccades_table),length(binBoundaries));

    for sb = 1:length(ripple_data.all_subjNames )
        for it = 1:params.trialNum
            if isempty(allRipples_across_channels_all{sb}{it}) || all(cellfun(@isempty, allRipples_across_channels_all{sb}{it})); continue; end
    
            sel = strcmp(saccades_table.subject_new, ripple_data.all_subjNames{sb}) & saccades_table.trial_number == it;
            curr_saccadeTiming =  saccades_table(sel, :).saccade_onset_time*500 + 2000;
    
            curr_rippleTiming_perm = cellfun(@median,  allRipples_across_channels_all{sb}{it});
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

save(fullfile(saveFolder,"Permutation.mat"),"saccRippTiming_perms")



%%% Compute permutation
P = size(saccRippTiming_perms,1);
perm_mean = mean(saccRippTiming_perms,1);
perm_std  = std(saccRippTiming_perms,[],1); perm_std(perm_std==0)=eps;
obs_diff = true_RipTiming - perm_mean;
perm_diff = bsxfun(@minus, saccRippTiming_perms, perm_mean);
z_map = obs_diff ./ perm_std;                                   % test statistic
p_bin = (1 + sum(abs(perm_diff) >= abs(obs_diff), 1)) ./ (1 + P); % empirical two-sided p per bin
perm_max = max(abs(perm_diff),[],2);
p_fwer = (1 + sum(perm_max >= abs(obs_diff), 1)) ./ (1 + P);    % FWER (max-stat) corrected p per bin
