%%% This script is used to make ripple statistics and their visualisations

%%% --- PARAMETERS  AND PREPARATION--- %%%

% Naming parameter and stimuli sizes
params.artifact_rejection_scheme_name =  'ArtifactRejection_alternative2';
params.ripple_detection_scheme_name   = "RippleDetection_franz_et_al_RipplePeak_manualCorr";
params.analysisName                   = "RippleDensity_Franz";
params.trialNum = 80;


% Analysis Progresion:
params.trialwise_ripple_count  = 1;
params.trialwise_ripple_length = 1;
params.saccade_ripple_timing   = 1;
params.fixation_ripple_timing  = 1;

% General analysis parameter
params.ripples_later = true ;% only get ripples that are not in the first 500ms of the trial beginning to account for increase in energy fallacy
params.trialwise_ripples2 = true; % Whether we do things trialwise_ripples2, or with inclusion ofchannels 
params.checks = true;    % Wrtie model checks 

%%% Visualisation parameters

%histogram vis parameters
params.bounds = 1;         
params.nbins = 100;                  



%%% Paths
addpath('D:\Documents_Dell\Predictive_Ripples_2025\analysis_codes\functions'); addpath("utilities\")
dataFolder =  fullfile("D:\Documents_Dell\Predictive_Ripples_2025\preprocessed",params.analysisName);
saveFolder =  fullfile("D:\Documents_Dell\Predictive_Ripples_2025\Results",params.analysisName);
mkdir(fullfile(saveFolder,"Results"));

%%% load data

data = load(fullfile(dataFolder,"rippleDensity.mat"));
ripple_table3 = data.ripple_table;
saccades_table = data.saccades_table;

ripple_table3 = movevars( ripple_table3, {'Subject','Trial','ChannelName','RippleNumber','RippleNumber_cleaned','RippleLength','RippleLength_cleaned'}, 'Before', 1);
ripple_table3 =  sortrows(ripple_table3, {'Subject','ChannelName','Trial'}, {'ascend','ascend','ascend'});
%%% Load ripple data:
ripple_data = load(fullfile("D:\Documents_Dell\Predictive_Ripples_2025\preprocessed",params.ripple_detection_scheme_name,"rippleData.mat"));


%%% Load Image Entropy data:
entropy_dat = load(fullfile("D:\Documents_Dell\Predictive_Ripples_2025\data\entropy_results.mat")).entropyTable;
aa = outerjoin(saccades_table,entropy_dat,'LeftKeys','name','RightKeys','image');
entropy_trials = unique(aa(:,[1,39,4,72:77]));
entropy_trials = sortrows(entropy_trials, {'subject_new','trial_number'}, {'ascend','ascend'});
 
numCols = [4,5,6,7,8,9];
g = findgroups(entropy_trials.subject_new);
newMat = nan(height(entropy_trials), numel(numCols));
for gg = 1:max(g)
    idx = find(g==gg);
    for k = 1:length(numCols)
        v = entropy_trials{idx, numCols(k)};
        newMat(idx,k) = [nan; diff(v)];           % current - previous
    end
end
for k = 1:numel(numCols)
    entropy_trials.([entropy_trials.Properties.VariableNames{numCols(k)} '_diff']) = newMat(:,k);
end

%%% PREPROCESSING

aaaa = ripple_table(~isnan(ripple_table.RippleLength),:);
ripple_table3 = aaaa ;
%%% Create Trialwise table:
if params.trialwise_ripples
    [~, ia, ~] = unique(ripple_table3(:, {'Subject', 'Trial'}), 'rows', 'first');
    
    % Wybierz wszystkie kolumny oprócz ChannelIdx, ChannelName, RippleNumber, RippleLength
    excludeCols = {'ChannelIdx', 'ChannelName', 'RippleNumber', 'RippleLength','RippleLength_cleaned','RippleNumber_cleaned'};
    includeCols = setdiff(ripple_table3.Properties.VariableNames, excludeCols);
    
    % Stwórz tabelę Out z pierwszych wierszy dla każdej grupy
    trialwise_ripples2 = ripple_table3(ia, includeCols);
    
    % Teraz agreguj tylko te 2 kolumny, które się zmieniają per kanał:
    G = findgroups(ripple_table3.Subject, ripple_table3.Trial);
    trialwise_ripples2.MeanRippleNumberr_cleaned = splitapply(@nanmean, ripple_table3.RippleNumber_cleaned, G);
    trialwise_ripples2.MeanRippleLength_cleaned  = splitapply(@nanmean, ripple_table3.RippleLength_cleaned, G);
    
    trialwise_ripples2.MeanRippleNumber         = splitapply(@nanmean, ripple_table3.RippleNumber, G);
    trialwise_ripples2.MeanRippleLength          = splitapply(@nanmean, ripple_table3.RippleLength, G);
 
    % Liczba kanałów
    trialwise_ripples2.NumChannels = splitapply(@numel, ripple_table3.RippleNumber, G);
    
    % Changing columns to make it
    preferredOrder = {'Subject', 'SubjectIdx', 'Trial', 'NumChannels', 'MeanRippleNumber', 'MeanRippleLength'};
    remainingCols = setdiff(trialwise_ripples2.Properties.VariableNames, preferredOrder);
    newOrder = [preferredOrder, remainingCols];
    trialwise_ripples2 = trialwise_ripples2(:, newOrder);
    trialwise_ripples2 = movevars( trialwise_ripples2, {'Subject','Trial','NumChannels','MeanRippleNumber','MeanRippleNumberr_cleaned','MeanRippleLength','MeanRippleLength_cleaned'}, 'Before', 1);
trialwise_ripples
end

% collapse ripple table and concatenate across channels: 
trialwise_ripples2 = outerjoin(trialwise_ripples2,entropy_trials,'LeftKeys',{'Subject','Trial'},'RightKeys',{'subject_new','trial_number'});
trialwise_ripples2 = trialwise_ripples2(~isnan(trialwise_ripples2.SubjectIdx),:);
%%% Global Filtering:


if params.ripples_later
    additional_mask =  saccades_table.closestRipple > 1 & saccades_table.saccade_onset_time >1;
else
    additional_mask = ones(height(saccades_table));
end

saccades_table_cleaned = saccades_table(~isnan(saccades_table.latency) & ...
                                            ~isnan(saccades_table.blinkSac) & ...
                                            saccades_table.closestRipple ~=0 & ...
                                            ~isnan(saccades_table.closestRipple) & ...
                                            additional_mask,:);


%%% Adding additonal relevant features
saccades_table_cleaned.closestRipple_TimingDiffAbs = abs(saccades_table_cleaned.closestRipple_TimingDiff);
saccades_table_cleaned.closest_rip_to_fix_beg_diffAbs = abs(saccades_table_cleaned.closest_rip_to_fix_beg_diff);
saccades_table_cleaned.closest_rip_to_fix_half_diffAbs = abs(saccades_table_cleaned.closest_rip_to_fix_half_diff);


%%%% ANALYSIS (1) Saccadic Featue and Ripple Numbers
saccadic_feat = [26];

%% %%% (1.1) Ripple Number

if params.trialwise_ripple_count
    filename = 'lme_outputs_trial.txt';
    fid = fopen(filename,'w+');  
    fprintf(fid, '--- %s ---\n\n',datestr(now,'yyyy-mm-dd HH:MM:SS'));
    for i = saccadic_feat
    
        %%% Sub-matrix for currenr predictions
        if params.trialwise_ripples
            if params.ripples_later; currentPredictor = "MeanRippleNumberr_cleaned"; else; currentPredictor = "MeanRippleNumber"; end
            currentDependant = trialwise_ripples2.Properties.VariableNames{i};    
            current_data = trialwise_ripples2(:, ["Subject","SubjectIdx","Trial","NumChannels",currentPredictor,"BlinkTrial",currentDependant]);
            varNames = {currentPredictor,"Trial","NumChannels","Subject"};
        else
            if params.ripples_later; currentPredictor = "RippleNumber_cleaned"; else; currentPredictor = "RippleNumber"; end
            currentDependant = ripple_table3.Properties.VariableNames{i};  
            current_data = ripple_table3(:, ["Subject","SubjectIdx","Trial","ChannelName",currentPredictor,"BlinkTrial",currentDependant]);
            varNames = {currentPredictor,"Trial","ChannelName","Subject"};
        end
    
        % Outlier removing (by robust MAD z-score  (4) 
        robz_x = 0.6745 * ( current_data.(currentPredictor) - median(current_data.(currentPredictor),"omitnan")) ./ mad( current_data.(currentPredictor), 1);
        current_data = current_data(abs(robz_x) <= 4, :); % Próg około 4 robust-σ od mediany
        
    
        %%% Model fiting
        formula = sprintf('%s ~ %s +Trial +%s+ (1|Subject)',currentDependant,currentPredictor,varNames{3}); % defining model structure
        lme = fitlme(current_data, formula); % fit LME
        txt = evalc('disp(lme)');                            
        fprintf(fid, '\n %s ~ %s \n -----------------------------\n\n %s\n\n',currentDependant ,currentPredictor, txt);
        
        
        %%% Plots: Linear Model and its checks 
        figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
        plot_linear_model(current_data,lme,currentDependant,currentPredictor,varNames,1,0.5)
        sgtitle(sprintf("LME_%s_vs_%s",currentDependant,currentPredictor))
        print(gcf,fullfile(saveFolder,"Results","LME_Ripple_Saccadic.ps"),'-dpsc','-append','-fillpage')
        %exportgraphics(gcf,fullfile(saveFolder,"Results",sprintf("LME_%s_vs_%s.png",currentDependant,currentPredictor)))
        close(gcf)
    
        if params.checks
            figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
            plot_diagnostic_lm(current_data, lme, currentDependant,currentPredictor)
            sgtitle(sprintf("LME_%s_vs_%s",currentDependant,currentPredictor))
            print(gcf,fullfile(saveFolder,"Results","LME_Diagnostic_Ripple_Saccadic.ps"),'-dpsc','-append','-fillpage')
            close(gcf)
        end
    
    end
    psDir = fullfile(saveFolder,"Results","LME_Ripple_Saccadic.ps");
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
    
    psDir = fullfile(saveFolder,"Results","LME_Diagnostic_Ripple_Saccadic.ps");
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
end

%% %%% (1.2) - Ripple Length


if params.trialwise_ripple_length
    
    filename = 'lme_outputs_trial_length.txt';
    fid = fopen(filename,'w+');  
    fprintf(fid, '--- %s ---\n\n',datestr(now,'yyyy-mm-dd HH:MM:SS'));
    for i = saccadic_feat
    
        %%% Sub-matrix for currenr predictions
        if params.trialwise_ripples2
            if params.ripples_later; currentPredictor = "MeanRippleLength_cleaned"; else; currentPredictor = "MeanRippleLength"; end
            currentDependant = trialwise_ripples2.Properties.VariableNames{i};    
            current_data = trialwise_ripples2(:, ["Subject","SubjectIdx","Trial","NumChannels",currentPredictor,"BlinkTrial",currentDependant]);
            varNames = {currentPredictor,"Trial","NumChannels","Subject"};
        else
            if params.ripples_later; currentPredictor = "RippleLength_cleaned"; else; currentPredictor = "RippleLength"; end
            currentDependant = ripple_table3.Properties.VariableNames{i};  
            current_data = ripple_table3(:, ["Subject","SubjectIdx","Trial","ChannelName",currentPredictor,"BlinkTrial",currentDependant]);
            varNames = {currentPredictor,"Trial","ChannelName","Subject"};
        end
    
        % Outlier removing (by robust MAD z-score  (4) 
        robz_x = 0.6745 * ( current_data.(currentPredictor) - median(current_data.(currentPredictor),'omitnan')) ./ mad( current_data.(currentPredictor), 1);
        current_data = current_data(abs(robz_x) <= 4, :); % Próg około 4 robust-σ od mediany
        
    
        %%% Model fiting
        formula = sprintf('%s ~ %s +Trial +%s+ (1|Subject)',currentDependant,currentPredictor,varNames{3}); % defining model structure
        lme = fitlme(current_data, formula); % fit LME
        txt = evalc('disp(lme)');                            
        fprintf(fid, '\n %s ~ %s \n -----------------------------\n\n %s\n\n',currentDependant ,currentPredictor, txt);
    
        %%% Plots: Linear Model and its checks 
        figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
        plot_linear_model(current_data,lme,currentDependant,currentPredictor,varNames,1)

        sgtitle(sprintf("LME_%s_vs_%s",currentDependant,currentPredictor))
        print(gcf,fullfile(saveFolder,"Results","LME_Ripple_Length_Saccadic.ps"),'-dpsc','-append','-fillpage')
        %exportgraphics(gcf,fullfile(saveFolder,"Results",sprintf("LME_%s_vs_%s.png",currentDependant,currentPredictor)))
        close(gcf)
    
        if params.checks
            figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
            plot_diagnostic_lm(current_data, lme, currentDependant,currentPredictor)
            sgtitle(sprintf("LME_%s_vs_%s",currentDependant,currentPredictor))
            print(gcf,fullfile(saveFolder,"Results","LME_Diagnostic_Ripple_Length_Saccadic.ps"),'-dpsc','-append','-fillpage')
            close(gcf)
        end
    
    end
    psDir = fullfile(saveFolder,"Results","LME_Ripple_Length_Saccadic.ps");
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
    
    psDir = fullfile(saveFolder,"Results","LME_Diagnostic_Ripple_Length_Saccadic.ps");
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
end
%% (2) ANALYSIS RIPPLES PRECEDENCE 

%% (2.1) SACCADE and FIXATION RIPPLE PRECEDENCE HISTOGRAMS

% This is specifically to be sure that we are not observing Saccadic Spike
%ut_proximity_histograms

%% (2.2)  PERI-SACCADIC COUNT (-1 to +1) 


figure("Visible","On","PaperOrientation","portrait","Units","normalized","Position",[0 0 1 1]);

cleaned_saccRippTiming =data.saccRippTiming(~isnan(saccades_table.latency) & ~isnan(saccades_table.blinkSac) & ...
                                            saccades_table.closestRipple ~=0 & ~isnan(saccades_table.closestRipple) & additional_mask,:);
binBoundaries = [sort(-(13:round(nanmean(saccades_table.duration/2)):501)) 13:round(nanmean( saccades_table.duration/2)):501];

subplot(3,1,3)
[G, ~] = findgroups(saccades_table_cleaned.trial_number);   % G: grupa dla każdego wiersza
perTrialCounts = splitapply(@(I) sum(cleaned_saccRippTiming(I,:),1), (1:size(cleaned_saccRippTiming,1))', G);
counts = double(splitapply(@numel,(1:height(saccades_table_cleaned))',findgroups(saccades_table_cleaned.trial_number)));
stdshade(perTrialCounts./counts, 0.5,[0.8500 0.3250 0.0980]); xline(38); xline(40)
xticks([1,39,length(binBoundaries)]); xticklabels({"-1s", "Sac" , "1s"})
title("TrialWise Ripple near saccade distribution")

subplot(3,1,2)
[G, ~] = findgroups(saccades_table_cleaned.subject_new);   % G: grupa dla każdego wiersza
perSubjectCounts = splitapply(@(I) sum(cleaned_saccRippTiming(I,:),1), (1:size(cleaned_saccRippTiming,1))', G);
counts = double(splitapply(@numel,(1:height(saccades_table_cleaned))',findgroups(saccades_table_cleaned.subject_new)));
stdshade(perSubjectCounts./counts, 0.5,[0.8500 0.3250 0.0980]); xline(38); xline(40)
xticks([1,39,length(binBoundaries)]); xticklabels({"-1s", "Sac" , "1s"})
title("SubjectWise Ripple near saccade distribution")

subplot(3,1,1)
plot(mean(cleaned_saccRippTiming,1), "Color",[0.8500 0.3250 0.0980]); xline(38); xline(40)
xticks([1,39,length(binBoundaries)]); xticklabels({"-1s", "Sac" , "1s"})
title("All Ripple near saccade distribution")

%% (2.3) Barplot of close ripple timings

figure();
subplot(1,2,1)
bar(table2array(sum(saccades_table_cleaned(:,[58:63]),1)))
title("Ripple numbers in mean saccade Distance")
xticks([1,2,3,4,5,6])
xticklabels(["before",'Sac','Fix1','Fix2','Fix3','Fix4'])

subplot(1,2,2)
bar(table2array(sum(saccades_table_cleaned(:,[64:69]),1)))
title("For Mean Saccade Distance")
xticks([1,2,3,4,5,6])
xticklabels(["before",'Sac','Fix1','Fix2','Fix3','Fix4'])

%% (2.4) Median Split alongsie Saccadic Features 

% The feature include Amplitude, peakVelocity, Duration and all! 

%%% Warning! The Timing difference is not a good indicator because it will
%%% always have a concentraned near zero distribution!

% However, using COUNT! that changes things! 

mspl_sacc = saccades_table_cleaned;

for i = [6,7,20,21,50]
    currentDependent = saccades_table_cleaned.Properties.VariableNames{i};
    mspl_sacc.(currentDependent+"_mspl") = saccades_table_cleaned.(currentDependent) >= median(saccades_table_cleaned.(currentDependent));
    plot_barplot_sacc_timing_med_split(mspl_sacc,saveFolder,currentDependent,"Lower_Upper_Rip_Sac_bar")
    plot_histogram_sacc_timing_msplit(mspl_sacc,saveFolder,currentDependent,"closestRipple_TimingDiff","Lower_Upper_Rip_Sac_hist")
    plot_histogram_sacc_timing_msplit(mspl_sacc,saveFolder,currentDependent,"closest_rip_to_fix_half_diff","Lower_Upper_Rip_Fix_hist")
end
   
psDir = fullfile(saveFolder,"Results","Lower_Upper_Rip_Sac_bar.ps");
system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 

psDir = fullfile(saveFolder,"Results","Lower_Upper_Rip_Fix_hist.ps");
system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 

psDir = fullfile(saveFolder,"Results","Lower_Upper_Rip_Sac_hist.ps");
system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 


% %%
% %%% barplots of values
% 
% %%
% saccades_table_cleaned.random_ripple_timing = saccades_table_cleaned.random_sac_timing_index -saccades_table_cleaned.closestRipple_random;
% 
% 
% figure()
% histogram(saccades_table_cleaned.random_ripple_timing)
%%

%%%% ANALYSIS (3) Saccade Entrained Features

,
saccadic_features = [6,7,20,21,12,13,22,23,27,28,29,33];
ripple_features = [41,42,54,71,73];
binwidth = [1,0,1,1,1,0,0,0,0,0];
filename = 'lme_outputs_saccade.txt';
fid = fopen(filename,'w+');  
fprintf(fid, '--- %s ---\n\n',datestr(now,'yyyy-mm-dd HH:MM:SS'));

for j = ripple_features
    for i = saccadic_features
        %%% Sub-matrix for currenr predictions
        currentPredictor = saccades_table_cleaned.Properties.VariableNames{j};
        currentDependant = saccades_table_cleaned.Properties.VariableNames{i};    

        sprintf("Current Analysis: %s ~%s",currentDependant,currentPredictor);
        current_data = saccades_table_cleaned(:, ["subject_new","trial_number",currentPredictor,"blinkSac",currentDependant]);
        varNames = {currentPredictor,"trial_number","blinkSac","subject_new"};
    
    
        % Outlier removing (by robust MAD z-score  (4) 
        robz_x = 0.6745 * ( current_data.(currentPredictor) - median(current_data.(currentPredictor))) ./ mad( current_data.(currentPredictor), 1);
         current_data = current_data(abs(robz_x) <= 4, :); % Próg około 4 robust-σ od mediany
       
        
        %%% Model fiting
        formula = sprintf('%s ~ %s +%s+ (1|%s)',currentDependant,currentPredictor,varNames{2},varNames{4}); % defining model structure
        lme = fitlme(current_data, formula); % fit LME
        txt = evalc('disp(lme)');                            
        fprintf(fid, '\n %s ~ %s \n -----------------------------\n\n %s\n\n',currentDependant ,currentPredictor, txt);
       
        %%% Plots: Linear Model and its checks 

        figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
        plot_linear_model(current_data,lme,currentDependant,currentPredictor,varNames,1,binwidth)
    
        sgtitle(sprintf("LME_%s_vs_%s",currentDependant,currentPredictor))
        print(gcf,fullfile(saveFolder,"Results",sprintf("LME_Saccade_Ripple_%s.ps",currentPredictor)),'-dpsc','-append','-fillpage')
        %exportgraphics(gcf,fullfile(saveFolder,"Results",sprintf("LME_%s_vs_%s.png",currentDependant,currentPredictor)))
        close(gcf)
    
        if params.checks
            figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
            plot_diagnostic_lm(current_data, lme, currentDependant,currentPredictor)
            sgtitle(sprintf("LME_%s_vs_%s",currentDependant,currentPredictor))
            print(gcf,fullfile(saveFolder,"Results",sprintf("LME_Diagnostics_Saccade_Ripple_%s.ps",currentPredictor)),'-dpsc','-append','-fillpage')
            close(gcf)
        end
    
    end
    psDir = fullfile(saveFolder,"Results",sprintf("LME_Saccade_Ripple_%s.ps",currentPredictor));
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
    
    psDir = fullfile(saveFolder,"Results",sprintf("LME_Diagnostics_Saccade_Ripple_%s.ps",currentPredictor));
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 

end


%%

for j = [58:60,64:66]
    for i = saccadic_features
        %%% Sub-matrix for currenr predictions
        currentPredictor = saccades_table_cleaned.Properties.VariableNames{j};
        currentDependant = saccades_table_cleaned.Properties.VariableNames{i};    

        sprintf("Current Analysis: %s ~%s",currentDependant,currentPredictor);
        current_data = saccades_table_cleaned(:, ["subject_new","trial_number",currentPredictor,"blinkSac",currentDependant]);
        varNames = {currentPredictor,"trial_number","blinkSac","subject_new"};
    
        % 
        % % Outlier removing (by robust MAD z-score  (4) 
        % robz_x = 0.6745 * ( current_data.(currentPredictor) - median(current_data.(currentPredictor))) ./ mad( current_data.(currentPredictor), 1);
        %  current_data = current_data(abs(robz_x) <= 4, :); % Próg około 4 robust-σ od mediany
       
        
        %%% Model fiting
        formula = sprintf('%s ~ %s +%s+ (1|%s)',currentDependant,currentPredictor,varNames{2},varNames{4}); % defining model structure
        lme = fitlme(current_data, formula); % fit LME
        txt = evalc('disp(lme)');                            
        fprintf(fid, '\n %s ~ %s \n -----------------------------\n\n %s\n\n',currentDependant ,currentPredictor, txt);
       
        %%% Plots: Linear Model and its checks 

        figure("Visible","on","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
        plot_linear_model(current_data,lme,currentDependant,currentPredictor,varNames,1,binwidth)
    
        sgtitle(sprintf("LME_%s_vs_%s",currentDependant,currentPredictor))
        print(gcf,fullfile(saveFolder,"Results",sprintf("LME_Saccade_Ripple_%s.ps",currentPredictor)),'-dpsc','-append','-fillpage')
        %exportgraphics(gcf,fullfile(saveFolder,"Results",sprintf("LME_%s_vs_%s.png",currentDependant,currentPredictor)))
        close(gcf)
    
        % if params.checks
        %     figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
        %     plot_diagnostic_lm(current_data, lme, currentDependant,currentPredictor)
        %     sgtitle(sprintf("LME_%s_vs_%s",currentDependant,currentPredictor))
        %     print(gcf,fullfile(saveFolder,"Results",sprintf("LME_Diagnostics_Saccade_Ripple_%s.ps",currentPredictor)),'-dpsc','-append','-fillpage')
        %     close(gcf)
        % end
    
    end
    psDir = fullfile(saveFolder,"Results",sprintf("LME_Saccade_Ripple_%s.ps",currentPredictor));
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
    
    psDir = fullfile(saveFolder,"Results",sprintf("LME_Diagnostics_Saccade_Ripple_%s.ps",currentPredictor));
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 

end
