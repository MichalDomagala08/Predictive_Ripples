%%% This script is used to make ripple statistics and their visualisations

%%% --- PARAMETERS  AND PREPARATION--- %%%

% Naming parameter and stimuli sizes
params.artifact_rejection_scheme_name =  'ArtifactRejection_alternative2';
params.ripple_detection_scheme_name   = "RippleDetection_franz_et_al_RipplePeak_manualCorr";
params.analysisName                   = "RippleDensity_Franz_reocurring_ripple";
params.trialNum = 80;

% Analysis Progresion:
params.trialwise_ripple        = 1;   % Analysis of Per-trial Features predicted by ripple count,length and other features 
params.saccade_ripple_timing   = 1;   % Analysis of Timing of ripples during saccades (Complete with permutation testing) 
params.fixation_ripple_timing  = 1;   % Analysis of Timing of ripples during fixations
params.saccadewise_ripples     = 1;   % Testing saccading feataures in relation to ripples

% General analysis parameter
params.ripples_later     = true ;% only get ripples that are not in the first 500ms of the trial beginning to account for increase in energy fallacy
params.trialwise_ripples = true; % Whether we do things trialwise_ripples, or with inclusion ofchannels 
params.checks            = true;    % Wrtie model checks 
params.ripple_no_coincidence = true; % Remove from counting ripples that are coinciding with oneanother
%%% Visualisation parameters

%histogram vis parameters
params.bounds = 1;         
params.nbins = 100;                  



%%% Paths
addpath('D:\Documents_Dell\Predictive_Ripples_2025\analysis_codes\functions'); addpath("utilities\")
dataFolder =  fullfile("D:\Documents_Dell\Predictive_Ripples_2025\preprocessed",params.analysisName);
saveFolder =  fullfile("D:\Documents_Dell\Predictive_Ripples_2025\Results",params.analysisName);
mkdir(fullfile(saveFolder,"Results"));

% load saccadic wise data  (and ripple concat data) 
data = load(fullfile(dataFolder,"rippleDensity.mat"));
saccades_table = data.saccades_table;


dd = load("D:\Documents_Dell\Predictive_Ripples_2025\data\saccadic_data\subject_viewing_distance.mat");
% Load pure ripple data:
ripple_data = load(fullfile("D:\Documents_Dell\Predictive_Ripples_2025\preprocessed",params.ripple_detection_scheme_name,"rippleData.mat"));
% CHECK - if the loaded ripple data is coinciding wth alternative concatenation scheme 
% ripple_table2 = load(fullfile("D:\Documents_Dell\Predictive_Ripples_2025\preprocessed\RippleDensity_Franz","rippleDensity.mat")).ripple_table;
% ripple_table2 = sortrows(ripple_table2, {'Subject','ChannelName','Trial'}, {'ascend','ascend','ascend'});
% ripple_table2 = movevars( ripple_table2, {'Subject','ChannelName','Trial','RippleNumber','RippleNumber_cleaned','RippleLength','RippleLength_cleaned'}, 'Before', 1);

% Load Image Entropy data:
entropy_table = load(fullfile("D:\Documents_Dell\Predictive_Ripples_2025\data\entropy_results.mat")).entropyTable;

%Load ripplewise table
ripplewise_table        = load(fullfile(dataFolder,"ripple_table.mat")).ripplewise_table;
ripplewise_table.len    = ripplewise_table.end - ripplewise_table.beg; % compute ripple length
ripplewise_table.is_beg = ripplewise_table.median > 2250; %|  ripplewise_table.median < 1750; % get an indicator of riple outside of boundaries


%%% --- PREPROCESSING ---

%%% Entropy Preprocessing
aa = outerjoin(saccades_table,entropy_table,'LeftKeys','name','RightKeys','image');

entropy_vars  = {'name','subject_new','trial_number','edges','luminance','aws','dg2ncb','dg2cb','meaning'} ;  % wklej swoje nazwy
entropy_table = unique(aa(:, entropy_vars));
entropy_table = sortrows(entropy_table, {'subject_new','trial_number'}, {'ascend','ascend'});


 
% compute the diff between images of entropy
numCols = [4,5,6,7,8,9]; g = findgroups(entropy_table.subject_new); newMat = nan(height(entropy_table), numel(numCols));
for gg = 1:max(g)
    idx = find(g==gg);
    for k = 1:length(numCols)
        v = entropy_table{idx, numCols(k)};
        newMat(idx,k) = [nan; diff(v)];           % current - previous
    end
end
for k = 1:numel(numCols)
    entropy_table.([entropy_table.Properties.VariableNames{numCols(k)} '_diff']) = newMat(:,k);
end



%%% Ripple Preprocessing

% Get out the Reocuffing

if params.ripple_no_coincidence; ripplewise_table = ripplewise_table(logical(ripplewise_table.is_repeated),:); end
% Create Channelwise Ripple Data:
[g, s, c, t] = findgroups(ripplewise_table.subject, ripplewise_table.channel, ripplewise_table.trial);
counts = splitapply(@numel, ripplewise_table.subject, g);     % count rows per group
meanLen = splitapply(@(len) mean(len,'omitnan'),ripplewise_table.len , g); 
crit_count_per_ch = splitapply(@(x) sum(x), ripplewise_table.is_beg, g);
meanLen_if_beg = splitapply(@(len,mask) mean(len(mask),'omitnan'), ripplewise_table.len, ripplewise_table.is_beg, g);

channelwise_ripples = table(s, c, t, counts, crit_count_per_ch, meanLen, meanLen_if_beg, ...
    'VariableNames', {'subject','channel','trial','count','count_clean','mean_length','mean_length_clean'});

% fill out missing channel counts to account for no ripples in given channel x trial x subject:
ut_bad_channels
for sb = 1:length(ripple_data.all_hipothetical_ripples_cluster)
    for chan =  1:length(ripple_data.all_hipothetical_ripples_cluster{sb})

        if isempty(ripple_data.all_chanNames{sb}{chan});  continue; end
        if ~isempty(artfi_subj_chan{sb})
            if sum(ismember(ripple_data.all_chanNames{sb}{chan},artfi_subj_chan{sb})); continue;  end; 
        end
        for it =  1:length(ripple_data.all_hipothetical_ripples_cluster{sb}{chan})
            if isempty(ripple_data.all_hipothetical_ripples_cluster{sb}{chan}{it})
                Trow = table(string(ripple_data.all_subjNames{sb}), string(ripple_data.all_chanNames {sb}{chan}), it, 0, 0, nan, nan, ...
                            'VariableNames', channelwise_ripples.Properties.VariableNames);
                channelwise_ripples = [channelwise_ripples; Trow];
            end
        end
    end
end



% Create Trial wise ripple data:
[g, s, t] = findgroups(channelwise_ripples.subject,  channelwise_ripples.trial);
chanCount = splitapply(@numel, channelwise_ripples.subject, g);     % count rows per group
count = splitapply(@sum, channelwise_ripples.count, g);     % count rows per group
count_clean = splitapply(@sum, channelwise_ripples.count_clean, g);     % count rows per group
mean_count = splitapply(@(x) mean(x,'omitnan'),channelwise_ripples.count, g);  
mean_count_clean = splitapply(@(x) mean(x,'omitnan'),channelwise_ripples.count_clean, g);  
mean_length = splitapply(@(x) mean(x,'omitnan'), channelwise_ripples.mean_length, g);     % count rows per group
mean_length_clean = splitapply(@(x) mean(x,'omitnan'), channelwise_ripples.mean_length_clean, g);     % count rows per group

tralwise_ripples = table(s, t, chanCount,count, count_clean, mean_count, mean_count_clean, mean_length, mean_length_clean, ...
    'VariableNames', {'subject','trial','channel_num', 'count','count_clean','mean_count','mean_count_clean','mean_length','mean_length_clean'});


%%% Saccadic Preprocessing

if params.ripples_later % Create a mask excluding saccaes happeninh in the first second of recording
    additional_mask =  saccades_table.closestRipple > 1 & saccades_table.saccade_onset_time >1;
else
    additional_mask = ones(height(saccades_table));
end
saccades_table_cleaned = saccades_table(~isnan(saccades_table.latency) & ...
                                            ~isnan(saccades_table.blinkSac) & ...
                                            saccades_table.closestRipple ~=0 & ...
                                            ~isnan(saccades_table.closestRipple) & ...
                                            additional_mask,:);
saccades_table_cleaned.closestRipple_TimingDiffAbs = abs(saccades_table_cleaned.closestRipple_TimingDiff);
saccades_table_cleaned.closest_rip_to_fix_beg_diffAbs = abs(saccades_table_cleaned.closest_rip_to_fix_beg_diff);
saccades_table_cleaned.closest_rip_to_fix_half_diffAbs = abs(saccades_table_cleaned.closest_rip_to_fix_half_diff);


% Create Trial-wise saccadic data: (Agregating by mean and std-s) 
[g, s, t,img,blk] = findgroups(saccades_table_cleaned.subject_new, saccades_table_cleaned.trial_number, saccades_table_cleaned.name,  saccades_table_cleaned.blinkTrial);

sac_gr_data.subject = s; sac_gr_data.trial = t; sac_gr_data.image = img; sac_gr_data.blinkTrial = blk;
sac_gr_data.saccade_number = splitapply(@numel, saccades_table_cleaned.subject, g);
for col = ["latency","duration","amplitude",'peakVelocity','Xpx_original','Ypx_original','aws_r50px','aws_r50px_diff','lum_r50px','lum_r50px_diff']
    sac_gr_data.(sprintf("%s_mean",col)) = splitapply(@(x) mean(x,"omitnan"), saccades_table_cleaned.(col), g);     
    sac_gr_data.(sprintf("%s_std",col))  = splitapply(@(x) std(x,"omitnan"), saccades_table_cleaned.(col), g); 
end
saccade_group_table = movevars( struct2table(sac_gr_data), {'subject','trial','image','blinkTrial'}, 'Before', 1);

%%% join saccade and ripple trialwise tables:: 
keys = {'subject','trial'};
ripple_wise_trial = outerjoin(tralwise_ripples,entropy_table,'LeftKeys',{'subject','trial'},'RightKeys',{'subject_new','trial_number'},'MergeKeys', true);
ripple_wise_trial = renamevars(ripple_wise_trial, {'subject_subject_new','trial_trial_number'}, {'subject','trial'});
ripple_wise_trial = outerjoin(ripple_wise_trial,saccade_group_table,'Keys', keys, 'MergeKeys', true, 'RightVariables', setdiff(saccade_group_table.Properties.VariableNames, keys, 'stable'));
ripple_wise_trial  = ripple_wise_trial(~isnan(ripple_wise_trial.count) & ripple_wise_trial.trial ~=0 ,:); %&~isnan(ripple_wise_trial.aws_diff) & ~isnan(ripple_wise_trial.saccade_number)
[ripple_wise_trial.subject_id, ~] = findgroups(ripple_wise_trial.subject);




%%
% 
ripplewise_table_s = func_saccade_locked_to_ripples(saccades_table,ripplewise_table,0);
ripplewise_table_m = func_saccade_locked_to_ripples(saccades_table,ripplewise_table,1);
% ripplewise_table_ps = func_saccade_locked_to_ripples(saccades_table,ripplewise_table,2);

% Create Channelwise Ripple Data:

% For is Fixation

conds = {"is_before_trial",'is_after_trial','is_during_trial','is_saccades','is_fixation'};
channelwise_ripples.anname = repmat("all",height(channelwise_ripples),1);
channelwise_ripples_all = channelwise_ripples;

for cond = 1:5
    condition = conds{cond};

    curr_table_s = ripplewise_table_m(logical(ripplewise_table_m.(condition)),:);
    [g, s, c, t] = findgroups(curr_table_s.subject, curr_table_s.channel, curr_table_s.trial);

    counts_s = splitapply(@numel, curr_table_s.subject, g);     % count rows per group
    meanLen_s = splitapply(@(len) mean(len,'omitnan'),curr_table_s.len , g); 
    crit_count_per_ch_s = splitapply(@(x) sum(x), curr_table_s.is_beg, g);
    meanLen_if_beg_s = splitapply(@(len,mask) mean(len(mask),'omitnan'), curr_table_s.len,...
                                                                         curr_table_s.is_beg, g);

    channelwise_ripples_s = table(s, c, t, counts_s, crit_count_per_ch_s, meanLen_s, meanLen_if_beg_s, ...
    'VariableNames', {'subject','channel','trial','count','count_clean','mean_length','mean_length_clean'});

    % fill out missing channel counts to account for no ripples in given channel x trial x subject:
    for sb = 1:length(ripple_data.all_hipothetical_ripples_cluster)
        for chan =  1:length(ripple_data.all_hipothetical_ripples_cluster{sb})
    
            if isempty(ripple_data.all_chanNames{sb}{chan});  continue; end
            if ~isempty(artfi_subj_chan{sb})
                if sum(ismember(ripple_data.all_chanNames{sb}{chan},artfi_subj_chan{sb})); continue;  end; 
            end
            for it =  1:length(ripple_data.all_hipothetical_ripples_cluster{sb}{chan})
                if isempty(ripple_data.all_hipothetical_ripples_cluster{sb}{chan}{it})
                    Trow = table(string(ripple_data.all_subjNames{sb}), string(ripple_data.all_chanNames {sb}{chan}), it, 0, 0, nan, nan, ...
                                'VariableNames', channelwise_ripples_s.Properties.VariableNames);
                    channelwise_ripples_s = [channelwise_ripples_s; Trow];
                end
            end
        end
    end
    channelwise_ripples_s.anname = repmat(condition,height(channelwise_ripples_s),1);
    channelwise_ripples_all = [channelwise_ripples_all; channelwise_ripples_s];
end

% Create Trial wise ripple data:
[g, a, s, t] = findgroups(channelwise_ripples_all.anname,channelwise_ripples_all.subject,  channelwise_ripples_all.trial);
chanCount = splitapply(@numel, channelwise_ripples_all.subject, g);     % count rows per group
count = splitapply(@sum, channelwise_ripples_all.count, g);     % count rows per group
count_clean = splitapply(@sum, channelwise_ripples_all.count_clean, g);     % count rows per group
mean_count = splitapply(@(x) mean(x,'omitnan'),channelwise_ripples_all.count, g);  
mean_count_clean = splitapply(@(x) mean(x,'omitnan'),channelwise_ripples_all.count_clean, g);  
mean_length = splitapply(@(x) mean(x,'omitnan'), channelwise_ripples_all.mean_length, g);     % count rows per group
mean_length_clean = splitapply(@(x) mean(x,'omitnan'), channelwise_ripples_all.mean_length_clean, g);     % count rows per group

tralwise_ripples_s = table(a,s, t, chanCount,count, count_clean, mean_count, mean_count_clean, mean_length, mean_length_clean, ...
    'VariableNames', {'anname','subject','trial','channel_num', 'count','count_clean','mean_count','mean_count_clean','mean_length','mean_length_clean'});


%%% join saccade and ripple trialwise tables:: 
keys = {'subject','trial'};
ripple_wise_trial_s = outerjoin(tralwise_ripples_s,entropy_table,'LeftKeys',{'subject','trial'},'RightKeys',{'subject_new','trial_number'},'MergeKeys', true);
ripple_wise_trial_s = renamevars(ripple_wise_trial_s, {'subject_subject_new','trial_trial_number'}, {'subject','trial'});
ripple_wise_trial_s = outerjoin(ripple_wise_trial_s,saccade_group_table,'Keys', keys, 'MergeKeys', true, 'RightVariables', setdiff(saccade_group_table.Properties.VariableNames, keys, 'stable'));
ripple_wise_trial_s  = ripple_wise_trial_s(~isnan(ripple_wise_trial_s.count) & ripple_wise_trial_s.trial ~=0 ,:); %&~isnan(ripple_wise_trial.aws_diff) & ~isnan(ripple_wise_trial.saccade_number)
[ripple_wise_trial_s.subject_id, ~] = findgroups(ripple_wise_trial_s.subject);



%%



ripplewise_ent = outerjoin(ripplewise_table_m,entropy_table,'LeftKeys',{'subject','trial'},'RightKeys',{'subject_new','trial_number'});

for i = (45:56)
    figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
    plot_ripple_entropy_2d_hist(ripplewise_ent.Properties.VariableNames{i},ripplewise_ent)
    title("Ripple Count vs %s",ripplewise_ent.Properties.VariableNames{i});
    print(gcf,fullfile(saveFolder,"Results","Peri_ripple_hist_M.ps"),'-dpsc','-append','-fillpage')
    close(gcf)

end

psDir = fullfile(saveFolder,"Results","Peri_ripple_hist_M.ps");
system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
        
% Return variables:
% time_centers (1 x nb_time)
% aws_centers (1 x nb_aws)
% counts2d (nb_aws x nb_time)
% counts2d_colnorm (normalized)



%% Save Critical Tables in Excel:


aaa = ripplewise_table_m(ripplewise_table_m.subject == "NS127_02" & ripplewise_table_m.trial == 1 ,:)

writetable(ripplewise_table_m,"D:\Documents_Dell\Predictive_Ripples_2025\data\saccadic_data\ripplewise_table.xlsx")
writetable(saccades_table,"D:\Documents_Dell\Predictive_Ripples_2025\data\saccadic_data\saccades_table.xlsx")

%%



%%%% --- ANALYSIS (1) Saccadic Featue and Ripple Numbers ---
predictor_features = [7,4,5,6,8,9];
dependentFeatueres_sac = (25:45);
dependentFeatueres_img = (11:22);


if params.trialwise_ripple
    filename = 'lme_outputs_trial.txt';
    fid = fopen(filename,'w+');  
    fprintf(fid, '--- %s ---\n\n',datestr(now,'yyyy-mm-dd HH:MM:SS'));
    for i = predictor_features

    
        Coef_list = []; SE_list = []; T_list = []; P_list = []; CI_lo_list = []; CI_hi_list = []; Sig_list = [];
        Dep_list = {}; Pred_list = {};
        for j   = [dependentFeatueres_img dependentFeatueres_sac]

    
            currentDependant = ripple_wise_trial.Properties.VariableNames{j};    
            currentPredictor = ripple_wise_trial.Properties.VariableNames{i};    
            current_data = ripple_wise_trial(:, ["subject","subject_id","trial","image","channel_num",currentPredictor,"blinkTrial",currentDependant]);
            varNames = {currentPredictor,"trial","channel_num","subject"};
            fprintf("%s ~  %s \n",currentDependant,currentPredictor)

    
            %current_data.(currentDependant) = sqrt(abs(current_data.(currentDependant)));
                % Outlier removing (by robust MAD z-score  (4) 
            robz_x = 0.6745 * ( current_data.(currentPredictor) - median(current_data.(currentPredictor),"omitnan")) ./ mad( current_data.(currentPredictor), 1);
            robz_y = 0.6745 * ( current_data.(currentDependant) - median(current_data.(currentDependant),"omitnan")) ./ mad( current_data.(currentDependant), 1);

            current_data = current_data(abs(robz_x) <= 4 & abs(robz_y) <= 4, :); % Próg około 4 robust-σ od mediany
         
               %%% Model fiting
            formula = sprintf('%s ~ %s +trial +%s+ (1|subject)',currentDependant,currentPredictor,varNames{3}); % defining model structure
            if all(current_data.(currentDependant) > 0)
                lme = fitglme(current_data, formula, 'Distribution','Gamma', 'Link','log'); % fit LME
            else
                lme = fitlme(current_data, formula); % fit LME
            end
            txt = evalc('disp(lme)');                            
            fprintf(fid, '\n %s ~ %s \n -----------------------------\n\n %s\n\n',currentDependant ,currentPredictor, txt);
    


            %%% Save model coefs in Excel Table:
            idx_coeffs = find(strcmp(lme.CoefficientNames, currentPredictor), 1);
            Coef_list(end+1) = lme.Coefficients.Estimate(idx_coeffs); T_list(end+1) = lme.Coefficients.tStat(idx_coeffs); P_list(end+1) = lme.Coefficients.pValue(idx_coeffs);
            [ciLo, ciHi] = coefCI(lme); CI_lo_list(end+1) = ciLo(idx_coeffs);CI_hi_list(end+1) = ciHi(idx_coeffs);
            Dep_list{end+1} = currentDependant; Pred_list{end+1} = currentPredictor;

    
            %%% Plots: Linear Model and its checks 
            figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
            plot_linear_model(current_data,lme,currentDependant,currentPredictor,varNames,1,(max(current_data.(currentPredictor)) - min(current_data.(currentPredictor)))/10)
            sgtitle(sprintf("LME_%s_vs_%s_bin(%d)",currentDependant,currentPredictor,(max(current_data.(currentPredictor)) - min(current_data.(currentPredictor)))/10))
            print(gcf,fullfile(saveFolder,"Results",sprintf("LME_log_ripplewise_%s.ps",currentPredictor)),'-dpsc','-append','-fillpage')
            %exportgraphics(gcf,fullfile(saveFolder,"Results",sprintf("LME_%s_vs_%s.png",currentDependant,currentPredictor)))
            close(gcf)
        
            if params.checks
                figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
                plot_diagnostic_lm(current_data, lme, currentDependant,currentPredictor)
                sgtitle(sprintf("LME_%s_vs_%s",currentDependant,currentPredictor))
                print(gcf,fullfile(saveFolder,"Results",sprintf("LME_log_Diagnostic_ripplewise_%s.ps",currentPredictor)),'-dpsc','-append','-fillpage')
                close(gcf)
            end
            
        end

        %%% Save Model data to external excel:
        T = table(string(Dep_list(:)), string(Pred_list(:)), Coef_list(:), T_list(:), P_list(:), CI_lo_list(:), CI_hi_list(:), ...
        'VariableNames', {'Dependent','Predictor','Coef','tStat','pValue','CI_low','CI_high'});
        outfn = fullfile(saveFolder,'LME_log_ripple.xlsx');

        fn = fullfile(saveFolder,'LME_ripple.xlsx'); sheet = matlab.lang.makeValidName(char(currentPredictor));
        writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
        ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
        sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
        wb.Save(); wb.Close(false); ex.Quit(); delete(ex);

        psDir = fullfile(saveFolder,"Results",sprintf("LME_log_ripplewise_%s.ps",currentPredictor));
        system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
        
        psDir = fullfile(saveFolder,"Results",sprintf("LME_log_Diagnostic_ripplewise_%s.ps",currentPredictor));
        system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
    end
end



%% 
%%%% --- ANALYSIS (1.2) Saccadic Featue and Ripple Numbers --- EXTENDED  to Ripple before, after and wahtnot
predictor_features = [7,5,6,8,9,10];
dependentFeatueres_sac = (26:46);
dependentFeatueres_img = (12:23);
conds = {"is_before_trial",'is_after_trial','is_during_trial','is_saccades','is_fixation'};


if params.trialwise_ripple
    filename = 'lme_outputs_trial_d.txt';
    fid = fopen(filename,'w+');  
    fprintf(fid, '--- %s ---\n\n',datestr(now,'yyyy-mm-dd HH:MM:SS'));

    for c = 1:length(conds)
        ripple_wise_trial_t = ripple_wise_trial_s(ripple_wise_trial_s.anname == conds{c},:);
        for i = predictor_features
    
        
            Coef_list = []; SE_list = []; T_list = []; P_list = []; CI_lo_list = []; CI_hi_list = []; Sig_list = [];
            Dep_list = {}; Pred_list = {};
            for j   = [dependentFeatueres_img dependentFeatueres_sac]
    
        
                currentDependant = ripple_wise_trial_t.Properties.VariableNames{j};    
                currentPredictor = ripple_wise_trial_t.Properties.VariableNames{i};    
                current_data = ripple_wise_trial_t(:, ["subject","subject_id","trial","image","channel_num",currentPredictor,"blinkTrial",currentDependant]);
                varNames = {currentPredictor,"trial","channel_num","subject"};
                fprintf("%s ~  %s \n",currentDependant,currentPredictor)
    
        
                %current_data.(currentDependant) = sqrt(abs(current_data.(currentDependant)));
                    % Outlier removing (by robust MAD z-score  (4) 
                robz_x = 0.6745 * ( current_data.(currentPredictor) - median(current_data.(currentPredictor),"omitnan")) ./ mad( current_data.(currentPredictor), 1);
                robz_y = 0.6745 * ( current_data.(currentDependant) - median(current_data.(currentDependant),"omitnan")) ./ mad( current_data.(currentDependant), 1);
    
                current_data = current_data(abs(robz_x) <= 4 & abs(robz_y) <= 4, :); % Próg około 4 robust-σ od mediany
             
                if isempty(current_data); continue; end

                   %%% Model fiting
                formula = sprintf('%s ~ %s +trial +%s+ (1|subject)',currentDependant,currentPredictor,varNames{3}); % defining model structure
                if all(current_data.(currentDependant) > 0)
                    lme = fitglme(current_data, formula, 'Distribution','Gamma', 'Link','log'); % fit LME
                else
                    lme = fitlme(current_data, formula); % fit LME
                end
                txt = evalc('disp(lme)');                            
                fprintf(fid, '\n %s ~ %s \n -----------------------------\n\n %s\n\n',currentDependant ,currentPredictor, txt);
        
    
    
                %%% Save model coefs in Excel Table:
                idx_coeffs = find(strcmp(lme.CoefficientNames, currentPredictor), 1);
                Coef_list(end+1) = lme.Coefficients.Estimate(idx_coeffs); T_list(end+1) = lme.Coefficients.tStat(idx_coeffs); P_list(end+1) = lme.Coefficients.pValue(idx_coeffs);
                [ciLo, ciHi] = coefCI(lme); CI_lo_list(end+1) = ciLo(idx_coeffs);CI_hi_list(end+1) = ciHi(idx_coeffs);
                Dep_list{end+1} = currentDependant; Pred_list{end+1} = currentPredictor;
    
        
                %%% Plots: Linear Model and its checks 
                figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
                plot_linear_model(current_data,lme,currentDependant,currentPredictor,varNames,1,(max(current_data.(currentPredictor)) - min(current_data.(currentPredictor)))/10)
                sgtitle(sprintf("LME_%s_vs_%s_bin(%d)",currentDependant,currentPredictor,(max(current_data.(currentPredictor)) - min(current_data.(currentPredictor)))/10))
                print(gcf,fullfile(saveFolder,"Results","Extended_Trialwise",sprintf("LME_log_ripplewise_%s_%s_M.ps",conds{c},currentPredictor)),'-dpsc','-append','-fillpage')
                %exportgraphics(gcf,fullfile(saveFolder,"Results",sprintf("LME_%s_vs_%s.png",currentDependant,currentPredictor)))
                close(gcf)
            
                if params.checks
                    figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
                    plot_diagnostic_lm(current_data, lme, currentDependant,currentPredictor)
                    sgtitle(sprintf("LME_%s_vs_%s",currentDependant,currentPredictor))
                    print(gcf,fullfile(saveFolder,"Results","Extended_Trialwise",sprintf("LME_log_Diagnostic_ripplewise_%s_%s_M.ps",conds{c},currentPredictor)),'-dpsc','-append','-fillpage')
                    close(gcf)
                end
                
            end
    
            %%% Save Model data to external excel:
            T = table(string(Dep_list(:)), string(Pred_list(:)), Coef_list(:), T_list(:), P_list(:), CI_lo_list(:), CI_hi_list(:), ...
            'VariableNames', {'Dependent','Predictor','Coef','tStat','pValue','CI_low','CI_high'});
    
            fn = fullfile(saveFolder,"Results","Extended_Trialwise",sprintf('LME_ripple_%s_M.xlsx',conds{c})); sheet = matlab.lang.makeValidName(char(currentPredictor));
            writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
            ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
            sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
            wb.Save(); wb.Close(false); ex.Quit(); delete(ex);
    
            psDir = fullfile(saveFolder,"Results","Extended_Trialwise",sprintf("LME_log_ripplewise_%s_%s_M.ps",conds{c},currentPredictor));
            system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
            
            psDir = fullfile(saveFolder,"Results","Extended_Trialwise",sprintf("LME_log_Diagnostic_ripplewise_%s_%s_M.ps",conds{c},currentPredictor));
            system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
        end
    end
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

ripsac_n1 = {'rippleNumBefore_saccade','rippleNumIn_saccade','rippleNumAfter_saccade','rippleNumAfter2_saccade','rippleNumAfter3_saccade','rippleNumAfter4_saccade'};
ripsac_n2 = {'rippleNumBefore_saccade_M','rippleNumIn_saccade_M','rippleNumAfter_saccade_M','rippleNumAfter2_saccade_M','rippleNumAfter3_saccade_M','rippleNumAfter4_saccade_M'};

figure();
subplot(1,2,1)
bar(table2array(sum(saccades_table_cleaned(:,ripsac_n1),1)))
title("Ripple numbers in mean saccade Distance")
xticks([1,2,3,4,5,6])
xticklabels(["before",'Sac','Fix1','Fix2','Fix3','Fix4'])

subplot(1,2,2)
bar(table2array(sum(saccades_table_cleaned(:,ripsac_n2),1)))
title("For Mean Saccade Distance")
xticks([1,2,3,4,5,6])
xticklabels(["before",'Sac','Fix1','Fix2','Fix3','Fix4'])

%% (2.4) Median Split alongsie Saccadic Features 

% The feature include Amplitude, peakVelocity, Duration and all! 

%%% Warning! The Timing difference is not a good indicator because it will
%%% always have a concentraned near zero distribution!

% However, using COUNT! that changes things! 

mspl_sacc = saccades_table_cleaned;
%[6,7,20,21,50]
sac_names = {"latency","duration","amplitude","peakVelocity","fixation_len"};
for i = sac_names
    currentDependent = i{1};
    mspl_sacc.(currentDependent+"_mspl") = logical(saccades_table_cleaned.(currentDependent) >= median(saccades_table_cleaned.(currentDependent)));
    plot_barplot_sacc_timing_med_split(mspl_sacc,saveFolder,currentDependent,"Lower_Upper_Rip_Sac_bar",params)
    plot_histogram_sacc_timing_msplit(mspl_sacc,saveFolder,currentDependent,"closestRipple_TimingDiff","Lower_Upper_Rip_Sac_hist",params)
    plot_histogram_sacc_timing_msplit(mspl_sacc,saveFolder,currentDependent,"closest_rip_to_fix_half_diff","Lower_Upper_Rip_Fix_hist",params)
end
   
psDir = fullfile(saveFolder,"Results","Lower_Upper_Rip_Sac_bar.ps");
system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 

psDir = fullfile(saveFolder,"Results","Lower_Upper_Rip_Fix_hist.ps");
system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 

psDir = fullfile(saveFolder,"Results","Lower_Upper_Rip_Sac_hist.ps");
system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 



%%

%%%% ANALYSIS (3) Saccade Entrained Features

saccadic_features = {'latency','duration','amplitude','peakVelocity','Xpx_original','Ypx_original','saccade_onset_time','saccade_offset_time','aws_r50px','aws_r50px_diff',...
                    'lum_r50px','aws_r20px','lum_r20px','aws_r20px_diff','lum_r50px_diff','lum_r20px_diff','lum_r20px_diffabs','aws_r20px_diffabs','aws_r50px_diffabs','lum_r50px_diffabs'};
ripple_features = {'closestRipple_TimingDiff','closestRipple_Length','closest_rip_to_fix_half_diff','rippleAfterSaccade','random_sac_timing_index'};

%saccadic_features = [6,7,20,21,12,13,22,23,27:38]; % ripple_features = [41,42,54,71,73];


filename = 'lme_outputs_saccade.txt';
fid = fopen(filename,'w+');  
fprintf(fid, '--- %s ---\n\n',datestr(now,'yyyy-mm-dd HH:MM:SS'));


saccades_table_cleaned.image_num  = findgroups(saccades_table_cleaned.name);


for j = ripple_features
    Coef_list = []; SE_list = []; T_list = []; P_list = []; CI_lo_list = []; CI_hi_list = []; Sig_list = [];
        Dep_list = {}; Pred_list = {};
    for i = saccadic_features
        %%% Sub-matrix for currenr predictions
        currentPredictor = saccades_table_cleaned.Properties.VariableNames{j};
        currentDependant = saccades_table_cleaned.Properties.VariableNames{i};    

        fprintf("Current Analysis: %s ~%s\n",currentDependant,currentPredictor);
        current_data = saccades_table_cleaned(:, ["subject_new","image_num",currentPredictor,"blinkSac",currentDependant]);
        varNames = {currentPredictor,"image_num","blinkSac","subject_new"};
    
    
        % Outlier removing (by robust MAD z-score  (4) - The same for our Dependenta variable
        robz_y = 0.6745 * ( current_data.(currentDependant) - median(current_data.(currentDependant),"omitnan")) ./ mad( current_data.(currentDependant), 1);
        robz_x = 0.6745 * ( current_data.(currentPredictor) - median(current_data.(currentPredictor),"omitnan")) ./ mad( current_data.(currentPredictor), 1);
        current_data = current_data(abs(robz_x) <= 4 & abs(robz_y) <= 4, :); % Próg około 4 robust-σ od mediany
       
        
        %%% Model fiting
        formula = sprintf('%s ~ %s +%s+ (1|%s)',currentDependant,currentPredictor,image_num{2},varNames{4}); % defining model structure
        lme = fitlme(current_data, formula); % fit LME
        txt = evalc('disp(lme)');                            
        fprintf(fid, '\n %s ~ %s \n -----------------------------\n\n %s\n\n',currentDependant ,currentPredictor, txt);

        
        %%% Save model coefs in Excel Table:
        idx_coeffs = find(strcmp(lme.CoefficientNames, currentPredictor), 1);
        Coef_list(end+1) = lme.Coefficients.Estimate(idx_coeffs); T_list(end+1) = lme.Coefficients.tStat(idx_coeffs); P_list(end+1) = lme.Coefficients.pValue(idx_coeffs);
        [ciLo, ciHi] = coefCI(lme); CI_lo_list(end+1) = ciLo(idx_coeffs);CI_hi_list(end+1) = ciHi(idx_coeffs);
        Dep_list{end+1} = currentDependant; Pred_list{end+1} = currentPredictor;


       
        %%% Plots: Linear Model and its checks 

        figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
        plot_linear_model(current_data,lme,currentDependant,currentPredictor,varNames,1,(max(current_data.(currentPredictor)) - min(current_data.(currentPredictor)))/10)
    
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

   %%% Save Model data to external excel:
    T = table(string(Dep_list(:)), string(Pred_list(:)), Coef_list(:), T_list(:), P_list(:), CI_lo_list(:), CI_hi_list(:), ...
    'VariableNames', {'Dependent','Predictor','Coef','tStat','pValue','CI_low','CI_high'});
    outfn = fullfile(saveFolder,'LME_ripple.xlsx');

    fn = fullfile(saveFolder,'LME_saccades.xlsx'); sheet = matlab.lang.makeValidName(char(currentPredictor));
    writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
    ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
    sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
    wb.Save(); wb.Close(false); ex.Quit(); delete(ex);


    psDir = fullfile(saveFolder,"Results",sprintf("LME_Saccade_Ripple_%s.ps",currentPredictor));
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
    
    psDir = fullfile(saveFolder,"Results",sprintf("LME_Diagnostics_Saccade_Ripple_%s.ps",currentPredictor));
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 

end


%%


%%% ANALYSIS 4 - Bootstrapping statistics of analysis 

%% 4.1 Bootstrapping analyses:

mkdir(fullfile(saveFolder,"Results","SaccdesInRippleLME"))
saccadic_features = {'latency','duration','amplitude','peakVelocity','Xpx_original','Ypx_original','saccade_onset_time','saccade_offset_time','aws_r50px','aws_r50px_diff',...
                    'lum_r50px','aws_r20px','lum_r20px','aws_r20px_diff','lum_r50px_diff','lum_r20px_diff','lum_r20px_diffabs','aws_r20px_diffabs','aws_r50px_diffabs','lum_r50px_diffabs','closestRipple_L'};

func_bootstrapping_statistics(saccades_table_cleaned,ripple_data.all_subjNames,saccadic_features,"rippleNumIn_saccade",10000,fullfile(saveFolder,"Results","SaccdesInRippleLME"))
func_bootstrapping_statistics(saccades_table_cleaned,ripple_data.all_subjNames,saccadic_features,"rippleNumBefore_saccade",10000,fullfile(saveFolder,"Results","SaccdesInRippleLME"))
func_bootstrapping_statistics(saccades_table_cleaned,ripple_data.all_subjNames,saccadic_features,"rippleNumAfter_saccade",10000,fullfile(saveFolder,"Results","SaccdesInRippleLME"))
func_bootstrapping_statistics(saccades_table_cleaned,ripple_data.all_subjNames,saccadic_features,"rippleNumIn_saccade_M",10000,fullfile(saveFolder,"Results","SaccdesInRippleLME"))
func_bootstrapping_statistics(saccades_table_cleaned,ripple_data.all_subjNames,saccadic_features,"rippleNumBefore_saccade_M",10000,fullfile(saveFolder,"Results","SaccdesInRippleLME"))
func_bootstrapping_statistics(saccades_table_cleaned,ripple_data.all_subjNames,saccadic_features,"rippleNumAfter_saccade_M",10000,fullfile(saveFolder,"Results","SaccdesInRippleLME"))
func_bootstrapping_statistics(saccades_table_cleaned,ripple_data.all_subjNames,saccadic_features,"rippleAfterSaccade",10000,fullfile(saveFolder,"Results","SaccdesInRippleLME"))
func_bootstrapping_statistics(saccades_table_cleaned,ripple_data.all_subjNames,saccadic_features,"rippleBeforeSaccade",10000,fullfile(saveFolder,"Results","SaccdesInRippleLME"))
func_bootstrapping_statistics(saccades_table_cleaned,ripple_data.all_subjNames,saccadic_features,"rippleAroundSaccade",10000,fullfile(saveFolder,"Results","SaccdesInRippleLME"))



%% 4.2 LME based analysis


%%%% --- ANALYSIS (1) Saccadic Featue and Ripple Numbers ---



saccades_table_cleaned.rippleAroundSaccade_bin = saccades_table_cleaned.rippleAroundSaccade ~=0;
saccades_table_cleaned.rippleAfterSaccade_bin = saccades_table_cleaned.rippleAfterSaccade ~=0;
saccades_table_cleaned.rippleBeforeSaccade_bin = saccades_table_cleaned.rippleBeforeSaccade ~=0;

% predictor_features = [77,78,79,58,59,60,63,64,65,66,69,70,71,72];
predictor_features = {'rippleAroundSaccade_bin','rippleAfterSaccade_bin','rippleBeforeSaccade_bin','rippleNumBefore_saccade','rippleNumIn_saccade','rippleNumAfter_saccade','rippleNumAfter4_saccade',...
                      'rippleNumBefore_saccade_M','rippleNumIn_saccade_M','rippleNumAfter_saccade_M','rippleNumAfter4_saccade_M','rippleAroundSaccade','rippleAfterSaccade','rippleBeforeSaccade'} ;
saccadic_features = {'latency','duration','amplitude','peakVelocity','Xpx_original','Ypx_original','saccade_onset_time','saccade_offset_time','aws_r50px','aws_r50px_diff',...
                    'lum_r50px','aws_r20px','lum_r20px','aws_r20px_diff','lum_r50px_diff','lum_r20px_diff','lum_r20px_diffabs','aws_r20px_diffabs','aws_r50px_diffabs','lum_r50px_diffabs','closestRipple_L'};





if params.trialwise_ripple
    filename = 'lme_outputs_rippleTiming.txt';
    fid = fopen(filename,'w+');  
    fprintf(fid, '--- %s ---\n\n',datestr(now,'yyyy-mm-dd HH:MM:SS'));
    for i = predictor_features

    
        Coef_list = []; SE_list = []; T_list = []; P_list = []; CI_lo_list = []; CI_hi_list = []; Sig_list = [];
        Dep_list = {}; Pred_list = {};
        for j   = saccadic_features

    
            currentDependant = saccades_table_cleaned.Properties.VariableNames{j};    
            currentPredictor = saccades_table_cleaned.Properties.VariableNames{i}; 
            current_data  = saccades_table_cleaned(saccades_table_cleaned.trial_number ~= 0,:);
            current_data = current_data(:, ["subject_new","trial_number",currentPredictor,"blinkSac",currentDependant]);
            varNames = {"trial_number",currentPredictor,"blinkSac","subject_new"};

            fprintf("%s ~  %s \n",currentDependant,currentPredictor)
            
            % Outlier removing (by robust MAD z-score  (4) 
            robz_y = 0.6745 * ( current_data.(currentDependant) - median(current_data.(currentDependant),"omitnan")) ./ mad( current_data.(currentDependant), 1);
            current_data = current_data( abs(robz_y) <= 4, :); % Próg około 4 robust-σ od mediany
         
               %%% Model fiting

            formula = sprintf('%s ~ %s +%s +(1|%s)',currentDependant,currentPredictor,"trial_number","subject_new"); % defining model structure
            if all(current_data.(currentDependant) > 0)
                lme = fitglme(current_data, formula, 'Distribution','Gamma', 'Link','log'); % fit LME
            else
                lme = fitlme(current_data, formula); % fit LME
            end
            txt = evalc('disp(lme)');                            
            fprintf(fid, '\n %s ~ %s \n -----------------------------\n\n %s\n\n',currentDependant ,currentPredictor, txt);
    


            %%% Save model coefs in Excel Table:
            idx_coeffs = find(contains(lme.CoefficientNames, currentPredictor), 1);
            Coef_list(end+1) = lme.Coefficients.Estimate(idx_coeffs); T_list(end+1) = lme.Coefficients.tStat(idx_coeffs); P_list(end+1) = lme.Coefficients.pValue(idx_coeffs);
            [ciLo, ciHi] = coefCI(lme); CI_lo_list(end+1) = ciLo(idx_coeffs);CI_hi_list(end+1) = ciHi(idx_coeffs);
            Dep_list{end+1} = currentDependant; Pred_list{end+1} = currentPredictor;

    
            %%% Plots: Linear Model and its checks 
            figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
            plot_linear_model_bin(current_data,lme,currentDependant,currentPredictor,varNames,1,5)
            %sgtitle(sprintf("LME_%s_vs_%s_bin(%d)",currentDependant,currentPredictor,(max(current_data.(currentPredictor)) - min(current_data.(currentPredictor)))/10))
            print(gcf,fullfile(saveFolder,"Results","SaccdesInRippleLME",sprintf("LME_sac_rip_%s.ps",currentPredictor)),'-dpsc','-append','-fillpage')
            close(gcf)
        
            if params.checks
                figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
                plot_diagnostic_lm(current_data, lme, currentDependant,currentPredictor)
                sgtitle(sprintf(strrep("LME_%s_vs_%s",currentDependant,currentPredictor),'_','-'))
                print(gcf,fullfile(saveFolder,"Results","SaccdesInRippleLME",sprintf("LME_sac_rip_Diagnostic_ripplewise_%s.ps",currentPredictor)),'-dpsc','-append','-fillpage')
                close(gcf)
            end
            
        end

        %%% Save Model data to external excel:
        T = table(string(Dep_list(:)), string(Pred_list(:)), Coef_list(:), T_list(:), P_list(:), CI_lo_list(:), CI_hi_list(:), ...
        'VariableNames', {'Dependent','Predictor','Coef','tStat','pValue','CI_low','CI_high'});

        fn = fullfile(saveFolder,'LME_RippleInSaccades.xlsx'); sheet = matlab.lang.makeValidName(char(currentPredictor));
        writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
        ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
        sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
        wb.Save(); wb.Close(false); ex.Quit(); delete(ex);

        psDir = fullfile(saveFolder,"Results","SaccdesInRippleLME",sprintf("LME_sac_rip_%s.ps",currentPredictor));
        system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
        
        psDir = fullfile(saveFolder,"Results","SaccdesInRippleLME",sprintf("LME_sac_rip_Diagnostic_ripplewise_%s.ps",currentPredictor));
        system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
    end
end

%% 4.3 Before vs After vs In Saccades Features 


zipped_vars = {["rippleNumBefore_saccade","rippleNumIn_saccade"],["rippleNumBefore_saccade","rippleNumAfter_saccade"],["rippleNumAfter_saccade","rippleNumIn_saccade"],...
               ["rippleNumBefore_saccade_M","rippleNumIn_saccade_M"],["rippleNumBefore_saccade_M","rippleNumAfter_saccade_M"],["rippleNumAfter_saccade_M","rippleNumIn_saccade_M"],...
               ["rippleBeforeSaccade","rippleAroundSaccade"],["rippleBeforeSaccade","rippleAfterSaccade"],["rippleAfterSaccade","rippleAroundSaccade"]};
saccadic_features = {'latency','duration','amplitude','peakVelocity','Xpx_original','Ypx_original','saccade_onset_time','saccade_offset_time','aws_r50px','aws_r50px_diff',...
                    'lum_r50px','aws_r20px','lum_r20px','aws_r20px_diff','lum_r50px_diff','lum_r20px_diff','lum_r20px_diffabs','aws_r20px_diffabs','aws_r50px_diffabs','lum_r50px_diffabs'};


for i = 1:length(zipped_vars)
    Coef_list = []; SE_list = []; T_list = []; P_list = []; CI_lo_list = []; CI_hi_list = []; Sig_list = [];
    Dep_list = {}; Pred_list = {};
    for j   = saccadic_features
        currData = saccades_table_cleaned((saccades_table_cleaned.(zipped_vars{i}(1)) | saccades_table_cleaned.(zipped_vars{i}(2))) &...
                                     ~(saccades_table_cleaned.(zipped_vars{i}(1)) & saccades_table_cleaned.(zipped_vars{i}(2))),:);

        
        currentDependant = currData.Properties.VariableNames{j};    
        currentPredictor = zipped_vars{i}(1); 
        varNames = {"trial_number",currentPredictor,"blinkSac","subject_new"};

        formula = sprintf('%s ~ %s +%s +(1|%s)',currentDependant,currentPredictor,"trial_number","subject_new"); % defining model structure
        if all(currData.(currentDependant) > 0)
            lme = fitglme(currData, formula, 'Distribution','Gamma', 'Link','log'); % fit LME
        else
            lme = fitlme(currData, formula); % fit LME
        end



        %%% Save model coefs in Excel Table:
        idx_coeffs = find(contains(lme.CoefficientNames, currentPredictor), 1);
        Coef_list(end+1) = lme.Coefficients.Estimate(idx_coeffs); T_list(end+1) = lme.Coefficients.tStat(idx_coeffs); P_list(end+1) = lme.Coefficients.pValue(idx_coeffs);
        [ciLo, ciHi] = coefCI(lme); CI_lo_list(end+1) = ciLo(idx_coeffs);CI_hi_list(end+1) = ciHi(idx_coeffs);
        Dep_list{end+1} = currentDependant; Pred_list{end+1} = currentPredictor;


        %%% Plots: Linear Model and its checks 
        figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
        plot_linear_model_bin(currData,lme,currentDependant,currentPredictor,varNames,1,5)
        sgtitle(strrep(sprintf("LME_%s_vs_%s",zipped_vars{i}(1),zipped_vars{i}(2)),'_','-'));
        ax = gca; p = ax.Position;
        ax.Position = [p(1) p(2) p(3) p(4)*0.97];  % zmień 0.9 na mniejszą (np. 0.85) żeby więcej miejsca   
        print(gcf,fullfile(saveFolder,"Results","BefvsAft_LME",sprintf("LME_bef_aft_%s_%s.ps",zipped_vars{i}(1),zipped_vars{i}(2))),'-dpsc','-append','-fillpage')
        close(gcf)
    
        if params.checks
            figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
            plot_diagnostic_lm(currData, lme, currentDependant,currentPredictor)
            sgtitle(sprintf(strrep("LME_%s_vs_%s",zipped_vars{i}(1),zipped_vars{i}(2)),'_','-'))
            print(gcf,fullfile(saveFolder,"Results","BefvsAft_LME",sprintf("LME_bef_aft_Diagnostic_ripplewise_%s_%s.ps",zipped_vars{i}(1),zipped_vars{i}(2))),'-dpsc','-append','-fillpage')
            close(gcf)
        end
            
         
    end

    %%% Save Model data to external excel:
    T = table(string(Dep_list(:)), string(Pred_list(:)), Coef_list(:), T_list(:), P_list(:), CI_lo_list(:), CI_hi_list(:), ...
    'VariableNames', {'Dependent','Predictor','Coef','tStat','pValue','CI_low','CI_high'});

    sheet_name = char(zipped_vars{i}(1)+"_"+zipped_vars{i}(2));
    fn = fullfile(saveFolder,'LME_BeforeVsAfter.xlsx'); sheet = matlab.lang.makeValidName(char(sheet_name(10:16)+string(sheet_name(25:end-9))));
    writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
    ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
    sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
    wb.Save(); wb.Close(false); ex.Quit(); delete(ex);

    mkdir( fullfile(saveFolder,"Results","BefvsAft_LME"));
    psDir = fullfile(saveFolder,"Results","BefvsAft_LME",sprintf("LME_bef_aft_%s_%s.ps",zipped_vars{i}(1),zipped_vars{i}(2)));
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
    
    psDir = fullfile(saveFolder,"Results","BefvsAft_LME",sprintf("LME_bef_aft_Diagnostic_ripplewise_%s_%s.ps",zipped_vars{i}(1),zipped_vars{i}(2)));
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 


end



%%

%%% ANALYSIS (4) Saccade Locked Ripple Features: Whether Gaze is different
%%% during saccades if they happen during Ripple 




saccadic_features = {'latency','duration','amplitude','peakVelocity','Xpx_original','Ypx_original','saccade_onset_time','saccade_offset_time','aws_r50px','aws_r50px_diff',...
                    'lum_r50px','aws_r20px','lum_r20px','aws_r20px_diff','lum_r50px_diff','lum_r20px_diff','lum_r20px_diffabs','aws_r20px_diffabs','aws_r50px_diffabs','lum_r50px_diffabs'};
ripple_features   = {'rippleNumBefore_saccade_M','rippleNumIn_saccade_M','rippleNumAfter_saccade_M','rippleAroundSaccade','rippleAfterSaccade','rippleBeforeSaccade'}; %[64,65,66,70,71,72];

filename = 'lme_outputs_saccade.txt';
fid = fopen(filename,'w+');  
fprintf(fid, '--- %s ---\n\n',datestr(now,'yyyy-mm-dd HH:MM:SS'));


saccades_table_cleaned.image_num  = findgroups(saccades_table_cleaned.name);


for j = ripple_features
    Coef_list = []; SE_list = []; T_list = []; P_list = []; CI_lo_list = []; CI_hi_list = []; Sig_list = [];
        Dep_list = {}; Pred_list = {};
    for i = saccadic_features
        %%% Sub-matrix for currenr predictions
        currentPredictor = saccades_table_cleaned.Properties.VariableNames{j};
        currentDependant = saccades_table_cleaned.Properties.VariableNames{i};    

        fprintf("Current Analysis: %s ~%s\n",currentDependant,currentPredictor);
        current_data = saccades_table_cleaned(:, ["subject_new","image_num",currentPredictor,"blinkSac",currentDependant]);
        varNames = {"image_num",currentPredictor,"blinkSac","subject_new"};
    
    
        % Outlier removing (by robust MAD z-score  (4) - The same for our Dependenta variable
        robz_y = 0.6745 * ( current_data.(currentDependant) - median(current_data.(currentDependant),"omitnan")) ./ mad( current_data.(currentDependant), 1);
        %robz_x = 0.6745 * ( current_data.(currentPredictor) - median(current_data.(currentPredictor),"omitnan")) ./ mad( current_data.(currentPredictor), 1);
        current_data = current_data( abs(robz_y) <= 4, :); % Próg około 4 robust-σ od mediany
       
        
        %%% Model fiting
        formula = sprintf('%s ~ %s +%s+ (1|%s)',currentDependant,currentPredictor,varNames{2},varNames{4}); % defining model structure
        lme = fitlme(current_data, formula); % fit LME
        txt = evalc('disp(lme)');                            
        fprintf(fid, '\n %s ~ %s \n -----------------------------\n\n %s\n\n',currentDependant ,currentPredictor, txt);

        
        %%% Save model coefs in Excel Table:
        idx_coeffs = find(strcmp(lme.CoefficientNames, currentPredictor), 1);
        Coef_list(end+1) = lme.Coefficients.Estimate(idx_coeffs); T_list(end+1) = lme.Coefficients.tStat(idx_coeffs); P_list(end+1) = lme.Coefficients.pValue(idx_coeffs);
        [ciLo, ciHi] = coefCI(lme); CI_lo_list(end+1) = ciLo(idx_coeffs);CI_hi_list(end+1) = ciHi(idx_coeffs);
        Dep_list{end+1} = currentDependant; Pred_list{end+1} = currentPredictor;


       
        %%% Plots: Linear Model and its checks 
        %current_data.(currentPredictor) = categorical(current_data.(currentPredictor) ~= 0);
        figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
        plot_linear_model_bin(current_data,lme,currentDependant,currentPredictor,varNames,1,(max(current_data.(varNames{1})) - min(current_data.(varNames{1})))/10)

        sgtitle(sprintf("LME_%s_vs_%s",currentDependant,currentPredictor))
        print(gcf,fullfile(saveFolder,"Results",sprintf("LME_Saccade_Ripple_%s_Coincidence.ps",currentPredictor)),'-dpsc','-append','-fillpage')
        %exportgraphics(gcf,fullfile(saveFolder,"Results",sprintf("LME_%s_vs_%s.png",currentDependant,currentPredictor)))
        close(gcf)
    
        if params.checks
            figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
            plot_diagnostic_lm(current_data, lme, currentDependant,currentPredictor)
            sgtitle(sprintf("LME_%s_vs_%s",currentDependant,currentPredictor))
            print(gcf,fullfile(saveFolder,"Results",sprintf("LME_Diagnostics_Saccade_Ripple_%s_Coincidence.ps",currentPredictor)),'-dpsc','-append','-fillpage')
            close(gcf)
        end
    
    end

    %%% Save Model data to external excel:
    T = table(string(Dep_list(:)), string(Pred_list(:)), Coef_list(:), T_list(:), P_list(:), CI_lo_list(:), CI_hi_list(:), ...
    'VariableNames', {'Dependent','Predictor','Coef','tStat','pValue','CI_low','CI_high'});
    outfn = fullfile(saveFolder,'LME_ripple.xlsx');

    fn = fullfile(saveFolder,'LME_saccades_Coincidence.xlsx'); sheet = matlab.lang.makeValidName(char(currentPredictor));
    writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
    ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
    sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
    wb.Save(); wb.Close(false); ex.Quit(); delete(ex);

    psDir = fullfile(saveFolder,"Results",sprintf("LME_Saccade_Ripple_%s_Coincidence.ps",currentPredictor));
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
    
    psDir = fullfile(saveFolder,"Results",sprintf("LME_Diagnostics_Saccade_Ripple_%s_Coincidence.ps",currentPredictor));
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 

end




%%

%%% ANALYSIS (5) Ripple Locked Saccades - Properties of Ripples

% In this analysis we look at the diffeences in Properties of ripple
% nbetween those locked and Not locked to saccades: 

cleanes_current_ripplewise = ripplewise_table_m(ripplewise_table_m.is_repeated ...
                                              & ripplewise_table_m.is_during_trial ...
                                             & ripplewise_table_m.is_beg ...
                                             & (ripplewise_table_m.is_saccades | ripplewise_table_m.is_fixation),:);


cleanes_current_ripplewise = outerjoin(cleanes_current_ripplewise,entropy_table(:,1:3),'LeftKeys',{'subject','trial'},'RightKeys',{'subject_new','trial_number'});
cleanes_current_ripplewise.image_num = findgroups(cleanes_current_ripplewise.name);

Coef_list = []; SE_list = []; T_list = []; P_list = []; CI_lo_list = []; CI_hi_list = []; Sig_list = [];
Dep_list = {}; Pred_list = {};
        
entropy_pred = {'aws_r50px','aws_r50px_diff','lum_r50px','aws_r20px','lum_r20px','aws_r20px_diff','lum_r50px_diff','lum_r20px_diff',...
                'lum_r20px_diffabs','aws_r20px_diffabs','aws_r50px_diffabs','lum_r50px_diffabs'};

for i = entropy_pred
    currentPredictor = "is_fixation";
    currentDependent = i{1};
    varNames = {"image_num",currentPredictor,"subject"};


    formula = sprintf("%s ~ %s +%s +(1|%s)",currentDependent,currentPredictor,"image_num","subject");
    if all(cleanes_current_ripplewise.(currentDependent) > 0)
        lme = fitglme(cleanes_current_ripplewise, formula, 'Distribution','Gamma', 'Link','log'); % fit LME
    else
        lme = fitlme(cleanes_current_ripplewise, formula); % fit LME
    end


    %%% Save model coefs in Excel Table:
    idx_coeffs = find(contains(lme.CoefficientNames, currentPredictor), 1);
    Coef_list(end+1) = lme.Coefficients.Estimate(idx_coeffs); T_list(end+1) = lme.Coefficients.tStat(idx_coeffs); P_list(end+1) = lme.Coefficients.pValue(idx_coeffs);
    [ciLo, ciHi] = coefCI(lme); CI_lo_list(end+1) = ciLo(idx_coeffs);CI_hi_list(end+1) = ciHi(idx_coeffs);
    Dep_list{end+1} = currentDependent; Pred_list{end+1} = currentPredictor;


    %%% Plots: Linear Model and its checks 
    figure("Visible","On","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
    plot_linear_model_bin(ripplewise_ent,lme,currentDependent,currentPredictor,varNames,1,5)
    % sgtitle(strrep(sprintf("LME_%s_vs_%s",zipped_vars{i}(1),zipped_vars{i}(2)),'_','-'));
    % ax = gca; p = ax.Position;
    % ax.Position = [p(1) p(2) p(3) p(4)*0.97];  % zmień 0.9 na mniejszą (np. 0.85) żeby więcej miejsca   
    print(gcf,fullfile(saveFolder,"Results","RippleLocked_Sac_ImageParts",sprintf("LME_RippleLocked.ps")),'-dpsc','-append','-fillpage')
    close(gcf)

end
 %%
%%%% ANALYSIS (4) Saccadic features against the image-wise
% filename = 'lme_outputs_saccade_img.txt';
% 
% for j = [58:60,64:66]
%      Coef_list = []; SE_list = []; T_list = []; P_list = []; CI_lo_list = []; CI_hi_list = []; Sig_list = [];
%         Dep_list = {}; Pred_list = {};
%     for i = saccadic_features
%         %%% Sub-matrix for currenr predictions
%         currentPredictor = saccades_table_cleaned.Properties.VariableNames{j};
%         currentDependant = saccades_table_cleaned.Properties.VariableNames{i};    
% 
%         sprintf("Current Analysis: %s ~%s",currentDependant,currentPredictor);
%         current_data = saccades_table_cleaned(:, ["subject_new","trial_number",currentPredictor,"blinkSac",currentDependant]);
%         varNames = {currentPredictor,"trial_number","blinkSac","subject_new"};
% 
%         % 
%         % % Outlier removing (by robust MAD z-score  (4) 
%         % robz_x = 0.6745 * ( current_data.(currentPredictor) - median(current_data.(currentPredictor))) ./ mad( current_data.(currentPredictor), 1);
%         %  current_data = current_data(abs(robz_x) <= 4, :); % Próg około 4 robust-σ od mediany
% 
% 
%         %%% Model fiting
%         formula = sprintf('%s ~ %s +%s+ (1|%s)',currentDependant,currentPredictor,varNames{2},varNames{4}); % defining model structure
%         lme = fitlme(current_data, formula); % fit LME
%         txt = evalc('disp(lme)');                            
%         fprintf(fid, '\n %s ~ %s \n -----------------------------\n\n %s\n\n',currentDependant ,currentPredictor, txt);
% 
% 
% 
%         %%% Save model coefs in Excel Table:
%         idx_coeffs = find(strcmp(lme.CoefficientNames, currentPredictor), 1);
%         Coef_list(end+1) = lme.Coefficients.Estimate(idx_coeffs); T_list(end+1) = lme.Coefficients.tStat(idx_coeffs); P_list(end+1) = lme.Coefficients.pValue(idx_coeffs);
%         [ciLo, ciHi] = coefCI(lme); CI_lo_list(end+1) = ciLo(idx_coeffs);CI_hi_list(end+1) = ciHi(idx_coeffs);
%         Dep_list{end+1} = currentDependant; Pred_list{end+1} = currentPredictor;
% 
% 
% 
%         %%% Plots: Linear Model and its checks 
% 
%         figure("Visible","on","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
%         plot_linear_model(current_data,lme,currentDependant,currentPredictor,varNames,1,(max(current_data.(currentPredictor)) - min(current_data.(currentPredictor)))/10)
% 
%         sgtitle(sprintf("LME_%s_vs_%s",currentDependant,currentPredictor))
%         print(gcf,fullfile(saveFolder,"Results",sprintf("LME_Saccade_Ripple_%s.ps",currentPredictor)),'-dpsc','-append','-fillpage')
%         %exportgraphics(gcf,fullfile(saveFolder,"Results",sprintf("LME_%s_vs_%s.png",currentDependant,currentPredictor)))
%         close(gcf)
% 
%         % if params.checks
%         %     figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
%         %     plot_diagnostic_lm(current_data, lme, currentDependant,currentPredictor)
%         %     sgtitle(sprintf("LME_%s_vs_%s",currentDependant,currentPredictor))
%         %     print(gcf,fullfile(saveFolder,"Results",sprintf("LME_Diagnostics_Saccade_Ripple_%s.ps",currentPredictor)),'-dpsc','-append','-fillpage')
%         %     close(gcf)
%         % end
% 
%     end
% 
% 
%    %%% Save Model data to external excel:
%     T = table(string(Dep_list(:)), string(Pred_list(:)), Coef_list(:), T_list(:), P_list(:), CI_lo_list(:), CI_hi_list(:), ...
%     'VariableNames', {'Dependent','Predictor','Coef','tStat','pValue','CI_low','CI_high'});
%     outfn = fullfile(saveFolder,'LME_sac_img.xlsx');
% 
%     fn = fullfile(saveFolder,'LME_saccades.xlsx'); sheet = matlab.lang.makeValidName(char(currentPredictor));
%     writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
%     ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
%     sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
%     wb.Save(); wb.Close(false); ex.Quit(); delete(ex);
% 
% 
%     psDir = fullfile(saveFolder,"Results",sprintf("LME_Saccade_Ripple_%s.ps",currentPredictor));
%     system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
% 
%     psDir = fullfile(saveFolder,"Results",sprintf("LME_Diagnostics_Saccade_Ripple_%s.ps",currentPredictor));
%     system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
% 
% end
