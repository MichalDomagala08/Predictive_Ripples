  
    
    %%% ----- PARAMETERS AND SETUP -----
    
    % General parametes
    addpath("functions")
    params.ripple_detection_scheme_name = "RippleDetection_franz_et_al_RipplePeak_manualCorr";
    params.analysisName                 = "RippleDensity_Franz_reocurring_ripple";
    params.artif_chan_rej               = true; % If we reject artefactual channels
    params.blink_saccades_rej           = true; % if we are rejecting blink saccades
    params.trialNum = 80; % How many trials are there
    params.viz = 0;
    params.exact_analysisName = "RSA_bha";
    params.make_repr = 0;               % to make representaitons actually 
    params.make_rsa  = 0;               % if we are making RSA tables or using exsisting ones 
    params.data_preprocessing = 0;      % if the data needs preprocessing or if we load preprocessed data directly
    
    %%% Parameters of Representation
    % bha parameters
    params.bha          = 1;            % If we are to add bha signal to our representations.
    params.lowpassfreq  = 70; 
    params.highpassfreq = 150;
    
    % timing parameters
    params.before_cond =  5;          % N. samples before Sacc Onset that are not coinciding with saccades
    params.after_cond  =  5;          % N. samples after Sacc Offset!!!  that are nto coinciding with saccade
    params.max_time    = 200;         % maximul n of samples in which we will count ripple (for selected analyses)
    
    % fixation parameters:
    params.fix_box_length       = 13;   % length of ms of a boxcar for getting samples: Default 26 ms
    params.fix_box_max          = 75;   % The max cutoff of the length of fixation samples 
    params.fix_box_min          = 5;    % The mn cutoff of the length of fixation samples 
    params.fix_rip_trial        = 0;    % do not compute fix - rip RSA if they are in the same trial 
    params.fix_rip_distance_min = 125;  % minimum distance from fixation to ripple (to further account for too much overlap
    
    % Refix parameters
    params.p_px = 100; % Threshold in pixels 
    
    % saccadic parameters:
    params.sac_rip_beginning   = 500; % if we limit ourselves only to ripples that happened during trials AND NOT DURING first N ms! 
    params.sac_rip_coincidence = 20;  % Making sure that there is no Saccade - ripple coincidence (as it will drive the Correlation for sure) 
                                              % - If numeric, then use to count  ONLY when distance between saccade and ripple is greater than N samples:
                                              % - If "exact", get exactly when Beg and End of saccades and ripple are divergent
    %%%  Artifact Rejection 
    
    %Spike detection parameters
    params.spikeCtsThresh  = 2;      % Zawsze 1 na czas debugowania!
    params.spikePeakWin    = 0.3;   % Zwiększ okno dopasowania (hp) do 150ms
    params.spikeZThresh    = 4;    % Obniż próg (z-score na dużych danych rzadko dobija do 4 dla rozlazłych fal)
    params.spikeAmpScale   = 2.5;    % Obniż sumaryczny wymóg Peak-to-Trough
    params.spikeMNegPeakW  = 150;    % KLUCZ: Pozwól negatywnej fazie trwać do 300ms
    params.spikeTrackPeaks = false;   % Szukaj od pozytywnego (tak jak na Twoim obrazku)
    params.spikeWindow     = 250;    % Zwiększ margines wycinania artefaktu wokół IED
    
    % Other artifact rejection parameters
    params.artPadding = 70; % how many samples should we pad the artifacts with
    params.cluster_tolerance = 20;  % How many samples is it "close" 
    params.range_threshold = 6; % How many stds from mean is considered an artifact 
    params.iqr_w  = 3; %How many IQR do we have to surpass to have an artefact
    params.artif_mad = 8;
    analysisName = 'ArtifactRejection_alternative2'; % for artefact rejctions
    
    
    
    %%% ----- PATHS and DATALOAD ------
    
    %%% Paths
    addpath('D:\Documents_Dell\Predictive_Ripples_2025\analysis_codes\functions'); addpath("utilities\")
    saveFolder =  fullfile("D:\Documents_Dell\Predictive_Ripples_2025\Results",params.analysisName); mkdir(saveFolder); 
    saveDataFolder = fullfile("D:\Documents_Dell\Predictive_Ripples_2025\preprocessed",params.analysisName); mkdir(saveDataFolder)
    imageFolder = "D:\Documents_Dell\Predictive_Ripples_2025\data\stimuli";
    saccadeFile = 'D:\Documents_Dell\Predictive_Ripples_2025\data\saccadic_data\content_all_subjects_improved.mat';
    
    data = load(fullfile(saveDataFolder,"rippleDensity.mat"));
    saccades_table = data.saccades_table;

    fixation_table =  readtable(fullfile("D:\Documents_Dell\Predictive_Ripples_2025\DL_fixation_analysis\fixation_data.csv"));
    ripplewise_data =  readtable(fullfile("D:\Documents_Dell\Predictive_Ripples_2025\data\saccadic_data\ripplewise_table.xlsx"));

    %%% Load ripple data:(Table with rowise ripple and columns in which
    %%% occular event are they being contained 

    ft = renamevars(fixation_table, {'subject_new','trial_number'}, {'subject','trial'});
    
    rip_clean = ripplewise_data(ripplewise_data.is_fixation & ripplewise_data.is_repeated ...
                                 & ripplewise_data.is_beg & ripplewise_data.is_during_trial, :);
    
    [gs, subj_u]  = findgroups(ripplewise_data.subject);
    channel_num_u = splitapply(@(c) numel(unique(c)), ripplewise_data.channel, gs);
    channel_tbl   = table(subj_u, channel_num_u, 'VariableNames', {'subject','channel_num'});
    
    joined = outerjoin(ft, rip_clean, 'Keys', {'subject','trial','indx_sac'}, 'MergeKeys', true);
    joined.has_rip = ~isnan(joined.rip_num);
    
    counts = groupsummary(joined, {'name','subject','trial','indx_sac'}, 'sum', 'has_rip');
    counts.Properties.VariableNames{'sum_has_rip'} = 'ripple_count';
    
    fixation_ripples = outerjoin(ft, counts, 'Keys', {'name','subject','trial','indx_sac'}, 'MergeKeys', true);
    fixation_ripples.ripple_count(isnan(fixation_ripples.ripple_count)) = 0;
    
    fixation_ripples = outerjoin(fixation_ripples, channel_tbl, 'Keys', 'subject', 'MergeKeys', true, 'Type', 'left');
    
    fixation_ripples.mean_ripple_channum = fixation_ripples.ripple_count ./ fixation_ripples.channel_num;
    fixation_ripples.mean_ripple_fixlen  = fixation_ripples.ripple_count ./ fixation_ripples.fixation_len;
    
    fixation_ripples = sortrows(fixation_ripples, {'subject','trial','indx_sac'});
    fixation_ripples.has_rip =  fixation_ripples.ripple_count ~= 0;


    %%% Load Saccade Ripple data with Fix Before and After indexes:

    ripplewise_fix_after = readtable("D:\Documents_Dell\Predictive_Ripples_2025\data\saccadic_data\sac_rip_against_fix.csv");
    sacrip_fix = innerjoin(ripplewise_fix_after,ft,'LeftKeys',{'subject','trial','fix_idx'},'RightKeys',{'subject','trial','indx_sac'});
    
    %% BHA Computation:

    %%% Load Necessery Data, Ripple Tables and get names of raw files
    dataFolder =  fullfile("D:\Documents_Dell\Predictive_Ripples_2025\data\preprocessed");
    files = {dir(fullfile(dataFolder, '*.mat')).name};

    ripple_data = load(fullfile("D:\Documents_Dell\Predictive_Ripples_2025\preprocessed",params.ripple_detection_scheme_name,"rippleData.mat"));
    ripple_table = load(fullfile(saveDataFolder,"ripple_table.mat")).ripplewise_table;
    ripple_table = ripple_table(logical(ripple_table.is_repeated),:);

    fixation_bha_table = saccades_table(:,{'name','subject_new','trial_number','indx_sac','firstInTrial','lastInTrial','amplitude',...
                                            'saccade_onset_time','saccade_offset_time','fixation_beg','fixation_end','fixation_len'}); % Main fixation_bha_table

    [fixation_bha_table_ret,fixation_bha_long] = func_compute_fixwise_bha(fixation_bha_table,ripple_table,ripple_data,dataFolder);
  

%%


%%%%%%%%%%%%%%%
%%% ANALYSES %%
%%%%%%%%%%%%%%%

%% (1) Analysis: Ripple emergence in NN Representation

% tests whether there is increase in NN fix-rep Sim when accounting for
% Binary ripple:

gr_names = {'incremental_sim_resnet50','incremental_sim_dino','incremental_sim_cr_resnet50','incremental_sim_cr_dino','pure_sim_resnet50','pure_sim_dino','pure_sim_cr_resnet50','pure_sim_cr_dino'};
dist_titl = "Fixation Length vs Representation:"; box_title= "Fixation NN Rep: With vs Without Ripple"; 
formula_tpl= "%s ~ has_rip + fixation_len + indx_sac + (1|subject) + (1|name)";
group_var = {"fixation_len","has_rip","indx_sac","subject","name"};

pspath =fullfile(saveFolder,params.exact_analysisName,'ripples_nn_fix_ripples.ps');

curr_est = []; curr_tStat = []; curr_p = []; nam = []; pred = []; anname = [];
[curr_est, curr_tStat, curr_p, nam, pred, anname] = func_binary_testing(fixation_ripples(fixation_ripples.fixation_len < 2 & fixation_ripples.fixation_len >0,:),formula_tpl , gr_names, (2:4),"Fix_NN_rep_ripples", "has_rip",...
                                                        {'Without','With'},group_var,[1,1,0.2], dist_titl, "All",box_title, pspath);


system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(pspath,".ps", ".pdf "),pspath));    delete(pspath); 

T = table(anname,nam',pred,curr_est,curr_tStat,curr_p, ...
'VariableNames', {'anName','Dependent','Predictor','Coef','tStat','pValue'});

fn = fullfile(saveFolder,params.exact_analysisName,char('LME_NN_Fixation.xlsx')); sheet = matlab.lang.makeValidName('fix_NN_ripple_emergence');
writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
wb.Save(); wb.Close(false); ex.Quit(); delete(ex);


%% (2) Analysis: Ripple Count vs NN Representation progressio

% Both computed on Earlier and Late Layers, Different models and Incremental as well as fix to fix Cosine Similarity:

gr_names = {'incremental_sim_resnet50','incremental_sim_dino','incremental_sim_cr_resnet50','incremental_sim_cr_dino','pure_sim_resnet50','pure_sim_dino','pure_sim_cr_resnet50','pure_sim_cr_dino'};
pred_names = {'ripple_count','mean_ripple_channum','mean_ripple_fixlen'};

pspath =  fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_nnrep_count.ps',params.exact_analysisName));
dist_titl = "Ripple Count ~ Representation";
formula_tpl= "%s ~ %s + indx_sac +fixation_len + (1|subject) + (1|name)";
group_var = {"PLACEHOLDER","indx_sac","fixation_len","subject","name"}; coef_range = (2:4);

curr_est = []; curr_tStat = []; curr_p = []; nam = []; pred = []; anname = [];
for j = 1:numel(pred_names)
    currentPredictor       = string(pred_names{j});
    anname_val = sprintf("%s_vs_Representation",currentPredictor); 

    for i = 1:numel(gr_names)
       
        currentDependant      = string(gr_names{i}); dep_disp = strrep(currentDependant, '_', '-');
        group_var_curr = num2cell([currentDependant string(group_var(2:5))]);

        current_data = fixation_ripples(fixation_ripples.has_rip,:);
        robz_y = 0.6745 * ( current_data.(currentDependant) - median(current_data.(currentDependant),"omitnan")) ./ mad( current_data.(currentDependant), 1);
        current_data = current_data(abs(robz_y) <=4, :); % Próg około 4 robust-σ od mediany
            
    
        % --- Fit LME --------------------------------------------- ---------
        if currentPredictor == "ripple_count"
            lme = fitglme(current_data, sprintf(formula_tpl,currentPredictor, currentDependant), 'Distribution','Poisson');
        else
            lme = fitglme(current_data, sprintf(formula_tpl,currentPredictor, currentDependant), 'Distribution','Gamma', 'Link','log');
        end
        % --- Wyciągnięcie współczynników ----------------------------------
        curr_est   = [curr_est;   lme.Coefficients.Estimate(coef_range)];
        curr_tStat = [curr_tStat; lme.Coefficients.tStat(coef_range)];
        curr_p     = [curr_p;     lme.Coefficients.pValue(coef_range)];
        pred       = [pred;       string(lme.Coefficients.Name(coef_range))];
        nam        = [nam,        repmat(currentDependant, 1, numel(coef_range))];
        anname     = [anname;     repmat(string(anname_val), numel(coef_range), 1)];
        main_idx = find(contains(lme.Coefficients.Name, currentDependant), 1);
    
        % --- Figura -------------------------------------------------------
        figure("Visible","Off","PaperOrientation","landscape", ...
               "Units","normalized","Position",[0 0 1 1]);
        sgtitle(sprintf("%s: %s: P: %2.3f T: %2.3f", dist_titl, dep_disp,  lme.Coefficients.pValue(main_idx),  lme.Coefficients.tStat(main_idx)));
    
          %%% Plots: Linear Model and its checks 
    
        plot_linear_model_impr(current_data,lme,currentPredictor,currentDependant,group_var_curr,1,20,0)
        print(gcf, pspath, '-dpsc', '-append', '-fillpage');
        close(gcf);
    end
end

system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(pspath,".ps", ".pdf "),pspath));    delete(pspath); 



T = table(anname,nam',pred,curr_est,curr_tStat,curr_p, ...
'VariableNames', {'anName','Dependent','Predictor','Coef','tStat','pValue'});

fn = fullfile(saveFolder,params.exact_analysisName,char('LME_NN_Fixation.xlsx')); sheet = matlab.lang.makeValidName('fix_NN_ripple_count');
writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
wb.Save(); wb.Close(false); ex.Quit(); delete(ex);

%% (3) Check  NN rep in Fixation before vs After  Ripples happening during saccades:

    
    
gr_names = {'incremental_sim_resnet50','incremental_sim_dino','incremental_sim_cr_resnet50','incremental_sim_cr_dino','pure_sim_resnet50','pure_sim_dino','pure_sim_cr_resnet50','pure_sim_cr_dino'};
fromul = "%s ~ after_2 + fixation_len_ft+ fix_idx+ gap_edge +(1|subject) ";
pspath =  fullfile(saveFolder,params.exact_analysisName,'sac_ripples_nn_fix_after.ps');
dist_titl = "Fixation Length"; box_title= "In Saccade Before vs After Ripple to Fix RSA";
model_covars = {'fixation_len_ft',"after_2",'fix_idx','gap_edge',"subject"};

[curr_est1, curr_tStat1, curr_p1, nam1, pred1, anname1] = func_binary_testing(sacrip_fix(sacrip_fix.gap_edge < 25 & sacrip_fix.fixation_len_ft < 1 ,:),fromul , gr_names, (2:3),"In_Saccades_After",...
                                                        "after_2",{'before','after'},model_covars,[1,1,0.1], dist_titl, "All",box_title, pspath);
[curr_est2, curr_tStat2, curr_p2, nam2, pred2, anname2] = func_binary_testing(sacrip_fix((sacrip_fix.is_before | sacrip_fix.is_after) & sacrip_fix.gap_edge < 25 & sacrip_fix.fixation_len_ft < 1  ,:),...
                                                        fromul, gr_names, (2:3),"In_Saccades_After_Restricted","after_2",{'before','after'}, model_covars,[1,1,0.1], dist_titl, "Restriced",box_title,pspath);
curr_est = [curr_est1; curr_est2]; curr_tStat = [curr_tStat1; curr_tStat2]; curr_p = [curr_p1; curr_p2]; nam = [nam1, nam2]; pred = [pred1; pred2];anname = [anname1; anname2];


system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(pspath,".ps", ".pdf "),pspath));    delete(pspath); 

T = table(anname,nam',pred,curr_est,curr_tStat,curr_p, ...
'VariableNames', {'anName','Dependent','Predictor','Coef','tStat','pValue'});

fn = fullfile(saveFolder,params.exact_analysisName,char('LME_NN_Fixation.xlsx')); sheet = matlab.lang.makeValidName('sac_ripples_nn_fix_after');
writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
wb.Save(); wb.Close(false); ex.Quit(); delete(ex);

%% Check whether Ripple emerging In saccades results in le

