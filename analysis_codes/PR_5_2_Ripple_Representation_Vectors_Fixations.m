    
    %%% This script computes representational vectors of a whole brain and
    %%% Saccade represenation
    
    
    
    %%% ----- PARAMETERS AND SETUP -----
    
    % General parametes
    addpath("functions")
    params.ripple_detection_scheme_name = "RippleDetection_franz_et_al_RipplePeak_manualCorr";
    params.analysisName                 = "RippleDensity_Franz_reocurring_ripple";
    params.artif_chan_rej               = true; % If we reject artefactual channels
    params.blink_saccades_rej           = true; % if we are rejecting blink saccades
    params.trialNum = 80; % How many trials are there
    params.viz = 0;
    params.exact_analysisName = "RSA_new_saccs";
    params.make_repr = 1;               % to make representaitons actually 
    params.make_rsa  = 1;               % if we are making RSA tables or using exsisting ones 
    params.data_preprocessing = 1;      % if the data needs preprocessing or if we load preprocessed data directly: 
                                        % 1 - Old Processing
                                        % 2 - Kasia's Processing
    params.data_save = 1; %Whether to save new data
    
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
    params.p_px = 75; % Threshold in pixels 
    
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
    dataFolder =  fullfile("D:\Documents_Dell\Predictive_Ripples_2025\data\reref");
    saveFolder =  fullfile("D:\Documents_Dell\Predictive_Ripples_2025\Results",params.analysisName); mkdir(saveFolder); 
    saveDataFolder = fullfile("D:\Documents_Dell\Predictive_Ripples_2025\preprocessed",params.analysisName); mkdir(saveDataFolder)
    imageFolder = "D:\Documents_Dell\Predictive_Ripples_2025\data\stimuli";
    saccadeFile = 'D:\Documents_Dell\Predictive_Ripples_2025\data\saccadic_data\content_all_subjects_improved.mat';
    
    
    data = load(fullfile(saveDataFolder,"rippleDensity.mat"));
    saccades_table = data.saccades_table;
    ut_subject_map
    [found, idx] = ismember(string(saccades_table.subject), string(subject_map(:,1)));
    saccades_table.subject_new = string(saccades_table.subject);             % domyślnie stara nazwa
    saccades_table.subject_new(found) = string(subject_map(idx(found), 3));  % podmiana wg mapy

    
    %%% Load ripple data:
    ripple_data = load(fullfile("D:\Documents_Dell\Predictive_Ripples_2025\preprocessed",params.ripple_detection_scheme_name,"rippleData.mat"));
    
    %Get bad-channels for artifact rejection purposes/
    ut_bad_channels
    
    % save parameters in txt file
    mkdir(fullfile(saveFolder,params.exact_analysisName))
    params_txt = evalc('disp(params)'); % Zamiana wyświetlenia struktury na tekst
    fileID = fopen(fullfile(saveFolder,params.exact_analysisName,'params.txt'), 'w'); % Zapis do pliku tekstowego
    fprintf(fileID, '%s', params_txt);
    fclose(fileID);
    
    
    %%% --- PREPROCESSING ---
    
    % Concatenate Ripples across channels withholding coincident ripples (for
    % both Peaks, Timing and Timecourse
    
    ripple_table = load(fullfile(saveDataFolder,"ripple_table.mat")).ripplewise_table;
    ripple_table = ripple_table(logical(ripple_table.is_repeated),:);
    
    %%% CREATE REPRESENTATIONAL DATA
    if params.make_repr 
        [subject_level_repr_rp, subject_level_repr_sac, subject_level_repr_fix] = func_create_representational_vectors(ripple_table,saccades_table,ripple_data,dataFolder,saveFolder,params);
    else
        if params.bha; reprTitle = "representational_data_bha.mat"; else; reprTitle = "representational_data.mat"; end
        repr = load(fullfile(saveFolder,params.exact_analysisName,reprTitle));
        subject_level_repr_rp = repr.subject_level_repr_rp;
        subject_level_repr_sac = repr.subject_level_repr_sac;
        subject_level_repr_fix = repr.subject_level_repr_fix;
    
    end
    
    %%% REPRESENTATIONAL SIMILATRITY: (1) Saccase X Ripple; (2) Ripple x Ripple; (3) Saccade x Saccade (4) Fixation X Ripple
    if params.make_rsa
        [T_all,T_rip,T_sac,T_fix,T_fix2] = func_create_representational_similarity_Rs(subject_level_repr_rp,subject_level_repr_sac,subject_level_repr_fix,saccades_table,ripple_table,ripple_data,saveFolder,params);
    else
        if params.bha; reprTitle = "RSA_tables_bha.mat"; else; reprTitle = "RSA_tables.mat"; end
    
        rsa = load(fullfile(saveFolder,params.exact_analysisName,reprTitle));
        T_all = rsa.T_all;
        T_fix = rsa.T_fix;
        T_rip = rsa.T_rip;
        T_sac = rsa.T_sac;
    end
    
    
    
    %%% GET FIXATION PRECEDENCE 
    saccades_table_refix = func_compute_refixation2(saccades_table,params);
    saccades_table_refix.fixation_beg = round((saccades_table_refix.fixation_beg *500)+2000);
    saccades_table_refix.fixation_end = round((saccades_table_refix.fixation_end *500)+2000);

    cols  = {'subject_new','name','trial_number','indx_sac','latency','duration','Xpx_original','Ypx_original','blinkSac','blinkTrial','firstInTrial','lastInTrial',...
             'amplitude','peakVelocity','saccade_onset_time','saccade_offset_time','aws_r50px','aws_r50px_diff','lum_r50px','aws_r20px','lum_r20px'...
             'aws_r20px_diff','lum_r50px_diff','lum_r20px_diff','lum_r20px_diffabs','aws_r20px_diffabs','aws_r50px_diffabs','lum_r50px_diffabs','fixation_beg','fixation_end'};
  
    fix_cols = {'subject','trial_number','indx_sac','Xpx_original','Ypx_original','subject_new','fixation_beg','fixation_end','fixation_len','is_refix','precursor_idx',...
                'precursor_distance','is_valid_refix','local_precursor_idx','local_lag','local_precursor_distance','precursor_lag','is_return_saccade'};

    T_fix_prefix = outerjoin(T_fix,saccades_table_refix(:,fix_cols),'LeftKeys',{'subject','trial','fix_beg','fix_end'},'RightKeys',{'subject_new','trial_number','fixation_beg','fixation_end'});
    T_fix_prefix(:,[88,89,90,93,94,95]) = [];
    T_fix_prefix = renamevars(T_fix_prefix,["subject_T_fix"],["subject"]);
    T_fix_prefix = T_fix_prefix(~isnan(T_fix_prefix.trial),:);
    
    %%% GET FIXATION VARIABLES 
    [T_fix_cleaned, T_fix_average, T_fix_tall]   = func_RSA_preprocessng_fixation(T_fix_prefix,params,1);
    
    
    %%% GET SACCADE VARIABLES 
    [sac_overlap,T_all_cleaned,Tall_cleaned,R_table] = func_RSA_preprocessing_saccades(T_all,ripple_table,params);
    
    %%% GET CLOSEST RIPPLE TO SACCADES
    R_unified= func_RSA_preprocessing_nearest_sac_fix(T_fix_prefix,T_all,ripple_table,params);
    
    %%% GET CLOSEST RIPPLE TO SACCADES (segments)
    R_unified_segments = func_RSA_preprocessing_nearest_sac_fix_segment(T_fix_tall,T_all,ripple_table,params);
    
    % %%% GET EQUALLY SPACED  (it takes time...)
    % R_unified_reverse  = func_create_representational_vectors_lite(R_unified,ripple_data,ripple_table,dataFolder,params);
    % R_unified_reverse.R_reverse_abs = abs(R_unified_reverse.R_reverse); 
    % R_unified_reverse.R_reverse_fisher = atanh(min(max(R_unified_reverse.R_reverse,-0.9999),0.9999)); 
    % R_unified_reverse(isnan(R_unified_reverse.R_reverse),:).R_reverse_fisher = nan(sum(isnan(R_unified_reverse.R_reverse)),1);
    % R_unified_reverse.R_reverse_fisher_abs = abs(R_unified_reverse.R_reverse_fisher);
    % 
    % % 1) Get Long format R_unified reverse for lienar model testing. 
    % has = ~isnan(R_unified_reverse.R_reverse);
    % T_fwd = R_unified_reverse; T_fwd.reversed = false(height(T_fwd),1); T_fwd.R = R_unified_reverse.R;
    % T_rev = R_unified_reverse(has,:); T_rev.reversed = true(height(T_rev),1); T_rev.R = R_unified_reverse.R_reverse(has);
    % R_unified_long = [T_fwd; T_rev];
    % R_unified_long.R_abs = abs(R_unified_long.R); R_unified_long.R_fisher = atanh(min(max(R_unified_long.R,-0.9999),0.9999)); R_unified_long.R_fisher_abs = abs(R_unified_long.R_fisher);
    
    
    %%% GET IN-SACCADES INFORMATION FOR RIPPLES: 
    saccades_table = load(fullfile("D:\Documents_Dell\Predictive_Ripples_2025\preprocessed",params.analysisName,"rippleDensity.mat")).saccades_table;
    ripplewise_table        = load(fullfile("D:\Documents_Dell\Predictive_Ripples_2025\preprocessed",params.analysisName,"ripple_table.mat")).ripplewise_table;
    ripplewise_table.len    = ripplewise_table.end - ripplewise_table.beg; % compute ripple length
    ripplewise_table.is_beg = ripplewise_table.median > 2250; %|  ripplewise_table.median < 1750; % get an indicator of riple outside of boundaries
    
    ripplewise_table_s = func_saccade_locked_to_ripples(saccades_table,ripplewise_table,0);

    ripplewise_segments = outerjoin(R_unified_segments,ripplewise_table_s,'LeftKeys',{'subject','channel','trial','rip_beg'},'RightKeys',{'subject','channel','trial','beg'});
    ripplewise = outerjoin(R_unified,ripplewise_table_s(ripplewise_table_s.is_repeated & ripplewise_table_s.is_beg,:),'LeftKeys',{'subject','channel','trial','rip_beg'},'RightKeys',{'subject','channel','trial','beg'});
    ripplewise = renamevars(ripplewise, {'subject_R_unified','channel_R_unified','trial_R_unified','rip_num_R_unified'}, {'subject','trial','channel','rip_num'});
    
    %Fixations: filter out Saccades, and Non prime elements:
    ripplewise_fix = ripplewise(~isnan(ripplewise.is_repeated) & ripplewise.is_fix& ripplewise.is_abs_closest,:);
    ripplewise_fix_segments = ripplewise_segments(~isnan(ripplewise_segments.is_repeated) & ripplewise_segments.is_fixseg& ripplewise_segments.is_abs_closest,:);
    
    % Saccades:filter out fixations and Non prime elements 
    ripplewise_sac = ripplewise(~isnan(ripplewise.is_repeated) &~ripplewise.is_fix& ripplewise.is_abs_closest,:);
    
    
    % Fixations with R of Ripple vs After and Before Fixation
    [~, T_fix_average_uncleared, T_fix_tall_uncleared]   = func_RSA_preprocessng_fixation(T_fix_prefix,params,0);
    
    % Compute absolute distance
    ripplewise2 = outerjoin(T_fix_average_uncleared,ripplewise_table_s(ripplewise_table_s.is_repeated & ripplewise_table_s.is_beg,:),'LeftKeys',{'subject','channel','trial','rip_beg'},'RightKeys',{'subject','channel','trial','beg'});
    ripplewise2 = renamevars(ripplewise2,  {'subject_T_fix_average_uncleared','trial_T_fix_average_uncleared','channel_T_fix_average_uncleared'}, {'subject','trial','channel'});
    ripplewise2 = ripplewise2(ripplewise2.is_saccades,:);
       
    % --- 1) before_2 / after_2 jako proste maski (overlap OK) ---
    ripplewise2.before_2 = ripplewise2.fix_beg < ripplewise2.rip_beg; ripplewise2.after_2  = ripplewise2.fix_end > ripplewise2.rip_end;
    
    % --- 2) Spójna odległość krawędziowa (zastępuje niespójne rip_fix_distance_edge) ---
    ripplewise2.gap_edge = max([ripplewise2.rip_beg - ripplewise2.fix_end, ripplewise2.fix_beg - ripplewise2.rip_end, zeros(height(ripplewise2),1)], [], 2);
    
    % --- 3) is_closest: per ripple (subject,trial,channel,rip_beg) osobno before / after ---
    ripplewise2.is_closest = false(height(ripplewise2),1);
    [g, ~] = findgroups(ripplewise2.subject, ripplewise2.trial, ripplewise2.channel, ripplewise2.rip_beg);
    for gg = 1:max(g)
        idx = find(g == gg);
        if isempty(idx), continue, end
    
        bidx = idx(ripplewise2.before_2(idx));        % wiersze "przed" w tym ripple
        if ~isempty(bidx)
            [~, loc] = min(ripplewise2.gap_edge(bidx));
            ripplewise2.is_closest(bidx(loc(1))) = true;
        end
    
        aidx = idx(ripplewise2.after_2(idx));         % wiersze "po" w tym ripple
        if ~isempty(aidx)
            [~, loc] = min(ripplewise2.gap_edge(aidx));
            ripplewise2.is_closest(aidx(loc(1))) = true;
        end
    end
    
    ripplewise_fix_after = ripplewise2(ripplewise2.is_closest,:);
    
    % The same but for segments: Getting Closest - Non Overlapping Segments to
    % the ripple at hand 
    
    
    ripplewise_segments = outerjoin(T_fix_tall_uncleared,ripplewise_table_s(ripplewise_table_s.is_repeated & ripplewise_table_s.is_beg,:),'LeftKeys',{'subject','channel','trial','rip_beg'},'RightKeys',{'subject','channel','trial','beg'});
    ripplewise_segments = renamevars(ripplewise_segments, {'subject_T_fix_tall_uncleared','trial_T_fix_tall_uncleared','channel_T_fix_tall_uncleared'}, {'subject','trial','channel'});
    ripplewise_segments = ripplewise_segments(ripplewise_segments.is_saccades,:);
       
    % --- 1) before_2 / after_2 jako proste maski (overlap OK) ---
    ripplewise_segments.before_2 = ripplewise_segments.seg_beg < ripplewise_segments.rip_beg; ripplewise_segments.after_2  = ripplewise_segments.seg_end > ripplewise_segments.rip_end;
    
    % --- 2) Spójna odległość krawędziowa (zastępuje niespójne rip_fix_distance_edge) ---
    ripplewise_segments.gap_edge = max([ripplewise_segments.rip_beg - ripplewise_segments.seg_end, ripplewise_segments.seg_beg - ripplewise_segments.rip_end, zeros(height(ripplewise_segments),1)], [], 2);
    
    % --- 3) is_closest: per ripple (subject,trial,channel,rip_beg) osobno before / after ---
    ripplewise_segments.is_closest = false(height(ripplewise_segments),1);
    [g, ~] = findgroups(ripplewise_segments.subject, ripplewise_segments.trial, ripplewise_segments.channel, ripplewise_segments.rip_beg);
    for gg = 1:max(g)
        idx = find(g == gg);
        if isempty(idx), continue, end
    
        bidx = idx(ripplewise_segments.before_2(idx));        % wiersze "przed" w tym ripple
        if ~isempty(bidx)
            [~, loc] = min(ripplewise_segments.gap_edge(bidx));
            ripplewise_segments.is_closest(bidx(loc(1))) = true;
        end
    
        aidx = idx(ripplewise_segments.after_2(idx));         % wiersze "po" w tym ripple
        if ~isempty(aidx)
            [~, loc] = min(ripplewise_segments.gap_edge(aidx));
            ripplewise_segments.is_closest(aidx(loc(1))) = true;
        end
    end
    
    ripplewise_seg_after = ripplewise_segments(ripplewise_segments.is_closest,:);
%%

    %%% Save All intermediate Tables in one excel:

    
    writetable(ripplewise,"D:\Documents_Dell\Predictive_Ripples_2025\data\saccadic_data\representation_data_ripplewise.csv")
    
    writetable(ripplewise_seg_after,"D:\Documents_Dell\Predictive_Ripples_2025\data\saccadic_data\sac_rip_against_seg.csv")
    writetable(ripplewise_fix_after,"D:\Documents_Dell\Predictive_Ripples_2025\data\saccadic_data\sac_rip_against_fix.csv")
    
    writetable(saccades_table_refix,"D:\Documents_Dell\Predictive_Ripples_2025\data\saccadic_data\saccades_table.xlsx")

    writetable(R_unified_reverse,"D:\Documents_Dell\Predictive_Ripples_2025\data\saccadic_data\representational_data_unified_reverse.csv")
    writetable(R_unified_long,"D:\Documents_Dell\Predictive_Ripples_2025\data\saccadic_data\representational_data_unified_long.csv")
    writetable(R_unified,"D:\Documents_Dell\Predictive_Ripples_2025\data\saccadic_data\representational_data_unified.csv")
    
    %% ======= ANALYSIS ======= %%%
    
    %% %%%  ----- ANALYSIS CLOSEST FIXATION vs OPPOSITE SPACE ------- %%%
    
    % (1) Closest Fixation  vs Opposite the same length
    R_unified_reverse_fix  = R_unified_long(R_unified_long.is_fix & ~isnan(R_unified_long.R) & ~isnan(R_unified_long.R_reverse) & ...
                                            R_unified_long.is_abs_closest & R_unified_long.rip_event_distance_edge >0,:);
    
    gr_names = {'R','R_abs','R_fisher','R_fisher_abs'}; fromul = "%s ~ reversed +  event_rip_latency_abs +(1|subject) +(1|rip_num)";
    pspath =  fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_reverse.ps',params.exact_analysisName));
    dist_titl = "Absolute Rip-Event distance (ms)"; box_title= "Fix-Rip RSA vs Reverse"; model_covars ={'event_rip_latency_abs',"reversed","subject","rip_num"};
    
    [curr_est1, curr_tStat1, curr_p1, nam1, pred1, anname1] = func_binary_testing(R_unified_reverse_fix,fromul , gr_names, (2:3),"Fixation_Reverse",...
                                                            "reversed",{'Fixation','Not Fixation'},model_covars,[1,0,250], dist_titl, "All",box_title, pspath);
    [curr_est2, curr_tStat2, curr_p2, nam2, pred2, anname2] = func_binary_testing(R_unified_reverse_fix(R_unified_reverse_fix.restricted,:), fromul, gr_names, (2:3),"Fixation_Reverse_Restricted",...
                                                            "reversed",{'Fixation','Not Fixation'}, model_covars,[1,1,20], dist_titl, "Restriced",box_title,pspath);
    curr_est = [curr_est1; curr_est2]; curr_tStat = [curr_tStat1; curr_tStat2]; curr_p = [curr_p1; curr_p2]; nam = [nam1, nam2]; pred = [pred1; pred2];anname = [anname1; anname2];
    
    
    
    psDir =fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_reverse.ps',params.exact_analysisName));
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
    
    T = table(anname,nam',pred,curr_est,curr_tStat,curr_p, ...
    'VariableNames', {'anName','Dependent','Predictor','Coef','tStat','pValue'});
    
    fn = fullfile(saveFolder,params.exact_analysisName,char(sprintf('LME_%s_Restr_Fixations_%s.xlsx',params.exact_analysisName,string(params.max_time)))); sheet = matlab.lang.makeValidName('fix_reverse_sim');
    writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
    ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
    sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
    wb.Save(); wb.Close(false); ex.Quit(); delete(ex);
    
    %%
    
    % (1) Closest Saccade  vs Opposite the same length
    
    R_unified_reverse_fix  = R_unified_long(~R_unified_long.is_fix& ~isnan(R_unified_long.R) & ~isnan(R_unified_long.R_reverse) & ...
                                            R_unified_long.is_abs_closest & R_unified_long.rip_event_distance_edge >0,:);
    
    gr_names = {'R','R_abs','R_fisher','R_fisher_abs'}; fromul = "%s ~ reversed +  event_rip_latency_abs +(1|subject) +(1|rip_num)";
    pspath =  fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_sac_reverse.ps',params.exact_analysisName));
    dist_titl = "Absolute Rip-Event distance (ms)"; box_title= "Sac-Rip RSA vs Reverse"; model_covars ={'event_rip_latency_abs',"reversed","subject","rip_num"};
    
    [curr_est1, curr_tStat1, curr_p1, nam1, pred1, anname1] = func_binary_testing(R_unified_reverse_fix,fromul , gr_names, (2:3),"Saccade_Reverse",...
                                                            "reversed",{'Saccade','Not Saccade'},model_covars,[1,0,250], dist_titl, "All",box_title, pspath);
    [curr_est2, curr_tStat2, curr_p2, nam2, pred2, anname2] = func_binary_testing(R_unified_reverse_fix(R_unified_reverse_fix.restricted,:), fromul, gr_names, (2:3),"Saccade_Reverse_Restricted",...
                                                            "reversed",{'Saccade','Not Saccade'}, model_covars,[1,1,20], dist_titl, "Restriced",box_title,pspath);
    curr_est = [curr_est1; curr_est2]; curr_tStat = [curr_tStat1; curr_tStat2]; curr_p = [curr_p1; curr_p2]; nam = [nam1, nam2]; pred = [pred1; pred2];anname = [anname1; anname2];
    
    
    psDir =fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_sac_reverse.ps',params.exact_analysisName));
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
    
    T = table(anname,nam',pred,curr_est,curr_tStat,curr_p, ...
    'VariableNames', {'anName','Dependent','Predictor','Coef','tStat','pValue'});
    
    fn = fullfile(saveFolder,params.exact_analysisName,char(sprintf('LME_%s_Restr_Fixations_%s.xlsx',params.exact_analysisName,string(params.max_time)))); sheet = matlab.lang.makeValidName('sac_reverse_sim');
    writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
    ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
    sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
    wb.Save(); wb.Close(false); ex.Quit(); delete(ex);
    
    
    %% %%%  ----- ANALYSIS CLOSEST FIXATION vs  CLOSEST SACCADES ------- %%%
    
    % (1) Closeste Fixation and Closest Saccade
    mask = R_unified.is_abs_closest & R_unified.rip_event_distance_edge >0;
    gr_names = {'R','R_abs','R_fisher','R_fisher_abs'}; fromul = "%s ~ is_fix +  event_rip_latency_abs +(1|subject)";
    pspath =  fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_vs_sac.ps',params.exact_analysisName));
    dist_titl = "Absolute Rip-Event distance (ms)"; box_title= "Sac-Rip vs Fix-Rip RSA"; model_covars ={'event_rip_latency_abs',"is_fix","subject"};
    
    [curr_est1, curr_tStat1, curr_p1, nam1, pred1, anname1] = func_binary_testing(R_unified(mask,:),fromul , gr_names, (2:3),"Fixation_vs_Saccades",...
                                                            "is_fix",{'Saccade','Fixation'},model_covars,[1,0,250], dist_titl, "All",box_title, pspath);
    [curr_est2, curr_tStat2, curr_p2, nam2, pred2, anname2] = func_binary_testing(R_unified(mask & R_unified.restricted ,:), fromul, gr_names, (2:3),"Fixation_vs_Saccades_Restricted",...
                                                            "is_fix",{'Saccade','Fixation'}, model_covars,[1,1,20], dist_titl, "Restriced",box_title,pspath);
    curr_est = [curr_est1; curr_est2]; curr_tStat = [curr_tStat1; curr_tStat2]; curr_p = [curr_p1; curr_p2]; nam = [nam1, nam2]; pred = [pred1; pred2];anname = [anname1; anname2];
    
    
    psDir =fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_vs_sac.ps',params.exact_analysisName));
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
    
    T = table(anname,nam',pred,curr_est,curr_tStat,curr_p, ...
    'VariableNames', {'anName','Dependent','Predictor','Coef','tStat','pValue'});
    
    fn = fullfile(saveFolder,params.exact_analysisName,char(sprintf('LME_%s_Restr_Fixations_%s.xlsx',params.exact_analysisName,string(params.max_time)))); sheet = matlab.lang.makeValidName('fixall_vs_sac_sim');
    writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
    ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
    sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
    wb.Save(); wb.Close(false); ex.Quit(); delete(ex);
    
    
    %%
    
    
    
    % (2) Closeste Fixation SEGMENT and Closest Saccade
    mask = R_unified_segments.is_abs_closest & R_unified_segments.rip_event_distance_edge >0;
    gr_names = {'R','R_abs','R_fisher','R_fisher_abs'}; fromul = "%s ~ is_fixseg + event_idx +  event_rip_latency_abs +(1|subject)";
    pspath =  fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_vs_sac_seg.ps',params.exact_analysisName));
    dist_titl = "Absolute Rip-Event distance (ms)"; box_title= "Sac-Rip vs Fix-Rip RSA Segment"; model_covars ={'event_rip_latency_abs',"is_fixseg","event_idx","subject"};
    
    [curr_est1, curr_tStat1, curr_p1, nam1, pred1, anname1] = func_binary_testing(R_unified_segments(mask,:),fromul , gr_names, (2:4),"Fixation_vs_Saccades_Segment",...
                                                            "is_fixseg",{'Saccade','Fixation'},model_covars,[1,0,250], dist_titl, "All",box_title, pspath);
    [curr_est2, curr_tStat2, curr_p2, nam2, pred2, anname2] = func_binary_testing(R_unified_segments(mask & R_unified_segments.restricted ,:), fromul, gr_names, (2:3),"Fixation_vs_Saccades_Segment_Restricted",...
                                                            "is_fixseg",{'Saccade','Fixation'}, model_covars,[1,1,20], dist_titl, "Restriced",box_title,pspath);
    curr_est = [curr_est1; curr_est2]; curr_tStat = [curr_tStat1; curr_tStat2]; curr_p = [curr_p1; curr_p2]; nam = [nam1, nam2]; pred = [pred1; pred2];anname = [anname1; anname2];
    
    
    
    psDir =fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_vs_sac_seg.ps',params.exact_analysisName));
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
    
    T = table(anname,nam',pred,curr_est,curr_tStat,curr_p, ...
    'VariableNames', {'anName','Dependent','Predictor','Coef','tStat','pValue'});
    
    fn = fullfile(saveFolder,params.exact_analysisName,char(sprintf('LME_%s_Restr_Fixations_%s.xlsx',params.exact_analysisName,string(params.max_time)))); sheet = matlab.lang.makeValidName('fixseg_vs_sac_sim');
    writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
    ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
    sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
    wb.Save(); wb.Close(false); ex.Quit(); delete(ex);
    
    
    %% %%%  ----- ANALYSIS FIXATIONS ------- %%%
    
    
    %%
    
    %%% (1) Are close ripple After vs Before Fixations Exhibiting different degree of Similarity 
    
    % Analysis done on both restricted as well as Far away:
    
    fix_overlap  = T_fix_average.rip_fix_timing_middle_abs >params.fix_rip_distance_min;
    
    gr_names = {'R','R_abs','R_fisher','R_fisher_abs'}; dist_titl = "Absolute Rip-Fix distance (ms)"; box_title= "Fixation  Before vs After Ripple"; 
    pspath =  fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_aft_vs_bef.ps',params.exact_analysisName));
    
    [curr_est1, curr_tStat1, curr_p1, nam1, pred1, anname1] = func_binary_testing(T_fix_average((T_fix_average.is_before | T_fix_average.is_after) & fix_overlap,:),...
                                                            "%s ~is_after +  rip_fix_distance_edge +(1|subject)" , gr_names, (2:3),"Fix_Bef_Aft", "is_after",...
                                                            {'Before','After'},{'rip_fix_distance_edge',"is_after","subject"},[1,0,250], dist_titl, "All",box_title, pspath);
    [curr_est2, curr_tStat2, curr_p2, nam2, pred2, anname2] = func_binary_testing(T_fix_average((T_fix_average.is_before_restricted | T_fix_average.is_after_restricted) & fix_overlap,:),...
                                                            "%s ~is_after_restricted +  rip_fix_distance_edge +(1|subject)", gr_names, (2:3),"Fix_Bef_Aft_Restricted","is_after_restricted",...
                                                            {'Before','After'}, {'rip_fix_distance_edge',"is_after_restricted","subject"},[1,1,20], dist_titl, "Restriced",box_title,pspath);
    curr_est = [curr_est1; curr_est2]; curr_tStat = [curr_tStat1; curr_tStat2]; curr_p = [curr_p1; curr_p2]; nam = [nam1, nam2]; pred = [pred1; pred2];anname = [anname1; anname2];
    
    
    psDir =fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_aft_vs_bef.ps',params.exact_analysisName));
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
    
    T = table(anname,nam',pred,curr_est,curr_tStat,curr_p, ...
    'VariableNames', {'anName','Dependent','Predictor','Coef','tStat','pValue'});
    
    fn = fullfile(saveFolder,params.exact_analysisName,char(sprintf('LME_%s_Restr_Fixations_%s.xlsx',params.exact_analysisName,string(params.max_time)))); sheet = matlab.lang.makeValidName('fixations_aft_bef');
    writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
    ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
    sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
    wb.Save(); wb.Close(false); ex.Quit(); delete(ex);
    
    
    %%
    %%% (2) Are close ripple After vs Before Fixations Exhibiting different degree of Similarity  (But on segment level:
    
    fix_overlap  = T_fix_tall.seg_to_ripple_distance_mid_abs >0;
    
    gr_names = {'R','R_abs','R_fisher','R_fisher_abs'}; dist_titl = "Absolute Rip-Event distance (ms)"; box_title= "Fixation Segment Before vs After Ripple"; 
    pspath =  fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_aft_vs_bef_seg.ps',params.exact_analysisName));
    
    [curr_est1, curr_tStat1, curr_p1, nam1, pred1, anname1] = func_binary_testing(T_fix_tall((T_fix_tall.is_before | T_fix_tall.is_after) & fix_overlap,:),...
                                                            "%s ~ is_after+ segment+ seg_to_ripple_distance_mid_abs +(1|subject)" , gr_names, (2:4),"Fix_Seg_Bef_Aft", "is_after",...
                                                            {'Before','After'},{'seg_to_ripple_distance_mid_abs',"is_after","segment","subject"},[1,0,250], dist_titl, "All",box_title, pspath);
    [curr_est2, curr_tStat2, curr_p2, nam2, pred2, anname2] = func_binary_testing(T_fix_tall((T_fix_tall.is_before_restricted | T_fix_tall.is_after_restricted) & fix_overlap,:),...
                                                            "%s ~ is_after_restricted+ segment+ seg_to_ripple_distance_mid_abs +(1|subject)", gr_names, (2:4),"Fix_Seg_Bef_Aft_Restricted","is_after_restricted",...
                                                            {'Before','After'}, {'seg_to_ripple_distance_mid_abs',"is_after_restricted","segment","subject"},[1,0,20], dist_titl, "Restriced",box_title,pspath);
    curr_est = [curr_est1; curr_est2]; curr_tStat = [curr_tStat1; curr_tStat2]; curr_p = [curr_p1; curr_p2]; nam = [nam1, nam2]; pred = [pred1; pred2];anname = [anname1; anname2];
    
    
    psDir =fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_aft_vs_bef_seg.ps',params.exact_analysisName));
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
    
    T = table(anname,nam',pred,curr_est,curr_tStat,curr_p, ...
    'VariableNames', {'anName','Dependent','Predictor','Coef','tStat','pValue'});
    
    fn = fullfile(saveFolder,params.exact_analysisName,char(sprintf('LME_%s_Restr_Fixations_%s.xlsx',params.exact_analysisName,string(params.max_time)))); sheet = matlab.lang.makeValidName('fixations_aft_bef_seg');
    writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
    ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
    sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
    wb.Save(); wb.Close(false); ex.Quit(); delete(ex);
    
    
    %%
    
    %%% (3) Is Ripple-Fixation similarity particular in specific part of
    %%% saccades? 
    
    fix_overlap  = T_fix_tall.seg_to_ripple_distance_mid_abs >0;
    T_fix_tall.fix_phase_cat = categorical(T_fix_tall.fix_phase);
    
    gr_names = {'R','R_abs','R_fisher','R_fisher_abs'}; dist_titl = "Absolute Saccade to Ripple distance (ms)"; box_title_tpl = "After - Before Boxplot %s";
    pspath   = fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_phase.ps',params.exact_analysisName));
    formula  = "%s ~ fix_phase +  seg_to_ripple_distance_mid_abs +(1|subject)"; pairs    = {{1,2},{2,3},{3,4},{4,5}};
    covars   = {'seg_to_ripple_distance_mid_abs',"fix_phase","subject"}; labels   = {'early2','early1','middle','late1','late2'};
    
    mask_all   = (T_fix_tall.is_before | T_fix_tall.is_after) & fix_overlap;
    mask_restr = mask_all & T_fix_tall.seg_to_ripple_distance_mid_abs < 250;
    
    % --- All -----------------------------------------------------------
    [curr_est1, curr_tStat1, curr_p1, nam1, pred1, anname1] = func_multilevel_testing(T_fix_tall(mask_all,:), formula, gr_names, (2:3), "RipProx","fix_phase",...
                                                            labels, pairs, "fix_phase", covars, [1,0,250], dist_titl, "All", box_title_tpl, "", "fix_phase_cat", pspath,[]);
    
    % --- Restricted (boxplot na full data, LME na <250 subset) ---------
    [curr_est2, curr_tStat2, curr_p2, nam2, pred2, anname2] = func_multilevel_testing(  T_fix_tall(mask_all,:), formula, gr_names, (2:3), "RipProx_Restr",  "fix_phase",... 
                    labels, pairs, "fix_phase", covars, [1,0,20], dist_titl, "Restricted", box_title_tpl, "Restricted ", "fix_phase_cat", pspath, T_fix_tall(mask_restr,:));
    
    % --- Concatenate + save --------------------------------------------
    curr_est = [curr_est1; curr_est2]; curr_tStat = [curr_tStat1; curr_tStat2]; curr_p = [curr_p1; curr_p2];
    nam = [nam1, nam2]; pred = [pred1; pred2]; anname = [anname1; anname2];
    
    
    psDir =fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_phase.ps',params.exact_analysisName));
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
    
    T = table(anname,nam',pred,curr_est,curr_tStat,curr_p, ...
    'VariableNames', {'anName','Dependent','Predictor','Coef','tStat','pValue'});
    
    fn = fullfile(saveFolder,params.exact_analysisName,char(sprintf('LME_%s_Restr_Fixations_%s.xlsx',params.exact_analysisName,string(params.max_time)))); sheet = matlab.lang.makeValidName('fixations_phase');
    writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
    ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
    sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
    wb.Save(); wb.Close(false); ex.Quit(); delete(ex);
    
    
    %% %%% ---- SACCADES ---- %%%
    
    
    %% (1) Representational Similarity vs Before,After in 1st, 2nd ,3rd saccade since ripple 
    
    % this is made with simple progressve tests
    groups = [[17,13,9,21,25,29];[39,35,31,43,47,51];[40,36,32,44,48,52];[18,14,10,22,26,30];[41,37,33,45,49,53];[42,38,34,46,50,54]];
    gr_names = {'R','R_abs','R_abs_log','R_fisher','R_fisher_abs','R_fisher_abs_log'};
    testinf = {{'b',3,2},{'b',2,1},{{'a','b'},1,1},{'a',1,2},{'a',2,3}};
    
    all_an_names = []; vs1 = []; vs2 = []; all_p = []; all_tStat = []; all_est = [];
    for i = 1:6
        curr_r_table =R_table(:,[groups(i,:)]);
        curr_est = []; curr_tStat = []; curr_p = [];
        for j = 1:5
            aaaa = Tall_cleaned(ismember(Tall_cleaned.side,testinf{j}{1}) & (Tall_cleaned.rank == testinf{j}{2} | Tall_cleaned.rank == testinf{j}{3}),:);
            if j ==3;  lme = fitlme(aaaa,sprintf("%s ~ side + sac_rip_latency + (1|subject)",gr_names{i}));
            else;      lme = fitlme(aaaa,sprintf("%s ~ rank + sac_rip_latency + (1|subject)",gr_names{i})); end
    
            vs1 = [vs1,string(curr_r_table.Properties.VariableNames{j})]; vs2 = [vs2,string(curr_r_table.Properties.VariableNames{j+1})];
            all_an_names = [all_an_names; string(gr_names{i})];
    
            curr_est = [curr_est,lme.Coefficients.Estimate(2)]; curr_tStat = [curr_tStat,lme.Coefficients.tStat(2)]; curr_p = [curr_p,lme.Coefficients.pValue(2)];
        end
    
        all_est = [all_est,curr_est]; all_tStat = [all_tStat,curr_tStat]; all_p = [all_p, curr_p];
        figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
    
        boxplot(table2array(curr_r_table)); hold on;
        for j = 1:5
            if curr_p(j) < 0.05; star = ' *'; else; star = ''; end
            y = max(table2array(curr_r_table(:,j:j+1)),[],'all') + 0.05;
            plot([j j j+1 j+1],[y y+0.02 y+0.02 y],'k','LineWidth',1.2)
            text(mean([j j+1]), y+0.1, sprintf('p=%.3g%s', curr_p(j), star), 'HorizontalAlignment','center')
        end
        title(strrep(gr_names{i},'_','-'));
        ylim([min(min(table2array(curr_r_table))), 1.5*max(max(table2array(curr_r_table)))*sign(max(max(table2array(curr_r_table))))])
    
        print(gcf,fullfile(saveFolder,params.exact_analysisName,'RSA_histograms_saccades.ps'),'-dpsc','-append','-fillpage');
    end
    
    
    %%% Save Model data to external excel:
    T = table(all_an_names,vs1',vs2',all_est',all_tStat',all_p', ...
    'VariableNames', {'anName','FirstContr','SecPred','Coef','tStat','pValue'});
    
    fn = fullfile(saveFolder,params.exact_analysisName,char(sprintf('LME_%s_Restr_%s_Coinc_%s.xlsx',params.exact_analysisName,string(params.max_time),string(params.sac_rip_coincidence)))); sheet = matlab.lang.makeValidName('histograms_saccades');
    writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
    ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
    sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
    wb.Save(); wb.Close(false); ex.Quit(); delete(ex);
    
    
    psDir =fullfile(saveFolder,params.exact_analysisName,'RSA_histograms_saccades.ps');
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
          
    %%
    
    
    % (2) Side + Trial to see whether the effects of size is indpndend from trial:
    curr_est = []; curr_tStat = []; curr_p = []; nam = []; dependent = []; pred = [];
    
    for i = 1:6
        lme = fitlme(Tall_cleaned,sprintf("%s ~ side+rank_cat + sac_rip_latency  + (1|subject) ",gr_names{i}));
        curr_est = [curr_est;lme.Coefficients.Estimate(2:5)]; curr_tStat = [curr_tStat;lme.Coefficients.tStat(2:5)]; curr_p = [curr_p;lme.Coefficients.pValue(2:5)];
        pred = [pred; string(lme.Coefficients.Name(2:5))]; nam = [nam, string(gr_names{i}), string(gr_names{i}), string(gr_names{i}), string(gr_names{i})];
    end
    
    T = table(pred,nam',curr_est,curr_tStat,curr_p, ...
    'VariableNames', {'Dependent','anName','Coef','tStat','pValue'});
    
    fn = fullfile(saveFolder,params.exact_analysisName,char(sprintf('LME_%s_Restr_%s_Coinc_%s.xlsx',params.exact_analysisName,string(params.max_time),string(params.sac_rip_coincidence)))); sheet = matlab.lang.makeValidName('saccades_1vs1');
    writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
    ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
    sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
    wb.Save(); wb.Close(false); ex.Quit(); delete(ex);
    
    %%
    
    
    % (3) Interaction side x distance 
    
    curr_est = []; curr_tStat = []; curr_p = []; nam = []; pred = [];
    for i = 1:6
        lme = fitlme(Tall_cleaned,sprintf("%s ~ side*rank_cat + sac_rip_latency  + (1|subject) ",gr_names{i}));
        curr_est = [curr_est;lme.Coefficients.Estimate(2:7)]; curr_tStat = [curr_tStat;lme.Coefficients.tStat(2:7)]; curr_p = [curr_p;lme.Coefficients.pValue(2:7)];
        pred = [pred; string(lme.Coefficients.Name(2:7))]; nam = [nam, string(gr_names{i}), string(gr_names{i}), string(gr_names{i}), string(gr_names{i}), string(gr_names{i}), string(gr_names{i})];
    end
    
    T = table(nam',pred,curr_est,curr_tStat,curr_p, ...
    'VariableNames', {'Dependent','anName','Coef','tStat','pValue'});
    
    fn = fullfile(saveFolder,params.exact_analysisName,char(sprintf('LME_%s_Restr_%s_Coinc_%s.xlsx',params.exact_analysisName,string(params.max_time),string(params.sac_rip_coincidence)))); sheet = matlab.lang.makeValidName('saccades_all_models');
    writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
    ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
    sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
    wb.Save(); wb.Close(false); ex.Quit(); delete(ex);
    
    %%
    
    %%% (4) Are close ripple After vs Before Saccade Exhibiting different degree of Similarity 
    
    % Those analyses Are taking ALL SACCADES to Ripples and checking ON average
    % instead of taking the closest ones! 
    
    
    
    gr_names = {'R','R_abs','R_abs_log','R_fisher','R_fisher_abs','R_fisher_abs_log'}; dist_titl = "Absolute Rip-Sac distance (ms)"; box_title= "Saccade Before vs After Ripple"; 
    pspath =  fullfile(saveFolder,params.exact_analysisName,'RSA_ripples_aft_vs_bef');
    
    [curr_est1, curr_tStat1, curr_p1, nam1, pred1, anname1] = func_binary_testing(T_all_cleaned((T_all_cleaned.is_before | T_all_cleaned.is_after) & ~sac_overlap,:),...
                                                            "%s ~ is_after +  rip_sac_timing_abs +(1|subject) " , gr_names, (2:3),"Ripple_Saccade_Proximity",...
                                                            "is_after",{'After','Before'},{'rip_sac_timing_abs',"is_after","subject"},[1,0,250], dist_titl, "All",box_title, pspath);
    [curr_est2, curr_tStat2, curr_p2, nam2, pred2, anname2] = func_binary_testing(T_all_cleaned((T_all_cleaned.is_before_restricted | T_all_cleaned.is_after_restricted) & ~sac_overlap,:),...
                                                            "%s ~ is_after_restricted +  rip_sac_timing_abs +(1|subject) ", gr_names, (2:3),"Ripple_Saccade_Proximity_restricted",...
                                                            "is_after_restricted",{'Not in saccades','In saccades'},  {'rip_sac_timing_abs',"is_after_restricted","subject"},[1,1,20], dist_titl, "Restriced",box_title,pspath);
    curr_est = [curr_est1; curr_est2]; curr_tStat = [curr_tStat1; curr_tStat2]; curr_p = [curr_p1; curr_p2]; nam = [nam1, nam2]; pred = [pred1; pred2];anname = [anname1; anname2];
    
    
    
    psDir =fullfile(saveFolder,params.exact_analysisName,'RSA_ripples_aft_vs_bef.ps');
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
    
    T = table(anname,nam',pred,curr_est,curr_tStat,curr_p, ...
    'VariableNames', {'anName','Dependent','Predictor','Coef','tStat','pValue'});
    
    fn = fullfile(saveFolder,params.exact_analysisName,char(sprintf('LME_%s_Restr_%s_Coinc_%s.xlsx',params.exact_analysisName,string(params.max_time),string(params.sac_rip_coincidence)))); sheet = matlab.lang.makeValidName('saccades_aft_vs_bef');
    writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
    ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
    sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
    wb.Save(); wb.Close(false); ex.Quit(); delete(ex);
    
    
    %% %%% ---- RIPPLES IN SACCADES  ---- %%%
    
    
    %% % (1) Ripples In saccades vs Other
    
    
    gr_names = {'R','R_abs','R_fisher','R_fisher_abs'}; fromul = "%s ~ is_saccades +  event_rip_latency_abs +(1|subject) ";
    pspath =  fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_in_saccades_fix_rsa.ps',params.exact_analysisName));
    dist_titl = "Absolute Rip-Event distance (ms)"; box_title= "In vs Out Saccade Ripple to Fix RSA"; model_covars = {'event_rip_latency_abs',"is_saccades","subject"};
    
    [curr_est1, curr_tStat1, curr_p1, nam1, pred1, anname1] = func_binary_testing(ripplewise_fix,fromul , gr_names, (2:3),"Ripples_In_Saccades",...
                                                            "is_saccades",{'Not in saccades','In saccades'},model_covars,[1,0,250], dist_titl, "All",box_title, pspath);
    [curr_est2, curr_tStat2, curr_p2, nam2, pred2, anname2] = func_binary_testing(ripplewise_fix(ripplewise_fix.restricted,:), fromul, gr_names, (2:3),"Ripples_In_Saccades_Restricted",...
                                                            "is_saccades",{'Not in saccades','In saccades'}, model_covars,[1,1,20], dist_titl, "Restriced",box_title,pspath);
    curr_est = [curr_est1; curr_est2]; curr_tStat = [curr_tStat1; curr_tStat2]; curr_p = [curr_p1; curr_p2]; nam = [nam1, nam2]; pred = [pred1; pred2];anname = [anname1; anname2];
    
    
    psDir =fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_in_saccades_fix_rsa.ps',params.exact_analysisName));
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
    
    T = table(anname,nam',pred,curr_est,curr_tStat,curr_p, ...
    'VariableNames', {'anName','Dependent','Predictor','Coef','tStat','pValue'});
    
    fn = fullfile(saveFolder,params.exact_analysisName,char(sprintf('LME_%s_Restr_Fixations_%s.xlsx',params.exact_analysisName,string(params.max_time)))); sheet = matlab.lang.makeValidName('in_saccades_fix_rsa');
    writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
    ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
    sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
    wb.Save(); wb.Close(false); ex.Quit(); delete(ex);
    
    
    
    
    %% % (2) Ripples In saccades: Fixation Before vs Fixation After
    
    
    gr_names = {'R','R_abs','R_fisher','R_fisher_abs'}; fromul = "%s ~ after_2 +  gap_edge +(1|subject) ";
    pspath =  fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_in_saccades_fix_aft_bef.ps',params.exact_analysisName));
    dist_titl = "Absolute Rip-Event distance (ms)"; box_title= "In Saccade Before vs After Ripple to Fix RSA"; model_covars = {'gap_edge',"after_2","subject"};
    
    [curr_est1, curr_tStat1, curr_p1, nam1, pred1, anname1] = func_binary_testing(ripplewise_fix_after(ripplewise_fix_after.gap_edge < 25 ,:),fromul , gr_names, (2:3),"In_Saccades_After",...
                                                            "after_2",{'before','after'},model_covars,[1,1,5], dist_titl, "All",box_title, pspath);
    [curr_est2, curr_tStat2, curr_p2, nam2, pred2, anname2] = func_binary_testing(ripplewise_fix_after((ripplewise_fix_after.is_before | ripplewise_fix_after.is_after) & ripplewise_fix_after.gap_edge < 25  ,:), fromul, gr_names, (2:3),"In_Saccades_After_Restricted",...
                                                            "after_2",{'before','after'}, model_covars,[1,1,5], dist_titl, "Restriced",box_title,pspath);
    
    curr_est = [curr_est1; curr_est2]; curr_tStat = [curr_tStat1; curr_tStat2]; curr_p = [curr_p1; curr_p2]; nam = [nam1, nam2]; pred = [pred1; pred2];anname = [anname1; anname2];
    
    
    psDir =fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_in_saccades_fix_aft_bef.ps',params.exact_analysisName));
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
    
    T = table(anname,nam',pred,curr_est,curr_tStat,curr_p, ...
    'VariableNames', {'anName','Dependent','Predictor','Coef','tStat','pValue'});
    
    fn = fullfile(saveFolder,params.exact_analysisName,char(sprintf('LME_%s_Restr_Fixations_%s.xlsx',params.exact_analysisName,string(params.max_time)))); sheet = matlab.lang.makeValidName('in_saccades_fix_aft_bef');
    writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
    ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
    sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
    wb.Save(); wb.Close(false); ex.Quit(); delete(ex);
    
    
    %% % (3) Ripples In saccades: Fixation Before vs Fixation Aftet (Segments) 
    
    
    gr_names = {'R','R_abs','R_fisher','R_fisher_abs'}; fromul = "%s ~ after_2 +  gap_edge +(1|subject) ";
    pspath =  fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_in_saccades_seg_aft_bef.ps',params.exact_analysisName));
    dist_titl = "Absolute Rip-Event distance (ms)"; box_title= "In Saccade Before vs After Ripple to Seg RSA"; model_covars = {'gap_edge',"after_2","subject"};
    
    [curr_est1, curr_tStat1, curr_p1, nam1, pred1, anname1] = func_binary_testing(ripplewise_seg_after(ripplewise_seg_after.seg_to_ripple_distance_mid_abs < 25 ,:),fromul , gr_names, (2:3),"In_Saccades_After",...
                                                            "after_2",{'before','after'},model_covars,[1,1,5], dist_titl, "All",box_title, pspath);
    [curr_est2, curr_tStat2, curr_p2, nam2, pred2, anname2] = func_binary_testing(ripplewise_seg_after(ripplewise_seg_after.is_before | ripplewise_seg_after.is_after & ripplewise_seg_after.seg_to_ripple_distance_mid_abs < 25 ,:), fromul, gr_names, (2:3),"In_Saccades_After_Restricted",...
                                                            "after_2",{'before','after'}, model_covars,[1,1,5], dist_titl, "Restriced",box_title,pspath);
    curr_est = [curr_est1; curr_est2]; curr_tStat = [curr_tStat1; curr_tStat2]; curr_p = [curr_p1; curr_p2]; nam = [nam1, nam2]; pred = [pred1; pred2];anname = [anname1; anname2];
    
    
    psDir =fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_in_saccades_seg_aft_bef.ps',params.exact_analysisName));
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
    
    T = table(anname,nam',pred,curr_est,curr_tStat,curr_p, ...
    'VariableNames', {'anName','Dependent','Predictor','Coef','tStat','pValue'});
    
    fn = fullfile(saveFolder,params.exact_analysisName,char(sprintf('LME_%s_Restr_Fixations_%s.xlsx',params.exact_analysisName,string(params.max_time)))); sheet = matlab.lang.makeValidName('in_saccades_seg_aft_bef');
    writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
    ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
    sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
    wb.Save(); wb.Close(false); ex.Quit(); delete(ex);
    
    
    
    %% %%% ---- RIPPLES and RETURN FIXATIONS  ---- %%%
    
    
    %% 1 Are Ripple to Fixation similarity greater for when ripple happen during Return Fixations? 


    %%% Remember! It is not about closest but All to All in case of trial!
    %%% and whole similarity 
    
    fix_overlap  = T_fix_average.rip_fix_timing_middle_abs >params.fix_rip_distance_min;
    T_fix_average.is_return_light = T_fix_average.precursor_lag > 1;
    T_fix_average.is_refix_exact = T_fix_average.local_lag == 1;

    
    gr_names = {'R','R_abs','R_fisher','R_fisher_abs'}; dist_titl = "Absolute Rip-Fix distance (ms)"; box_title= "Fixation  Primary vs Return Ripple"; 
    pspath =  fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_prime_vs_refix.ps',params.exact_analysisName));
    
    [curr_est1, curr_tStat1, curr_p1, nam1, pred1, anname1] = func_binary_testing(T_fix_average((T_fix_average.is_before | T_fix_average.is_after) & fix_overlap,:),...
                                                            "%s ~is_refix +  rip_fix_distance_edge +(1|subject)" , gr_names, (2:3),"Fix_Prime_Refix", "is_refix",...
                                                            {'Primary Fix','Refix'},{'rip_fix_distance_edge',"is_refix","subject"},[1,0,250], dist_titl, "All",box_title, pspath);
    [curr_est2, curr_tStat2, curr_p2, nam2, pred2, anname2] = func_binary_testing(T_fix_average((T_fix_average.is_before_restricted | T_fix_average.is_after_restricted) & fix_overlap,:),...
                                                            "%s ~is_refix +  rip_fix_distance_edge +(1|subject)", gr_names, (2:3),"Fix_Prime_Refix_Restricted","is_refix",...
                                                            {'Primary Fix','Refix'}, {'rip_fix_distance_edge',"is_refix","subject"},[1,1,20], dist_titl, "Restriced",box_title,pspath);
    [curr_est3, curr_tStat3, curr_p3, nam3, pred3, anname3] = func_binary_testing(T_fix_average((T_fix_average.is_before | T_fix_average.is_before) & fix_overlap,:),...
                                                            "%s ~is_return_saccade +  rip_fix_distance_edge +(1|subject)", gr_names, (2:3),"Fix_Prime_Return","is_return_saccade",...
                                                            {'Primary Fix','Return'}, {'rip_fix_distance_edge',"is_return_saccade","subject"},[1,1,250], dist_titl, "All",box_title,pspath);
    [curr_est4, curr_tStat4, curr_p4, nam4, pred4, anname4] = func_binary_testing(T_fix_average((T_fix_average.is_before_restricted | T_fix_average.is_after_restricted) & fix_overlap,:),...
                                                        "%s ~is_return_saccade +  rip_fix_distance_edge +(1|subject)", gr_names, (2:3),"Fix_Prime_Return_Restricted","is_return_saccade",...
                                                        {'Primary Fix','Return'}, {'rip_fix_distance_edge',"is_return_saccade","subject"},[1,1,20], dist_titl, "Restriced",box_title,pspath);
    [curr_est5, curr_tStat5, curr_p5, nam5, pred5, anname5] = func_binary_testing(T_fix_average((T_fix_average.is_before | T_fix_average.is_before) & fix_overlap,:),...
                                                            "%s ~is_return_light +  rip_fix_distance_edge +(1|subject)", gr_names, (2:3),"Fix_Prime_Return_lite","is_return_light",...
                                                            {'Primary Fix','Return Light'}, {'rip_fix_distance_edge',"is_return_light","subject"},[1,1,250], dist_titl, "All",box_title,pspath);
    [curr_est6, curr_tStat6, curr_p6, nam6, pred6, anname6] = func_binary_testing(T_fix_average((T_fix_average.is_before_restricted | T_fix_average.is_after_restricted) & fix_overlap,:),...
                                                        "%s ~is_return_light +  rip_fix_distance_edge +(1|subject)", gr_names, (2:3),"Fix_Prime_Return_lite_Restricted","is_return_light",...
                                                        {'Primary Fix','Return Light'}, {'rip_fix_distance_edge',"is_return_light","subject"},[1,1,20], dist_titl, "Restriced",box_title,pspath);
    [curr_est7, curr_tStat7, curr_p7, nam7, pred7, anname7] = func_binary_testing(T_fix_average((T_fix_average.is_before | T_fix_average.is_before) & fix_overlap,:),...
                                                            "%s ~is_refix_exact +  rip_fix_distance_edge +(1|subject)", gr_names, (2:3),"Fix_Prime_refix_exact","is_refix_exact",...
                                                            {'Primary Fix','refix exact'}, {'rip_fix_distance_edge',"is_refix_exact","subject"},[1,1,250], dist_titl, "All",box_title,pspath);
    [curr_est8, curr_tStat8, curr_p8, nam8, pred8, anname8] = func_binary_testing(T_fix_average((T_fix_average.is_before_restricted | T_fix_average.is_after_restricted) & fix_overlap,:),...
                                                        "%s ~is_refix_exact +  rip_fix_distance_edge +(1|subject)", gr_names, (2:3),"Fix_Prime_refix_exact_Restricted","is_refix_exact",...
                                                        {'Primary Fix','refix exact'}, {'rip_fix_distance_edge',"is_refix_exact","subject"},[1,1,20], dist_titl, "Restriced",box_title,pspath);
    [curr_est9, curr_tStat9, curr_p9, nam9, pred9, anname9] = func_binary_testing(T_fix_average((T_fix_average.is_before | T_fix_average.is_before) & fix_overlap,:),...
                                                            "%s ~is_valid_refix +  rip_fix_distance_edge +(1|subject)", gr_names, (2:3),"Fix_Prime_is_valid_refix","is_valid_refix",...
                                                            {'Primary Fix','valid refix'}, {'rip_fix_distance_edge',"is_valid_refix","subject"},[1,1,250], dist_titl, "All",box_title,pspath);
    [curr_est10, curr_tStat10, curr_p10, nam10, pred10, anname10] = func_binary_testing(T_fix_average((T_fix_average.is_before_restricted | T_fix_average.is_after_restricted) & fix_overlap,:),...
                                                        "%s ~is_valid_refix +  rip_fix_distance_edge +(1|subject)", gr_names, (2:3),"Fix_Prime_is_valid_refix_Restricted","is_valid_refix",...
                                                        {'Primary Fix','valid refix'}, {'rip_fix_distance_edge',"is_valid_refix","subject"},[1,1,20], dist_titl, "Restriced",box_title,pspath);

    curr_est = [curr_est1; curr_est2;curr_est3;curr_est4;curr_est5;curr_est6; curr_est7; curr_est8; curr_est9; curr_est10]; 
    curr_tStat = [curr_tStat1; curr_tStat2;curr_tStat3;curr_tStat4;curr_tStat5;curr_tStat6; curr_tStat7; curr_tStat8; curr_tStat9; curr_tStat10]; curr_p = [curr_p1; curr_p2;curr_p3;curr_p4;curr_p5;curr_p6; curr_p7; curr_p8;curr_p9;curr_p10]; 
    nam = [nam1, nam2,nam3,nam4,nam5,nam6, nam7, nam8, nam9, nam10]; pred = [pred1; pred2;pred3;pred4;pred5;pred6; pred7; pred8; pred9; pred10]; anname = [anname1; anname2;anname3;anname4;anname5;anname6; anname7; anname8; anname9;anname10];
   
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(pspath,".ps", ".pdf "),pspath));    delete(pspath); 
    
    T = table(anname,nam',pred,curr_est,curr_tStat,curr_p, ...
    'VariableNames', {'anName','Dependent','Predictor','Coef','tStat','pValue'});
    
    fn = fullfile(saveFolder,params.exact_analysisName,char(sprintf('LME_%s_Restr_Fixations_%s.xlsx',params.exact_analysisName,string(params.max_time)))); sheet = matlab.lang.makeValidName('fixations_prime_refix');
    writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
    ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
    sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
    wb.Save(); wb.Close(false); ex.Quit(); delete(ex);
    
    
    %%
    
    fix_overlap  = T_fix_tall.seg_to_ripple_distance_mid_abs >0;
    
    T_fix_tall.is_return_light = T_fix_tall.precursor_lag > 1;
    T_fix_tall.is_refix_exact = T_fix_tall.local_lag == 1;

    gr_names = {'R','R_abs','R_fisher','R_fisher_abs'}; dist_titl = "Absolute Rip-Fix distance (ms)"; box_title= "Fixation  Primary vs Return Ripple"; 
    pspath =  fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_SEG_vs_refix.ps',params.exact_analysisName));
    
    [curr_est1, curr_tStat1, curr_p1, nam1, pred1, anname1] = func_binary_testing(T_fix_tall((T_fix_tall.is_before | T_fix_tall.is_after) & fix_overlap,:),...
                                                            "%s ~is_refix + segment+ seg_to_ripple_distance_mid_abs +(1|subject)" , gr_names, (2:4),"Fix_Prime_Refix", "is_refix",...
                                                            {'Primary Fix','Refix'},{"seg_to_ripple_distance_mid_abs","is_refix","segment","subject"},[1,0,250], dist_titl, "All",box_title, pspath);
    [curr_est2, curr_tStat2, curr_p2, nam2, pred2, anname2] = func_binary_testing(T_fix_tall((T_fix_tall.is_before_restricted | T_fix_tall.is_after_restricted) & fix_overlap,:),...
                                                            "%s ~is_refix + segment+ seg_to_ripple_distance_mid_abs +(1|subject)", gr_names, (2:4),"Fix_Prime_Refix_Restricted","is_refix",...
                                                            {'Primary Fix','Refix'}, {'seg_to_ripple_distance_mid_abs',"is_refix","segment","subject"},[1,1,20], dist_titl, "Restriced",box_title,pspath);
    [curr_est3, curr_tStat3, curr_p3, nam3, pred3, anname3] = func_binary_testing(T_fix_tall((T_fix_tall.is_before | T_fix_tall.is_before) & fix_overlap,:),...
                                                            "%s ~is_return_saccade + segment+ seg_to_ripple_distance_mid_abs +(1|subject)", gr_names, (2:4),"Fix_Prime_Return","is_return_saccade",...
                                                            {'Primary Fix','Return'}, {'seg_to_ripple_distance_mid_abs',"is_return_saccade","segment","subject"},[1,1,250], dist_titl, "All",box_title,pspath);
    [curr_est4, curr_tStat4, curr_p4, nam4, pred4, anname4] = func_binary_testing(T_fix_tall((T_fix_tall.is_before_restricted | T_fix_tall.is_after_restricted) & fix_overlap,:),...
                                                        "%s ~is_return_saccade + segment+ seg_to_ripple_distance_mid_abs +(1|subject)", gr_names, (2:4),"Fix_Prime_Return_Restricted","is_return_saccade",...
                                                        {'Primary Fix','Return'}, {'seg_to_ripple_distance_mid_abs',"is_return_saccade","segment","subject"},[1,1,20], dist_titl, "Restriced",box_title,pspath);
    [curr_est5, curr_tStat5, curr_p5, nam5, pred5, anname5] = func_binary_testing(T_fix_tall((T_fix_tall.is_before | T_fix_tall.is_before) & fix_overlap,:),...
                                                            "%s ~is_return_light+ segment+ seg_to_ripple_distance_mid_abs +(1|subject)", gr_names, (2:4),"Fix_Prime_Return_lite","is_return_light",...
                                                            {'Primary Fix','Return Light'}, {'seg_to_ripple_distance_mid_abs',"is_return_light","segment","subject"},[1,1,250], dist_titl, "All",box_title,pspath);
    [curr_est6, curr_tStat6, curr_p6, nam6, pred6, anname6] = func_binary_testing(T_fix_tall((T_fix_tall.is_before_restricted | T_fix_tall.is_after_restricted) & fix_overlap,:),...
                                                        "%s ~is_return_light + segment+ seg_to_ripple_distance_mid_abs +(1|subject)", gr_names, (2:4),"Fix_Prime_Return_lite_Restricted","is_return_light",...
                                                        {'Primary Fix','Return Light'}, {'seg_to_ripple_distance_mid_abs',"is_return_light","segment","subject"},[1,1,20], dist_titl, "Restriced",box_title,pspath);
    [curr_est7, curr_tStat7, curr_p7, nam7, pred7, anname7] = func_binary_testing(T_fix_tall((T_fix_tall.is_before | T_fix_tall.is_before) & fix_overlap,:),...
                                                            "%s ~is_refix_exact + segment+ seg_to_ripple_distance_mid_abs +(1|subject)", gr_names, (2:4),"Fix_Prime_refix_exact","is_refix_exact",...
                                                            {'Primary Fix','Refix Exact'}, {'seg_to_ripple_distance_mid_abs',"is_refix_exact","segment","subject"},[1,1,250], dist_titl, "All",box_title,pspath);
    [curr_est8, curr_tStat8, curr_p8, nam8, pred8, anname8] = func_binary_testing(T_fix_tall((T_fix_tall.is_before_restricted | T_fix_tall.is_after_restricted) & fix_overlap,:),...
                                                        "%s ~is_refix_exact + segment+ seg_to_ripple_distance_mid_abs +(1|subject)", gr_names, (2:4),"Fix_Prime_refix_exact_Restricted","is_refix_exact",...
                                                        {'Primary Fix','Refix Exact'}, {'seg_to_ripple_distance_mid_abs',"is_refix_exact","segment","subject"},[1,1,20], dist_titl, "Restriced",box_title,pspath);
    [curr_est9, curr_tStat9, curr_p9, nam9, pred9, anname9] = func_binary_testing(T_fix_tall((T_fix_tall.is_before | T_fix_tall.is_before) & fix_overlap,:),...
                                                            "%s ~is_valid_refix + segment+ seg_to_ripple_distance_mid_abs +(1|subject)", gr_names, (2:4),"Fix_Prime_valid_refix","is_valid_refix",...
                                                            {'Primary Fix','Refix Valid'}, {'seg_to_ripple_distance_mid_abs',"is_valid_refix","segment","subject"},[1,1,250], dist_titl, "All",box_title,pspath);
    [curr_est10, curr_tStat10, curr_p10, nam10, pred10, anname10] = func_binary_testing(T_fix_tall((T_fix_tall.is_before_restricted | T_fix_tall.is_after_restricted) & fix_overlap,:),...
                                                        "%s ~is_valid_refix + segment+ seg_to_ripple_distance_mid_abs +(1|subject)", gr_names, (2:4),"Fix_Prime_valid_refix_Restricted","is_valid_refix",...
                                                        {'Primary Fix','Valid'}, {'seg_to_ripple_distance_mid_abs',"is_valid_refix","segment","subject"},[1,1,20], dist_titl, "Restriced",box_title,pspath);

    curr_est = [curr_est1; curr_est2;curr_est3;curr_est4;curr_est5;curr_est6; curr_est7; curr_est8; curr_est9; curr_est10]; 
    curr_tStat = [curr_tStat1; curr_tStat2;curr_tStat3;curr_tStat4;curr_tStat5;curr_tStat6; curr_tStat7; curr_tStat8; curr_tStat9; curr_tStat10]; curr_p = [curr_p1; curr_p2;curr_p3;curr_p4;curr_p5;curr_p6; curr_p7; curr_p8;curr_p9;curr_p10]; 
    nam = [nam1, nam2,nam3,nam4,nam5,nam6, nam7, nam8, nam9, nam10]; pred = [pred1; pred2;pred3;pred4;pred5;pred6; pred7; pred8; pred9; pred10]; anname = [anname1; anname2;anname3;anname4;anname5;anname6; anname7; anname8; anname9;anname10];
      
    
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(pspath,".ps", ".pdf "),pspath));    delete(pspath); 
    
    T = table(anname,nam',pred,curr_est,curr_tStat,curr_p, ...
    'VariableNames', {'anName','Dependent','Predictor','Coef','tStat','pValue'});
    
    fn = fullfile(saveFolder,params.exact_analysisName,char(sprintf('LME_%s_Restr_Fixations_%s.xlsx',params.exact_analysisName,string(params.max_time)))); sheet = matlab.lang.makeValidName('fixations_prime_refix_seg');
    writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
    ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
    sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
    wb.Save(); wb.Close(false); ex.Quit(); delete(ex);
    
    
    
    
    
    %% 2 Does Ripple to Fixation simialrity differs between Immediate (REFIX) vs Later Return? 
    
    temp_d = T_fix_average(logical(T_fix_average.is_refix),:);
    temp_d.is_immidieate_refixation = temp_d.precursor_lag == 1 |  temp_d.precursor_lag == 2;
    lme = fitlme(temp_d,sprintf("%s ~precursor_lag + precursor_distance +rip_fix_distance_edge +(1|subject)","R_fisher_abs"));
    
    
    fix_overlap  = temp_d.rip_fix_timing_middle_abs >params.fix_rip_distance_min;
    
    gr_names = {'R','R_abs','R_fisher','R_fisher_abs'}; dist_titl = "Absolute Rip-Fix distance (ms)"; box_title= "Refix Later vs Immediate Ripple"; 
    pspath =  fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_immediate_vs_refix.ps',params.exact_analysisName));
    
    [curr_est1, curr_tStat1, curr_p1, nam1, pred1, anname1] = func_binary_testing(temp_d((temp_d.is_before | temp_d.is_after) & fix_overlap,:),...
                                                            "%s ~is_immidieate_refixation +  rip_fix_distance_edge +(1|subject)" , gr_names, (2:3),"Refix_Immediate", "is_immidieate_refixation",...
                                                            {'Later Refix','Immediate'},{'rip_fix_distance_edge',"is_immidieate_refixation","subject"},[1,0,250], dist_titl, "All",box_title, pspath);
    [curr_est2, curr_tStat2, curr_p2, nam2, pred2, anname2] = func_binary_testing(temp_d((temp_d.is_before_restricted | temp_d.is_after_restricted) & fix_overlap,:),...
                                                            "%s ~is_immidieate_refixation +  rip_fix_distance_edge +(1|subject)", gr_names, (2:3),"Refix_Immediate_Restricted","is_immidieate_refixation",...
                                                            {'Later Refix','Immediate'}, {'rip_fix_distance_edge',"is_immidieate_refixation","subject"},[1,1,20], dist_titl, "Restriced",box_title,pspath);
    curr_est = [curr_est1; curr_est2]; curr_tStat = [curr_tStat1; curr_tStat2]; curr_p = [curr_p1; curr_p2]; nam = [nam1, nam2]; pred = [pred1; pred2];anname = [anname1; anname2];
    
    
    psDir =fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_immediate_vs_refix.ps',params.exact_analysisName));
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
    
    T = table(anname,nam',pred,curr_est,curr_tStat,curr_p, ...
    'VariableNames', {'anName','Dependent','Predictor','Coef','tStat','pValue'});
    
    fn = fullfile(saveFolder,params.exact_analysisName,char(sprintf('LME_%s_Restr_Fixations_%s.xlsx',params.exact_analysisName,string(params.max_time)))); sheet = matlab.lang.makeValidName('fixations_immediate_refix');
    writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
    ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
    sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
    wb.Save(); wb.Close(false); ex.Quit(); delete(ex);
    
    
    
%% 3 Are return Saccades tied to ripple emergence? 
    
% Analysis of Time to nearest Return Saccade 

saccades_table_refix = renamevars(saccades_table_refix,{'trial_number'},{'trial'});
saccades_table_refix.subject = saccades_table_refix.subject_new;

saccades_table_refix.is_return_light = saccades_table_refix.precursor_lag > 1;
saccades_table_refix.is_refix_exact = saccades_table_refix.local_lag == 1;
sac_refix_ripple = outerjoin(ripplewise_table_s,saccades_table_refix,'keys',{'subject','trial','indx_sac'},'MergeKeys', true, 'Type', 'left');

mask_restr = sac_refix_ripple.occular_duration > 0 & sac_refix_ripple.occular_duration < 1;



curr_est = []; curr_tStat = []; curr_p = []; nam = []; pred = []; anname = []; coef_range = (2:3);
gr_names = {"is_refix","is_return_saccade","is_return_light", "is_refix_exact", "is_valid_refix"};

for i = 1:length(gr_names)
    glme_all = fitglme(sac_refix_ripple(~isnan(sac_refix_ripple.is_refix),:), ...
        sprintf("is_saccades ~ %s +  indx_sac + (1|subject)",gr_names{i}), ...
        "Distribution", "binomial", "Link", "logit");
    
    glme_restr = fitglme(sac_refix_ripple(~isnan(sac_refix_ripple.is_refix) & mask_restr,:), ...
        sprintf("is_saccades ~ %s +  indx_sac + (1|subject)",gr_names{i}), ...
        "Distribution", "binomial", "Link", "logit");

    curr_est   = [curr_est;   glme_all.Coefficients.Estimate(coef_range); glme_restr.Coefficients.Estimate(coef_range);];
    curr_tStat = [curr_tStat; glme_all.Coefficients.tStat(coef_range); glme_restr.Coefficients.tStat(coef_range);];
    curr_p     = [curr_p;     glme_all.Coefficients.pValue(coef_range); glme_restr.Coefficients.pValue(coef_range);];
    pred       = [pred;       string(glme_all.Coefficients.Name(coef_range)); string(glme_restr.Coefficients.Name(coef_range));];
    nam        = [nam,        repmat("is_saccades", 1,  numel(coef_range)), repmat("is_saccades", 1, numel(coef_range))];
    anname     = [anname;     repmat("Return_Saccades",  numel(coef_range), 1); repmat("Restricted_Return_Saccades",  numel(coef_range), 1)];

end

T = table(anname,nam',pred,curr_est,curr_tStat,curr_p, ...
'VariableNames', {'anName','Dependent','Predictor','Coef','tStat','pValue'});

fn = fullfile(saveFolder,params.exact_analysisName,char(sprintf('LME_%s_Restr_Fixations_%s.xlsx',params.exact_analysisName,string(params.max_time)))); sheet = matlab.lang.makeValidName('return_saccades_ripples');
writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
wb.Save(); wb.Close(false); ex.Quit(); delete(ex);


%%
saccades_table_refix.ripple_is_in_sac = saccades_table_refix.rippleNumIn_saccade ~=0;



curr_est = []; curr_tStat = []; curr_p = []; nam = []; pred = []; anname = []; coef_range = (2:3);
gr_names = {"is_refix","is_return_saccade","is_return_light", "is_refix_exact" , "is_valid_refix"};

for i = 1:length(gr_names)
     glme_all = fitglme(saccades_table_refix(~isnan(sac_refix_ripple.is_refix),:), ...
    sprintf("ripple_is_in_sac ~ %s +  amplitude + (1|subject)",gr_names{i}), ...
    "Distribution", "binomial", "Link", "logit");

    
    glme_restr = fitglme(saccades_table_refix(~isnan(sac_refix_ripple.is_refix) & mask_restr,:), ...
    sprintf("ripple_is_in_sac ~ %s +  amplitude + (1|subject)",gr_names{i}), ...
    "Distribution", "binomial", "Link", "logit");
    curr_est   = [curr_est;   glme_all.Coefficients.Estimate(coef_range); glme_restr.Coefficients.Estimate(coef_range);];
    curr_tStat = [curr_tStat; glme_all.Coefficients.tStat(coef_range); glme_restr.Coefficients.tStat(coef_range);];
    curr_p     = [curr_p;     glme_all.Coefficients.pValue(coef_range); glme_restr.Coefficients.pValue(coef_range);];
    pred       = [pred;       string(glme_all.Coefficients.Name(coef_range)); string(glme_restr.Coefficients.Name(coef_range));];
    nam        = [nam,        repmat("is_saccades", 1,  numel(coef_range)), repmat("is_saccades", 1, numel(coef_range))];
    anname     = [anname;     repmat("Return_Saccades",  numel(coef_range), 1); repmat("Restricted_Return_Saccades",  numel(coef_range), 1)];

end

T = table(anname,nam',pred,curr_est,curr_tStat,curr_p, ...
'VariableNames', {'anName','Dependent','Predictor','Coef','tStat','pValue'});

fn = fullfile(saveFolder,params.exact_analysisName,char(sprintf('LME_%s_Restr_Fixations_%s.xlsx',params.exact_analysisName,string(params.max_time)))); sheet = matlab.lang.makeValidName('return_saccades_ripples_better');
writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
wb.Save(); wb.Close(false); ex.Quit(); delete(ex);


%% 4 Are there more ripples in Return vs Refixation vs First fixations?

fixation_table_refix = saccades_table_refix; %(:,{"name","subject"}[1,39,4,5,10,11,40:42,46:51,59,74:82]);
%fixation_table_refix= renamevars(fixation_table_refix,{'subject_new'},{'subject'});
fixation_table_refix.ripples_in_fixation_len = fixation_table_refix.ripples_in_fixation./fixation_table_refix.fixation_len;

[gs, subj_u]  = findgroups(ripplewise_table_s.subject);
channel_num_u = splitapply(@(c) numel(unique(c)), ripplewise_table_s.channel, gs);
channel_tbl   = table(subj_u, channel_num_u, 'VariableNames', {'subject','channel_num'});
fixation_table_refix = outerjoin(fixation_table_refix, channel_tbl, 'Keys', 'subject', 'MergeKeys', true, 'Type', 'left');
fixation_table_refix.ripples_in_fixation_chan = fixation_table_refix.ripples_in_fixation./fixation_table_refix.channel_num;
fixation_table_refix.is_return_light = fixation_table_refix.precursor_lag > 1;
fixation_table_refix.ripples_in_fixation_bin = fixation_table_refix.ripples_in_fixation ~=0;
fixation_table_refix = fixation_table_refix(~isnan(fixation_table_refix.fixation_len) & fixation_table_refix.fixation_len < 1 & fixation_table_refix.fixation_len > 0 ,:);
fixation_table_refix.is_refix_exact = fixation_table_refix.local_lag == 1;

gr_names = {'ripples_in_fixation_bin','ripples_in_fixation','ripples_in_fixation_len','ripples_in_fixation_chan'};
lme_par = {'binom','poisson','',''};

psDir =fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_in_fix_Refix.ps',params.exact_analysisName));

curr_est = []; curr_tStat = []; curr_p = []; nam = []; pred = []; anname = [];

   
[curr_est1, curr_tStat1, curr_p1, nam1, pred1, anname1] = func_binary_testing2(fixation_table_refix,...
                                                        "%s ~ is_refix +  fixation_len+ indx_sac +(1|subject) +(1|name)" , gr_names, (2:4),"Refix", lme_par, "is_refix",...
                                                        {'Later Refix','Immediate'},{'fixation_len',"is_refix","indx_sac","subject","name"},[1,1,0.1],  "Inital Fix vs Refix: is refix", "All", "Inital Fix vs Refix: is refix", psDir);
[curr_est2, curr_tStat2, curr_p2, nam2, pred2, anname2] = func_binary_testing2(fixation_table_refix,...
                                                        "%s ~ is_return_saccade +  fixation_len+ indx_sac +(1|subject) +(1|name)" , gr_names, (2:4),"Return", lme_par,"is_return_saccade",...
                                                        {'Later Refix','Immediate'},{'fixation_len',"is_return_saccade","indx_sac","subject","name"},[1,1,0.1],  "Inital Fix vs Refix: is return", "All", "Inital Fix vs Refix: is return", psDir);
[curr_est3, curr_tStat3, curr_p3, nam3, pred3, anname3] = func_binary_testing2(fixation_table_refix,...
                                                        "%s ~ is_return_light +  fixation_len+ indx_sac +(1|subject) +(1|name)" , gr_names, (2:4),"Return_light", lme_par, "is_return_light",...
                                                        {'Later Refix','Immediate'},{'fixation_len',"is_return_light","indx_sac","subject","name"},[1,1,0.1], "Inital Fix vs Refix: is return light", "All", "Inital Fix vs Refix: is  return light", psDir);
[curr_est4, curr_tStat4, curr_p4, nam4, pred4, anname4] = func_binary_testing2(fixation_table_refix,...
                                                        "%s ~ is_refix_exact +  fixation_len+ indx_sac +(1|subject) +(1|name)" , gr_names, (2:4),"Return_light", lme_par, "is_refix_exact",...
                                                        {'Later Refix','Immediate'},{'fixation_len',"is_refix_exact","indx_sac","subject","name"},[1,1,0.1], "Inital Fix vs Refix: is refix exect", "All", "Inital Fix vs Refix: is  refix exect", psDir);
[curr_est5, curr_tStat5, curr_p5, nam5, pred5, anname5] = func_binary_testing2(fixation_table_refix,...
                                                        "%s ~ is_valid_refix +  fixation_len+ indx_sac +(1|subject) +(1|name)" , gr_names, (2:4),"Return_light", lme_par, "is_valid_refix",...
                                                        {'Later Refix','Immediate'},{'fixation_len',"is_valid_refix","indx_sac","subject","name"},[1,1,0.1], "Inital Fix vs Refix: is valid refix", "All", "Inital Fix vs Refix: is  valid refix", psDir);


curr_est = [curr_est1; curr_est2;curr_est3;curr_est4;curr_est5]; curr_tStat = [curr_tStat1; curr_tStat2;curr_tStat3;curr_tStat4;curr_tStat5]; nam = [nam1, nam2,nam3,nam4,nam5];
curr_p = [curr_p1; curr_p2;curr_p3;curr_p4;curr_p5];  pred = [pred1; pred2;pred3;pred4;pred5];anname = [anname1; anname2;anname3;anname4;anname5];


system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
T = table(anname,nam',pred,curr_est,curr_tStat,curr_p, ...
'VariableNames', {'anName','Dependent','Predictor','Coef','tStat','pValue'});

fn = fullfile(saveFolder,params.exact_analysisName,char(sprintf('LME_%s_Restr_Fixations_%s.xlsx',params.exact_analysisName,string(params.max_time)))); sheet = matlab.lang.makeValidName('ripples_in_fix_refix');
writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
wb.Save(); wb.Close(false); ex.Quit(); delete(ex);

    
    
