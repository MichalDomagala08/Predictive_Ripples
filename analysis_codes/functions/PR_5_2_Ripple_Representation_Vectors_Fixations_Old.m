
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
params.exact_analysisName = "RSA";
params.make_repr = 0;               % to make representaitons actually 
params.make_rsa  = 0;               % if we are making RSA tables or using exsisting ones 
params.data_preprocessing = 1;      % if the data needs preprocessing or if we load preprocessed data directly

%%% Parameters of Representation
% bha parameters
params.bha          = 0;            % If we are to add bha signal to our representations.
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
saccadeFile = 'D:\Documents_Dell\Predictive_Ripples_2025\data\sacextr_ekm\content_all_subjects.mat';

data = load(fullfile(saveDataFolder,"rippleDensity.mat"));
saccades_table = data.saccades_table;

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
    repr = load(fullfile(saveFolder,reprTitle));
    subject_level_repr_rp = repr.subject_level_repr_rp;
    subject_level_repr_sac = repr.subject_level_repr_sac;
    subject_level_repr_fix = repr.subject_level_repr_fix;

end

%%% REPRESENTATIONAL SIMILATRITY: (1) Saccase X Ripple; (2) Ripple x Ripple; (3) Saccade x Saccade (4) Fixation X Ripple
if params.make_rsa
    [T_all,T_rip,T_sac,T_fix] = func_create_representational_similarity_Rs(subject_level_repr_rp,subject_level_repr_sac,subject_level_repr_fix,saccades_table,ripple_table,ripple_data,saveFolder,params);
else
    if params.bha; reprTitle = "RSA_tables_bha.mat"; else; reprTitle = "RSA_tables.mat"; end

    rsa = load(fullfile(saveFolder,reprTitle));
    T_all = rsa.T_all;
    T_fix = rsa.T_fix;
    T_rip = rsa.T_rip;
    T_sac = rsa.T_sac;
end


%%% GET FIXATION VARIABLES 
[T_fix_cleaned, T_fix_average, T_fix_tall]   = func_RSA_preprocessng_fixation(T_fix,params);


%%% GET SACCADE VARIABLES 
[sac_overlap,T_all_cleaned,Tall_cleaned,R_table] = func_RSA_preprocessing_saccades(T_all,ripple_table,params);

%%% GET CLOSEST RIPPLE TO SACCADES
R_unified= func_RSA_preprocessing_nearest_sac_fix(T_fix,T_all,ripple_table,params);

%%% GET CLOSEST RIPPLE TO SACCADES (segments)
R_unified_segments = func_RSA_preprocessing_nearest_sac_fix_segment(T_fix_tall,T_all,ripple_table,params);

%%% GET EQUALLY SPACED  (it takes time...)
R_unified_reverse  = func_create_representational_vectors_lite(R_unified,ripple_data,ripple_table,dataFolder,params);
R_unified_reverse.R_reverse_abs = abs(R_unified_reverse.R_reverse); 
R_unified_reverse.R_reverse_fisher = atanh(min(max(R_unified_reverse.R_reverse,-0.9999),0.9999)); 
R_unified_reverse(isnan(R_unified_reverse.R_reverse),:).R_reverse_fisher = nan(sum(isnan(R_unified_reverse.R_reverse)),1);
R_unified_reverse.R_reverse_fisher_abs = abs(R_unified_reverse.R_reverse_fisher);

% 1) Get Long format R_unified reverse for lienar model testing. 
has = ~isnan(R_unified_reverse.R_reverse);
T_fwd = R_unified_reverse; T_fwd.reversed = false(height(T_fwd),1); T_fwd.R = R_unified_reverse.R;
T_rev = R_unified_reverse(has,:); T_rev.reversed = true(height(T_rev),1); T_rev.R = R_unified_reverse.R_reverse(has);
R_unified_long = [T_fwd; T_rev];
R_unified_long.R_abs = abs(R_unified_long.R); R_unified_long.R_fisher = atanh(min(max(R_unified_long.R,-0.9999),0.9999)); R_unified_long.R_fisher_abs = abs(R_unified_long.R_fisher);


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



%% ======= ANALYSIS ======= %%%

%% %%%  ----- ANALYSIS CLOSEST FIXATION vs OPPOSITE SPACE ------- %%%

% (1) Closest Fixation  vs Opposite the same length
R_unified_reverse_fix  = R_unified_long(R_unified_long.is_fix & ~isnan(R_unified_long.R) & ~isnan(R_unified_long.R_reverse) & ...
                                        R_unified_long.is_abs_closest & R_unified_long.rip_event_distance_edge >0,:);


curr_est = []; curr_tStat = []; curr_p = []; nam = []; pred = [];anname = [];

gr_names = {'R','R_abs','R_fisher','R_fisher_abs'};
for i = 1:4
    
    lme = fitlme(R_unified_reverse_fix,sprintf("%s ~ reversed +  event_rip_latency_abs +(1|subject) +(1|rip_num)", string(gr_names{i})));
    curr_est = [curr_est;lme.Coefficients.Estimate(2:3)]; curr_tStat = [curr_tStat;lme.Coefficients.tStat(2:3)]; curr_p = [curr_p;lme.Coefficients.pValue(2:3)];
    pred = [pred; string(lme.Coefficients.Name(2:3))]; nam = [nam, string(gr_names{i}),string(gr_names{i})]; anname = [anname; "RipProx"; "RipProx"];

    figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
    sgtitle(sprintf("All: %s: P: %2.3f T: %2.3f",string(strrep(gr_names{i},'_','-')),lme.Coefficients.pValue(3),lme.Coefficients.tStat(3)))
    subplot(1,2,1)
    boxplot(R_unified_reverse_fix.(string(gr_names{i})), R_unified_reverse_fix.reversed, 'Labels', {'Fixation','Not Fixation'}); % upewnij się kolejność etykiet
    ylabel(string(strrep(gr_names{i},'_','-'))); title(sprintf("Restricted Fixation vs Reverse %s",string(strrep(gr_names{i},'_','-'))));
    title("Fixation vs Reverse to Ripple RSA");
    
    subplot(1,2,2)
    plot_linear_model_bin2(R_unified_reverse_fix,lme, string(gr_names{i}),'reversed',{'event_rip_latency_abs',"reversed","subject","rip_num"},1,1,250)
    xticks ([sort([xticks 25])])
    xticklabels(string(str2double(xticklabels) * 2));
    xlabel("Absolute Saccade to Ripple distance (ms)");    
    title("Closeness linear Model: ");
    
    print(gcf,fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_reverse.ps',params.exact_analysisName)),'-dpsc','-append','-fillpage');
    close(gcf);


    
    lme = fitlme(R_unified_reverse_fix(R_unified_reverse_fix.restricted,:),sprintf("%s ~  reversed +  event_rip_latency_abs +(1|subject) +(1|rip_num)", string(gr_names{i})));
    curr_est = [curr_est;lme.Coefficients.Estimate(2:3)];  curr_tStat = [curr_tStat;lme.Coefficients.tStat(2:3)]; curr_p = [curr_p;lme.Coefficients.pValue(2:3)];
    pred     = [pred; string(lme.Coefficients.Name(2:3))]; nam  = [nam, string(gr_names{i}),string(gr_names{i})]; anname = [anname; "RipProx_Restr"; "RipProx_Restr"];

    figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
    sgtitle(sprintf("Restricted: %s: P: %2.3f T: %2.3f",string(strrep(gr_names{i},'_','-')),lme.Coefficients.pValue(3),lme.Coefficients.tStat(3)))
    subplot(1,2,1)
    boxplot(R_unified_reverse_fix(R_unified_reverse_fix.restricted ,:).(string(gr_names{i})), ...
            R_unified_reverse_fix(R_unified_reverse_fix.restricted,:).reversed, 'Labels', {'Fixation','Not Fixation'}); % upewnij się kolejność etykiet
    ylabel(string(strrep(gr_names{i},'_','-'))); title(sprintf("Restricted Fixation vs Reverse %s",string(strrep(gr_names{i},'_','-'))));
    title("Fixation vs Reverse to Ripple RSA");
    
    subplot(1,2,2)
    plot_linear_model_bin2(R_unified_reverse_fix(R_unified_reverse_fix.restricted,:),lme, string(gr_names{i}),'reversed',{'event_rip_latency_abs',"reversed","subject","rip_num"},1,1,20)
    xticks ([sort([xticks 25])])
    xticklabels(string(str2double(xticklabels) * 2));
    xlabel("Absolute Fixation to Ripple distance (ms)");    
    title("Closeness linear Model: ");

    print(gcf,fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_reverse.ps',params.exact_analysisName)),'-dpsc','-append','-fillpage');
    close(gcf);


end


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


curr_est = []; curr_tStat = []; curr_p = []; nam = []; pred = [];anname = [];

gr_names = {'R','R_abs','R_fisher','R_fisher_abs'};
for i = 1:4
    
    lme = fitlme(R_unified_reverse_fix,sprintf("%s ~ reversed +  event_rip_latency_abs +(1|subject) +(1|rip_num)", string(gr_names{i})));
    curr_est = [curr_est;lme.Coefficients.Estimate(2:3)]; curr_tStat = [curr_tStat;lme.Coefficients.tStat(2:3)]; curr_p = [curr_p;lme.Coefficients.pValue(2:3)];
    pred = [pred; string(lme.Coefficients.Name(2:3))]; nam = [nam, string(gr_names{i}),string(gr_names{i})]; anname = [anname; "RipProx"; "RipProx"];

    figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
    sgtitle(sprintf("All: %s: P: %2.3f T: %2.3f",string(strrep(gr_names{i},'_','-')),lme.Coefficients.pValue(3),lme.Coefficients.tStat(3)))
    subplot(1,2,1)
    boxplot(R_unified_reverse_fix.(string(gr_names{i})), R_unified_reverse_fix.reversed, 'Labels', {'Fixation','Not Fixation'}); % upewnij się kolejność etykiet
    ylabel(string(strrep(gr_names{i},'_','-'))); title(sprintf("Restricted Fixation vs Reverse %s",string(strrep(gr_names{i},'_','-'))));
    title("Fixation vs Reverse to Ripple RSA");
    
    subplot(1,2,2)
    plot_linear_model_bin2(R_unified_reverse_fix,lme, string(gr_names{i}),'reversed',{'event_rip_latency_abs',"reversed","subject","rip_num"},1,1,250)
    xticks ([sort([xticks 25])])
    xticklabels(string(str2double(xticklabels) * 2));
    xlabel("Absolute Saccade to Ripple distance (ms)");    
    title("Closeness linear Model: ");
    
    print(gcf,fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_sac_reverse.ps',params.exact_analysisName)),'-dpsc','-append','-fillpage');
    close(gcf);


    
    lme = fitlme(R_unified_reverse_fix(R_unified_reverse_fix.restricted,:),sprintf("%s ~  reversed +  event_rip_latency_abs +(1|subject) +(1|rip_num)", string(gr_names{i})));
    curr_est = [curr_est;lme.Coefficients.Estimate(2:3)];  curr_tStat = [curr_tStat;lme.Coefficients.tStat(2:3)]; curr_p = [curr_p;lme.Coefficients.pValue(2:3)];
    pred     = [pred; string(lme.Coefficients.Name(2:3))]; nam  = [nam, string(gr_names{i}),string(gr_names{i})]; anname = [anname; "RipProx_Restr"; "RipProx_Restr"];

    figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
    sgtitle(sprintf("Restricted: %s: P: %2.3f T: %2.3f",string(strrep(gr_names{i},'_','-')),lme.Coefficients.pValue(3),lme.Coefficients.tStat(3)))
    subplot(1,2,1)
    boxplot(R_unified_reverse_fix(R_unified_reverse_fix.restricted ,:).(string(gr_names{i})), ...
            R_unified_reverse_fix(R_unified_reverse_fix.restricted,:).reversed, 'Labels', {'Fixation','Not Fixation'}); % upewnij się kolejność etykiet
    ylabel(string(strrep(gr_names{i},'_','-'))); title(sprintf("Restricted Fixation vs Reverse %s",string(strrep(gr_names{i},'_','-'))));
    title("Fixation vs Reverse to Ripple RSA");
    
    subplot(1,2,2)
    plot_linear_model_bin2(R_unified_reverse_fix(R_unified_reverse_fix.restricted,:),lme, string(gr_names{i}),'reversed',{'event_rip_latency_abs',"reversed","subject","rip_num"},1,1,20)
    xticks ([sort([xticks 25])])
    xticklabels(string(str2double(xticklabels) * 2));
    xlabel("Absolute Fixation to Ripple distance (ms)");    
    title("Closeness linear Model: ");

    print(gcf,fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_sac_reverse.ps',params.exact_analysisName)),'-dpsc','-append','-fillpage');
    close(gcf);


end


psDir =fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_sac_reverse.ps',params.exact_analysisName));
system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 

T = table(anname,nam',pred,curr_est,curr_tStat,curr_p, ...
'VariableNames', {'anName','Dependent','Predictor','Coef','tStat','pValue'});

fn = fullfile(saveFolder,params.exact_analysisName,char(sprintf('LME_%s_Restr_Fixations_%s.xlsx',params.exact_analysisName,string(params.max_time)))); sheet = matlab.lang.makeValidName('sac_reverse_sim');
writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
wb.Save(); wb.Close(false); ex.Quit(); delete(ex);



%% %%%  ----- ANALYSIS CLOSEST FIXATION vs  CLOSEST SACCAES ------- %%%

% (1) Closeste Fixation and Closest Saccade
curr_est = []; curr_tStat = []; curr_p = []; nam = []; pred = [];anname = [];
gr_names = {'R','R_abs','R_fisher','R_fisher_abs'};
for i = 1:4
    
    lme = fitlme(R_unified(R_unified.is_abs_closest & R_unified.rip_event_distance_edge >0,:),sprintf("%s ~ is_fix +  event_rip_latency_abs +(1|subject)", string(gr_names{i})));
    curr_est = [curr_est;lme.Coefficients.Estimate(2:3)]; curr_tStat = [curr_tStat;lme.Coefficients.tStat(2:3)]; curr_p = [curr_p;lme.Coefficients.pValue(2:3)];
    pred = [pred; string(lme.Coefficients.Name(2:3))]; nam = [nam, string(gr_names{i}),string(gr_names{i})]; anname = [anname; "RipProx"; "RipProx"];

    figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
    sgtitle(sprintf("All: %s: P: %2.3f T: %2.3f",string(strrep(gr_names{i},'_','-')),lme.Coefficients.pValue(2),lme.Coefficients.tStat(2)))
    subplot(1,2,1)
    boxplot(R_unified(R_unified.is_abs_closest & R_unified.rip_event_distance_edge >0,:).(string(gr_names{i})),...
            R_unified(R_unified.is_abs_closest & R_unified.rip_event_distance_edge >0,:).is_fix, 'Labels', {'Saccade','Fixation'}); % upewnij się kolejność etykiet
    ylabel(string(strrep(gr_names{i},'_','-'))); title(sprintf("Restricted Fixation vs Saccades %s",string(strrep(gr_names{i},'_','-'))));
    title("Fixation vs Saccades to Ripple RSA");
    
    subplot(1,2,2)
    plot_linear_model_bin2(R_unified(R_unified.is_abs_closest  & R_unified.rip_event_distance_edge >0,:),lme, string(gr_names{i}),'is_fix',{'event_rip_latency_abs',"is_fix","subject"},1,1,250)
    xticks ([sort([xticks 25])])
    xticklabels(string(str2double(xticklabels) * 2));
    xlabel("Absolute Saccade to Ripple distance (ms)");    
    title("Closeness linear Model: ");
    
    print(gcf,fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_vs_sac.ps',params.exact_analysisName)),'-dpsc','-append','-fillpage');
    close(gcf);


    
    lme = fitlme(R_unified(R_unified.is_abs_closest & R_unified.restricted & R_unified.rip_event_distance_edge >0,:),sprintf("%s ~ is_fix +  event_rip_latency_abs +(1|subject)", string(gr_names{i})));
    curr_est = [curr_est;lme.Coefficients.Estimate(2:3)]; curr_tStat = [curr_tStat;lme.Coefficients.tStat(2:3)]; curr_p = [curr_p;lme.Coefficients.pValue(2:3)];
    pred = [pred; string(lme.Coefficients.Name(2:3))]; nam = [nam, string(gr_names{i}),string(gr_names{i})]; anname = [anname; "RipProx_Restr"; "RipProx_Restr"];

    figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
    sgtitle(sprintf("Restricted: %s: P: %2.3f T: %2.3f",string(strrep(gr_names{i},'_','-')),lme.Coefficients.pValue(2),lme.Coefficients.tStat(2)))
    subplot(1,2,1)
    boxplot(R_unified(R_unified.is_abs_closest & R_unified.restricted & R_unified.rip_event_distance_edge >0,:).(string(gr_names{i})), ...
            R_unified(R_unified.is_abs_closest & R_unified.restricted & R_unified.rip_event_distance_edge >0,:).is_fix, 'Labels', {'Saccade','Fixation'}); % upewnij się kolejność etykiet
    ylabel(string(strrep(gr_names{i},'_','-'))); title(sprintf("Restricted Fixation vs Saccades %s",string(strrep(gr_names{i},'_','-'))));
    title("Fixation vs Saccades to Ripple RSA");
    
    subplot(1,2,2)
    plot_linear_model_bin2(R_unified(R_unified.is_abs_closest & R_unified.restricted & R_unified.rip_event_distance_edge >0,:),lme, string(gr_names{i}),'is_fix',{'event_rip_latency_abs',"is_fix","subject"},1,1,20)
    xticks ([sort([xticks 25])])
    xticklabels(string(str2double(xticklabels) * 2));
    xlabel("Absolute Saccade to Ripple distance (ms)");    
    title("Closeness linear Model: ");

    print(gcf,fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_vs_sac.ps',params.exact_analysisName)),'-dpsc','-append','-fillpage');
    close(gcf);


end


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
curr_est = []; curr_tStat = []; curr_p = []; nam = []; pred = [];anname = [];
gr_names = {'R','R_abs','R_fisher','R_fisher_abs'};
for i = 1:4
    
    lme = fitlme(R_unified_segments(R_unified_segments.is_abs_closest & R_unified_segments.rip_event_distance_edge >0,:),sprintf("%s ~ is_fixseg + event_idx+  event_rip_latency_abs +(1|subject)", string(gr_names{i})));
    curr_est = [curr_est;lme.Coefficients.Estimate(2:4)]; curr_tStat = [curr_tStat;lme.Coefficients.tStat(2:4)]; curr_p = [curr_p;lme.Coefficients.pValue(2:4)];
    pred = [pred; string(lme.Coefficients.Name(2:4))]; nam = [nam, string(gr_names{i}),string(gr_names{i}),string(gr_names{i})]; anname = [anname; "RipProx"; "RipProx"; "RipProx"];


    figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
    sgtitle(sprintf("All: %s: P: %2.3f T: %2.3f",string(strrep(gr_names{i},'_','-')),lme.Coefficients.pValue(2),lme.Coefficients.tStat(2)))
    subplot(1,2,1)
    plot_linear_model_bin2(R_unified_segments(R_unified_segments.is_abs_closest & R_unified_segments.rip_event_distance_edge >0,:),lme, string(gr_names{i}),'is_fixseg',{'event_idx',"is_fixseg","event_rip_latency_abs","subject"},1,0,5)
    xlabel("Saccade Segment");    
    title(sprintf("Segment linear Model: P: %2.3f | T: %2.3f",lme.Coefficients.pValue(3),lme.Coefficients.tStat(3)));

    subplot(1,2,2)
    plot_linear_model_bin2(R_unified_segments(R_unified_segments.is_abs_closest & R_unified_segments.rip_event_distance_edge >0,:),lme, string(gr_names{i}),'is_fixseg',{'event_rip_latency_abs',"is_fixseg","event_idx","subject"},1,0,250)
    xticks ([sort([xticks 25])])
    xticklabels(string(str2double(xticklabels) * 2));
    xlabel("Absolute Saccade to Ripple distance (ms)");    
    title(sprintf("Closeness linear Model: P: %2.3f | T: %2.3f",lme.Coefficients.pValue(4),lme.Coefficients.tStat(4)));

    print(gcf,fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_vs_sac_seg.ps',params.exact_analysisName)),'-dpsc','-append','-fillpage');
    close(gcf);

    % For restricted (-/+ params.max ms around Onset)
    lme = fitlme(R_unified_segments(R_unified_segments.is_abs_closest & R_unified_segments.rip_event_distance_edge >0 & R_unified_segments.restricted,:),sprintf("%s ~ is_fixseg+ event_idx +  event_rip_latency_abs +(1|subject) ", string(gr_names{i})));
    curr_est = [curr_est;lme.Coefficients.Estimate(2:4)]; curr_tStat = [curr_tStat;lme.Coefficients.tStat(2:4)]; curr_p = [curr_p;lme.Coefficients.pValue(2:4)];
    pred = [pred; string(lme.Coefficients.Name(2:4))]; nam = [nam, string(gr_names{i}),string(gr_names{i}),string(gr_names{i})]; anname = [anname; "RipProx_Restr"; "RipProx_Restr" ; "RipProx_Restr"];

    figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
    sgtitle(sprintf("Restricted (+/- %d ms): %s: P: %2.3f T: %2.3f",params.max_time/2,string(strrep(gr_names{i},'_','-')),lme.Coefficients.pValue(2),lme.Coefficients.tStat(2)))
    subplot(1,2,1)
    plot_linear_model_bin2(R_unified_segments(R_unified_segments.is_abs_closest & R_unified_segments.rip_event_distance_edge >0 & R_unified_segments.restricted,:),...
                            lme, string(gr_names{i}),'is_fixseg',{'event_idx',"is_fixseg","event_rip_latency_abs","subject"},1,0,5)
    xlabel("Segment");    
    title(sprintf("Segment linear Model: P: %2.3f | T: %2.3f",lme.Coefficients.pValue(2),lme.Coefficients.tStat(2)));


    subplot(1,2,2)
    plot_linear_model_bin2(R_unified_segments(R_unified_segments.is_abs_closest & R_unified_segments.rip_event_distance_edge >0 & R_unified_segments.restricted,:),...
                        lme, string(gr_names{i}),'is_fixseg',{'event_rip_latency_abs',"is_fixseg","event_idx","subject"},1,1,20)
    xticks ([sort([xticks 25])])
    xticklabels(string(str2double(xticklabels) * 2));
    xlabel("Absolute Saccade to Ripple distance (ms)");    
    title(sprintf("Closeness linear Model: P: %2.3f | T: %2.3f",lme.Coefficients.pValue(4),lme.Coefficients.tStat(4)));

    print(gcf,fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_vs_sac_seg.ps',params.exact_analysisName)),'-dpsc','-append','-fillpage');
    close(gcf);



end


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
% 
% for i = 1:4
% 
%     %%% For all Before vs After 
%     lme = fitlme(T_fix_average((T_fix_average.is_before | T_fix_average.is_after) & fix_overlap,:),sprintf("%s ~ is_after +  rip_fix_distance_edge +(1|subject)", string(gr_names{i})));
%     curr_est = [curr_est;lme.Coefficients.Estimate(2:3)]; curr_tStat = [curr_tStat;lme.Coefficients.tStat(2:3)]; curr_p = [curr_p;lme.Coefficients.pValue(2:3)];
%     pred = [pred; string(lme.Coefficients.Name(2:3))]; nam = [nam, string(gr_names{i}),string(gr_names{i})]; anname = [anname; "RipProx_Restr"; "RipProx_Restr"];
% 
% 
%     figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
%     sgtitle(sprintf("All: %s: P: %2.3f T: %2.3f",string(strrep(gr_names{i},'_','-')),lme.Coefficients.pValue(2),lme.Coefficients.tStat(2)))
%     subplot(1,2,1)
%     boxplot(T_fix_average((T_fix_average.is_before | T_fix_average.is_after)  & fix_overlap,:).(string(gr_names{i})), T_fix_average((T_fix_average.is_before | T_fix_average.is_after)  & fix_overlap,:).is_after, 'Labels', {'after','before'}); % upewnij się kolejność etykiet
%     ylabel(string(strrep(gr_names{i},'_','-'))); title(sprintf("Restricted Before vs After %s",string(strrep(gr_names{i},'_','-'))));
%     title("After - Before Boxplot");
% 
%     subplot(1,2,2)
%     plot_linear_model_bin2(T_fix_average( (T_fix_average.is_after | T_fix_average.is_before)  & fix_overlap,:),lme, string(gr_names{i}),'is_after',{'rip_fix_distance_edge',"is_after","subject"},1,0,250)
%     xticks ([sort([xticks 25])])
%     xticklabels(string(str2double(xticklabels) * 2));
%     xlabel("Absolute Saccade to Ripple distance (ms)");    
%     title("Closeness linear Model: ");
% 
%     print(gcf,fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_aft_vs_bef.ps',params.exact_analysisName)),'-dpsc','-append','-fillpage');
%     close(gcf);
% 
%     % For restricted (-/+ params.max ms around Onset)
%     lme = fitlme(T_fix_average((T_fix_average.is_before_restricted | T_fix_average.is_after_restricted )  & fix_overlap ,:),sprintf("%s ~ is_after_restricted +  rip_fix_distance_edge +(1|subject) ", string(gr_names{i})));
%     curr_est = [curr_est;lme.Coefficients.Estimate(2:3)]; curr_tStat = [curr_tStat;lme.Coefficients.tStat(2:3)]; curr_p = [curr_p;lme.Coefficients.pValue(2:3)];
%     pred = [pred; string(lme.Coefficients.Name(2:3))]; nam = [nam, string(gr_names{i}),string(gr_names{i})]; anname = [anname; "RipProx_Restr"; "RipProx_Restr"];
% 
%     figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
%     sgtitle(sprintf("Restricted (+/- %d ms): %s: P: %2.3f T: %2.3f",params.max_time/2,string(strrep(gr_names{i},'_','-')),lme.Coefficients.pValue(2),lme.Coefficients.tStat(2)))
%     subplot(1,2,1)
%     boxplot(T_fix_average((T_fix_average.is_before_restricted | T_fix_average.is_after_restricted )  & fix_overlap ,:).(string(gr_names{i})), T_fix_average((T_fix_average.is_before_restricted | T_fix_average.is_after_restricted)  & fix_overlap ,:).is_after_restricted, 'Labels', {'after','before'}); % upewnij się kolejność etykiet
%     ylabel(string(strrep(gr_names{i},'_','-'))); title(sprintf("Restricted Before vs After %s",string(strrep(gr_names{i},'_','-'))));
%     title("After - Before Boxplot");
% 
%     subplot(1,2,2)
%     plot_linear_model_bin2(T_fix_average((T_fix_average.is_after_restricted | T_fix_average.is_before_restricted)  & fix_overlap,:),lme, string(gr_names{i}),'is_after_restricted',{'rip_fix_distance_edge',"is_after_restricted","subject"},1,1,20)
%     xticks ([sort([xticks 25])])
%     xticklabels(string(str2double(xticklabels) * 2));
%     xlabel("Absolute Saccade to Ripple distance (ms)");    
%     title("Closeness linear Model: ");
% 
%     print(gcf,fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_aft_vs_bef.ps',params.exact_analysisName)),'-dpsc','-append','-fillpage');
%     close(gcf);
% 
% end


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
% 
% 
% 
% for i = 1:4
% 
%     %%% For all Before vs After 
%     lme = fitlme(T_fix_tall((T_fix_tall.is_before | T_fix_tall.is_after) & fix_overlap,:),sprintf("%s ~ is_after+ segment+ seg_to_ripple_distance_mid_abs +(1|subject)", string(gr_names{i})));
%     curr_est = [curr_est;lme.Coefficients.Estimate(2:4)]; curr_tStat = [curr_tStat;lme.Coefficients.tStat(2:4)]; curr_p = [curr_p;lme.Coefficients.pValue(2:4)];
%     pred = [pred; string(lme.Coefficients.Name(2:4))]; nam = [nam, string(gr_names{i}),string(gr_names{i}),string(gr_names{i})]; anname = [anname; "RipProx"; "RipProx"; "RipProx"];
% 
% 
%     figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
%     sgtitle(sprintf("All: %s: P: %2.3f T: %2.3f",string(strrep(gr_names{i},'_','-')),lme.Coefficients.pValue(2),lme.Coefficients.tStat(2)))
%     subplot(1,2,1)
%     plot_linear_model_bin2(T_fix_tall( T_fix_tall.is_after | T_fix_tall.is_before  & fix_overlap,:),lme, string(gr_names{i}),'is_after',{'segment',"is_after","seg_to_ripple_distance_mid_abs","subject"},1,0,5)
%     xlabel("Saccade Segment");    
%     title(sprintf("Segment linear Model: P: %2.3f | T: %2.3f",lme.Coefficients.pValue(2),lme.Coefficients.tStat(2)));
% 
%     subplot(1,2,2)
%     plot_linear_model_bin2(T_fix_tall( T_fix_tall.is_after | T_fix_tall.is_before  & fix_overlap,:),lme, string(gr_names{i}),'is_after',{'seg_to_ripple_distance_mid_abs',"is_after","segment","subject"},1,0,250)
%     xticks ([sort([xticks 25])])
%     xticklabels(string(str2double(xticklabels) * 2));
%     xlabel("Absolute Saccade to Ripple distance (ms)");    
%     title(sprintf("Closeness linear Model: P: %2.3f | T: %2.3f",lme.Coefficients.pValue(4),lme.Coefficients.tStat(4)));
% 
%     print(gcf,fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_aft_vs_bef_seg.ps',params.exact_analysisName)),'-dpsc','-append','-fillpage');
%     close(gcf);
% 
%     % For restricted (-/+ params.max ms around Onset)
%     lme = fitlme(T_fix_tall((T_fix_tall.is_before_restricted | T_fix_tall.is_after_restricted )  & fix_overlap ,:),sprintf("%s ~ is_after_restricted+ segment +  seg_to_ripple_distance_mid_abs +(1|subject) ", string(gr_names{i})));
%     curr_est = [curr_est;lme.Coefficients.Estimate(2:4)]; curr_tStat = [curr_tStat;lme.Coefficients.tStat(2:4)]; curr_p = [curr_p;lme.Coefficients.pValue(2:4)];
%     pred = [pred; string(lme.Coefficients.Name(2:4))]; nam = [nam, string(gr_names{i}),string(gr_names{i}),string(gr_names{i})]; anname = [anname; "RipProx_Restr"; "RipProx_Restr" ; "RipProx_Restr"];
% 
%     figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
%     sgtitle(sprintf("Restricted (+/- %d ms): %s: P: %2.3f T: %2.3f",params.max_time/2,string(strrep(gr_names{i},'_','-')),lme.Coefficients.pValue(2),lme.Coefficients.tStat(2)))
%     subplot(1,2,1)
%     plot_linear_model_bin2(T_fix_tall( T_fix_tall.is_after_restricted | T_fix_tall.is_before_restricted  & fix_overlap,:),lme, string(gr_names{i}),'is_after',{'segment',"is_after_restricted","seg_to_ripple_distance_mid_abs","subject"},1,0,5)
%     xlabel("Segment");    
%     title(sprintf("Segment linear Model: P: %2.3f | T: %2.3f",lme.Coefficients.pValue(2),lme.Coefficients.tStat(2)));
% 
% 
%     subplot(1,2,2)
%     plot_linear_model_bin2(T_fix_tall((T_fix_tall.is_after_restricted | T_fix_tall.is_before_restricted)  & fix_overlap,:),lme, string(gr_names{i}),'is_after_restricted',{'seg_to_ripple_distance_mid_abs',"is_after_restricted","segment","subject"},1,0,20)
%     xticks ([sort([xticks 25])])
%     xticklabels(string(str2double(xticklabels) * 2));
%     xlabel("Absolute Saccade to Ripple distance (ms)");    
%     title(sprintf("Closeness linear Model: P: %2.3f | T: %2.3f",lme.Coefficients.pValue(4),lme.Coefficients.tStat(4)));
% 
%     print(gcf,fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_aft_vs_bef_seg.ps',params.exact_analysisName)),'-dpsc','-append','-fillpage');
%     close(gcf);
% 
% end


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
nam = [nam1, nam2]; pred = [pred1, pred2]; anname = [anname1; anname2];



% 
% 
% curr_est = []; curr_tStat = []; curr_p = []; nam = []; pred = [];anname = [];
% gr_names = {'R','R_abs','R_fisher','R_fisher_abs'};
% fix_overlap  = T_fix_tall.seg_to_ripple_distance_mid_abs >0;
% T_fix_tall.fix_phase_cat = categorical(T_fix_tall.fix_phase);
% testinf = {{1,2},{2,3},{3,4},{4,5}};
% 
% 
% for i = 1:4
% 
%     %%% For all Before vs After 
%     lme = fitlme(T_fix_tall((T_fix_tall.is_before | T_fix_tall.is_after) & fix_overlap,:),sprintf("%s ~ fix_phase +  seg_to_ripple_distance_mid_abs +(1|subject)", string(gr_names{i})));
%     curr_est = [curr_est;lme.Coefficients.Estimate(2:3)]; curr_tStat = [curr_tStat;lme.Coefficients.tStat(2:3)]; curr_p = [curr_p;lme.Coefficients.pValue(2:3)];
%     pred = [pred; string(lme.Coefficients.Name(2:3))]; nam = [nam, string(gr_names{i}),string(gr_names{i})]; anname = [anname; "RipProx"; "RipProx"];
% 
%     % Progression model: 
% 
%     figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
%     sgtitle(sprintf("All: %s: P: %2.3f T: %2.3f",string(strrep(gr_names{i},'_','-')),lme.Coefficients.pValue(2),lme.Coefficients.tStat(2)))
%     subplot(1,2,1);     hold on;
%     boxplot(T_fix_tall((T_fix_tall.is_before | T_fix_tall.is_after)  & fix_overlap,:).(string(gr_names{i})), T_fix_tall((T_fix_tall.is_before | T_fix_tall.is_after)  & fix_overlap,:).fix_phase, 'Labels', {'early2','early1','middle','late1','late2'}); % upewnij się kolejność etykiet
%     for j = 1:4
% 
%         aaaa = T_fix_tall((T_fix_tall.is_before | T_fix_tall.is_after) & fix_overlap & (T_fix_tall.fix_phase == testinf{j}{1} | T_fix_tall.fix_phase == testinf{j}{2}),:);
%         lme_t = fitlme(aaaa,sprintf("%s ~ fix_phase +  seg_to_ripple_distance_mid_abs +(1|subject)", string(gr_names{i})));
% 
%         if lme_t.Coefficients.pValue(2) < 0.05; star = ' *'; else; star = ''; end
%         y = max(T_fix_tall((T_fix_tall.is_before | T_fix_tall.is_after)  & fix_overlap,:).(string(gr_names{i})),[],'all') - 0.05;
%         plot([j j j+1 j+1],[y y+0.02 y+0.02 y],'k','LineWidth',1.2)
%         text(mean([j j+1]), y+0.05, sprintf('p=%.3g%s',lme_t.Coefficients.pValue(2), star), 'HorizontalAlignment','center')
%     end
%     ylabel(string(strrep(gr_names{i},'_','-'))); title(sprintf("After - Before Boxplot %s",string(strrep(gr_names{i},'_','-'))));
% 
%     cmap = lines(5);
%     hBox = flipud(findobj(gca,'Tag','Box')); hMed = flipud(findobj(gca,'Tag','Median')); hOut = flipud(findobj(gca,'Tag','Outliers'));
%     for k=1:min(5,numel(hBox)), set(hBox(k),'Color',cmap(k,:)); set(hMed(k),'Color',max(0,cmap(k,:)-0.15)); if k<=numel(hOut), set(hOut(k),'MarkerEdgeColor',cmap(k,:)); end; end
% 
%     subplot(1,2,2)
%     plot_linear_model_multilevel(T_fix_tall( T_fix_tall.is_after | T_fix_tall.is_before  & fix_overlap,:),lme, string(gr_names{i}),'fix_phase',{'seg_to_ripple_distance_mid_abs',"fix_phase","subject"},1,0,250)
%     xticks ([sort([xticks 25])])
%     xticklabels(string(str2double(xticklabels) * 2));
%     xlabel("Absolute Saccade to Ripple distance (ms)");    
%     title(sprintf("%s Closeness:  P: %2.3f T: %2.3f",string(strrep(gr_names{i},'_','-')),lme.Coefficients.pValue(2),lme.Coefficients.tStat(2)));
% 
%     lme_c = fitlme(T_fix_tall((T_fix_tall.is_before | T_fix_tall.is_after) & fix_overlap,:),sprintf("%s ~ fix_phase_cat +  seg_to_ripple_distance_mid_abs +(1|subject)", string(gr_names{i})));
%     c_idx = find(contains(lme_c.Coefficients.Name,"fix_phase_cat"));
%     pvals_length5 = lme_c.Coefficients.pValue(c_idx);
%     lg = legend(gca); s = cellstr(lg.String); n=numel(s);
%     for k=1:4, pv=pvals_length5(k); if ~isnan(pv), s{n-4+k}=sprintf('%s (p=%.3g)', s{n-4+k}, pv); end; end
%     lg.String = s;
% 
%     print(gcf,fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_phase.ps',params.exact_analysisName)),'-dpsc','-append','-fillpage');
%     close(gcf);
% 
% 
%     %%% Restricted!
%     lme = fitlme(T_fix_tall((T_fix_tall.is_before | T_fix_tall.is_after) & fix_overlap & T_fix_tall.seg_to_ripple_distance_mid_abs  <250,:),sprintf("%s ~ fix_phase +  seg_to_ripple_distance_mid_abs +(1|subject)", string(gr_names{i})));
%     curr_est = [curr_est;lme.Coefficients.Estimate(2:3)]; curr_tStat = [curr_tStat;lme.Coefficients.tStat(2:3)]; curr_p = [curr_p;lme.Coefficients.pValue(2:3)];
%     pred = [pred; string(lme.Coefficients.Name(2:3))]; nam = [nam, string(gr_names{i}),string(gr_names{i})]; anname = [anname; "RipProx_Restr"; "RipProx_Restr"];
% 
%     % Progression model: 
% 
%     figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
%     sgtitle(sprintf("Restricted: %s: P: %2.3f T: %2.3f",string(strrep(gr_names{i},'_','-')),lme.Coefficients.pValue(2),lme.Coefficients.tStat(2)))
%     subplot(1,2,1);     hold on; curr_p2 = [];
%     boxplot(T_fix_tall((T_fix_tall.is_before | T_fix_tall.is_after)  & fix_overlap,:).(string(gr_names{i})), T_fix_tall((T_fix_tall.is_before | T_fix_tall.is_after)  & fix_overlap,:).fix_phase, 'Labels', {'early2','early1','middle','late1','late2'}); % upewnij się kolejność etykiet
%     for j = 1:4
%         aaaa = T_fix_tall((T_fix_tall.is_before | T_fix_tall.is_after) & fix_overlap & T_fix_tall.seg_to_ripple_distance_mid_abs  <250 & (T_fix_tall.fix_phase == testinf{j}{1} | T_fix_tall.fix_phase == testinf{j}{2}),:);
%         lme_t = fitlme(aaaa,sprintf("%s ~ fix_phase +  seg_to_ripple_distance_mid_abs +(1|subject)", string(gr_names{i})));
%         if lme_t.Coefficients.pValue(2) < 0.05; star = ' *'; else; star = ''; end
%         y = max(T_fix_tall((T_fix_tall.is_before | T_fix_tall.is_after)  & fix_overlap,:).(string(gr_names{i})),[],'all') - 0.05;
%         plot([j j j+1 j+1],[y y+0.02 y+0.02 y],'k','LineWidth',1.2)
%         text(mean([j j+1]), y+0.05, sprintf('p=%.3g%s',lme_t.Coefficients.pValue(2), star), 'HorizontalAlignment','center')
%     end
% 
%     cmap = lines(5);
%     hBox = flipud(findobj(gca,'Tag','Box')); hMed = flipud(findobj(gca,'Tag','Median')); hOut = flipud(findobj(gca,'Tag','Outliers'));
%     for k=1:min(5,numel(hBox)), set(hBox(k),'Color',cmap(k,:)); set(hMed(k),'Color',max(0,cmap(k,:)-0.15)); if k<=numel(hOut), set(hOut(k),'MarkerEdgeColor',cmap(k,:)); end; end
% 
%     subplot(1,2,2)
%     plot_linear_model_multilevel(T_fix_tall(( T_fix_tall.is_after | T_fix_tall.is_before)  & fix_overlap & T_fix_tall.seg_to_ripple_distance_mid_abs <250,:),lme, string(gr_names{i}),'fix_phase',{'seg_to_ripple_distance_mid_abs',"fix_phase","subject"},1,0,20)
%     xticks ([sort([xticks 25])])
%     xticklabels(string(str2double(xticklabels) * 2));
%     xlabel("Absolute Saccade to Ripple distance (ms)");    
%     title(sprintf("Restricted %s Closeness:  P: %2.3f T: %2.3f",string(strrep(gr_names{i},'_','-')),lme.Coefficients.pValue(2),lme.Coefficients.tStat(2)));
% 
%     lme_c = fitlme(T_fix_tall((T_fix_tall.is_before | T_fix_tall.is_after) & fix_overlap & T_fix_tall.seg_to_ripple_distance_mid_abs  <250,:),sprintf("%s ~ fix_phase_cat +  seg_to_ripple_distance_mid_abs +(1|subject)", string(gr_names{i})));
%     c_idx = find(contains(lme_c.Coefficients.Name,"fix_phase_cat"));
%     pvals_length5 = lme_c.Coefficients.pValue(c_idx);
%     lg = legend(gca); s = cellstr(lg.String); n=numel(s);
%     for k=1:4, pv=pvals_length5(k); if ~isnan(pv), s{n-4+k}=sprintf('%s (p=%.3g)', s{n-4+k}, pv); end; end
%     lg.String = s;
% 
%     print(gcf,fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_phase.ps',params.exact_analysisName)),'-dpsc','-append','-fillpage');
%     close(gcf);
% 
% end



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


%%% (1) Ripples In saccades vs Other


gr_names = {'R','R_abs','R_fisher','R_fisher_abs'}; fromul = "%s ~ is_saccades +  event_rip_latency_abs +(1|subject) ";
pspath =  fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_in_saccades_fix_rsa.ps',params.exact_analysisName));
dist_titl = "Absolute Rip-Event distance (ms)"; box_title= "In vs Out Saccade Ripple to Fix RSA"; model_covars = {'event_rip_latency_abs',"is_saccades","subject"};

[curr_est1, curr_tStat1, curr_p1, nam1, pred1, anname1] = func_binary_testing(ripplewise_fix,fromul , gr_names, (2:3),"Ripples_In_Saccades",...
                                                        "is_saccades",{'Not in saccades','In saccades'},model_covars,[1,0,250], dist_titl, "All",box_title, pspath);
[curr_est2, curr_tStat2, curr_p2, nam2, pred2, anname2] = func_binary_testing(ripplewise_fix(ripplewise_fix.restricted,:), fromul, gr_names, (2:3),"Ripples_In_Saccades_Restricted",...
                                                        "is_saccades",{'Not in saccades','In saccades'}, model_covars,[1,1,20], dist_titl, "Restriced",box_title,pspath);
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

