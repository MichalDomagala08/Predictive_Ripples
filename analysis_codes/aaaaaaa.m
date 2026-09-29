function [classif_table_iters,ica2stat] = sc_5_ds_classify_meg(infile, ds, rc, chan, ica, seg, currentFolder, progressfile,Data_volume,  classification_threshold, Mode, Unsure_only, Validation, Verbosity, model_type)

    %%% Load Model info and ica stat
    classif_table_iters = {};
    if isfield(ica, 'ICA2_classModel')
        if ~strcmp(ica.ICA2_classModel(end-3:end), '.mat')
            ica.ICA2_classModel = strcat(ica.ICA2_classModel,'.mat');
        end
        load (ica.ICA2_classModel);
    else
        warning('Model for classification of components not specifiied. Trying to fit the existing one, check the classification results !!')
    end
    
    
        % load file

    [~, dsname, ~] = fileparts(infile{ds}); % File parts return path, name and extension, we want only name
    disp(strcat(' [', num2str( round(100*ds/length(infile) )), '%] Loading artifree ICA2 decomposed DS: ', dsname, ' - ', infile{ds}));
    load_file = infile{ds};
    % tmpdata = load(load_file,'data_ica2');
    % data_ica2 =  tmpdata.data_ica2;
    % clear tmpdata
    load(load_file);

    whole_g_ic = [];
    
    % It reads biosemi EEG system : 70 channels 64 eeg, 4 electric, 2 ref

    %Extracting data from all_ica for specific subject: ICA, events, frequency, ica_parameters
    ICA = data_ica2.ICA;
    event = data_ica2.event;
    Fsnew = data_ica2.fsample;
    ncond = size(event,2);
    ica_par = data_ica2.ica_par;
    dsname = data_ica2.hdr.dsname;
    saveFolderds = fullfile(rc.analysisFolder, '6_ICA2', dsname);
    classif_table_all.("ds_"+dsname) = struct();

    disp(sprintf(strcat('\nClassifying ICA components; DS: ', dsname)));


    %Getting electric channels info
    n_elec = length(chan.elec); %number of electric channels
    
    chan2icanr = [find(strcmp(data_ica2.hdr.actualchantype, 'eeg'))' find(strcmp(data_ica2.hdr.actualchantype, 'mag'))' find(strcmp(data_ica2.hdr.actualchantype, 'meg'))' find(strcmp(data_ica2.hdr.actualchantype, 'meggrad'))' find(strcmp(data_ica2.hdr.actualchantype, 'megmag'))'];
    chan2ica = data_ica2.label(chan2icanr);

    if rc.eeg
        elecloc_ica = sc_trim_sens(chan2ica, data_ica2.hdr.elec);
    else
        elecloc_ica = sc_trim_sens(chan2ica, data_ica2.hdr.grad);  %provides labels of chans that survived AR and were used for ICA
    end

    if strcmp(Data_volume, 'all')
        num_iters = size(ICA.iter,2);
    elseif strcmp(Data_volume,'best') % Only Best N features (5) 
        [~,I] = sort([ICA.iter.total_ic_number],'descend');
        I = I(1:5);
        num_iters = 5;
    else
        num_iters =1;
    end

    trim_folder = fullfile (rc.analysisFolder, '6_ICA2', dsname);

    for indeks = 1:num_iters

            if isfield(rc.simul, 'setallICasbrain')
                if rc.simul.setallICasbrain % When this parameter is ON, omit classification by model

                    if rc.simul.setallICasbrain   %settinga all brain as IC - for sim purposes only!!
                        for it =  length(data_ica2.ICA.iter)
                            nofic(it)= data_ica2.ICA.iter.total_ic_number;
                        end
                        [icnobestiter,bestindex] = max(nofic);
                        data_ica2.ICA.best_iter.index = bestindex;
                        data_ica2.ICA.iter(bestindex).brain_ic = 1:icnobestiter;
                        data_ica2.ICA.iter(bestindex).brain_ic_number = data_ica2.ICA.iter(bestindex).total_ic_number;
                    end

                    continue
                end
            end
        

        disp([' [' num2str( round(100*((ds-1)/length(infile) + indeks/num_iters/length(infile)))) '%] Classification progress: DS: ', dsname] );

        if strcmp(Data_volume, 'bestiter')
            indeks = ICA.best_iter.index;
        elseif strcmp(Data_volume,'best')
            indeks = I(indeks)
        end
        %
        fprintf('ICA iteration no. %d',indeks);
        disp('')
        disp('    IC identification ...');

        %Getting best channels for analysis as well as targets for validation
        [nIC, nptIC] = size(ICA.iter(indeks).IC_sig);
        %For those that there is brain_signal give One (in case that there exist previous classification.)

        if isempty(ICA.iter(indeks).unmixing)
            continue
        end


        %%----------------------------------------
        %------ Eye peaks detection Feature  -----
        % Finds the number of Peaks in data 
        clear peaks_per_min
        %%% MD -  Experimental signal Feature - Work in progress
        for comps = 1:size(ICA.iter(indeks).IC_sig,1)
            sig_temp =ICA.iter(indeks).IC_sig(comps,:);
            num_peaks_in_window = zeros(round(length(sig_temp)/(5*data_ica2.hdr.Fs)),1);
            for windows = 1:round(length(sig_temp)/(5*data_ica2.hdr.Fs))
                try
                    current_window = sig_temp((windows-1)*5*data_ica2.hdr.Fs +1:windows*5*data_ica2.hdr.Fs);
                catch
                    current_window = sig_temp((windows-1)*5*data_ica2.hdr.Fs +1:end);
                end
                [peaks,LOCS,width,prominence] = findpeaks(current_window,'MinPeakHeight',std(current_window)*5);
                [peaks2,LOCS2,width2,prominence2] = findpeaks(-current_window,'MinPeakHeight',std(current_window)*5);
                %                 [peaks,~,~,~] = findpeaks(abs(current_window),'MinPeakDistance',Fsnew*0.5,'MinPeakHeight',std(current_window)*2)   ;
                peaks = [peaks peaks2];
                warn = warning('query','last');
                if ~isempty(warn)
                    warn_id = warn.identifier;
                    warning('off',warn_id)
                end
                num_peaks_in_window(windows) = length(peaks);
            end

            peaks_per_min(comps) = sum(num_peaks_in_window)/(length(data_ica2.ica_par.time_ica)/(60*data_ica2.hdr.Fs));
            clear num_peaks_in_window
        end


        if (strcmp(Mode,'class') || strcmp(Mode,'review')) %||strcmp(Mode,'label'))


            %% %%%% FEATURE LOADING %%%%%%%%%%%

            %%%% LOADING IMAGE DATA AND TARGET DATA TO VECOTRS

            % variable initialization

            disp('    Extracting correlational variables...');
            % iter = ICA.iter(indeks); % ICA iteration
            temp_Y_vec = 1:ICA.iter(indeks).total_ic_number;

            Y = ismember(temp_Y_vec,ICA.iter(indeks).brain_ic);             % our target vector
            elc_signal_corr = ICA.iter(indeks).elc_signal_correlation;      % correlation with signal of electric electrodes
            elc_power_corr = ICA.iter(indeks).elc_power_correlation;        % correlation with signal of electric electrodes
            elc_spectrum_corr = ICA.iter(indeks).elc_spectrum_correlation;  % correlation with signal of electric electrodes
            spect_one_over_f = ICA.iter(indeks).spectrum_one_over_f;        % similarity to 1/f
            spect_flat = ICA.iter(indeks).spectrum_flat;                    % flat spectrum similarity
            kurtosis = ICA.iter(indeks).time_kurtosis;                      % kurtosis

            %Creating Correlational Data Table

            % For further Fixing -- FIXME - not doable by chan.elec
            if isfield(chan,"addMockEOG")
                eogModel = 0;
            elseif ~isempty(chan.elec) 
                eogModel = 1;
            else
                eogModel = 1;
            end
            if eogModel  % EOG present %%% FIXME - defininf NOT by chan.elec, but by used model...
            
                unmx_table = array2table(normalize(elc_signal_corr'),'VariableNames',{'ELC_signal_corr'});
                unmx_table.Elc_spectrum_corr = normalize(elc_spectrum_corr');
                unmx_table.Spect_one_over_f = normalize(spect_one_over_f');
                unmx_table.Spect_flat = normalize(spect_flat');
                unmx_table.Kurtosis = normalize(kurtosis');
                unmx_table.ELC_power_corr = normalize(elc_power_corr');

                %%% Currently Opared ERP Variability
                if strcmp(seg.studyType,'ER')
                    temp_erp_var = ICA.iter(indeks).T(1,:).^2;
                    for ts = 2:size(ICA.iter(indeks).T,1)
                        temp_erp_var =temp_erp_var + ICA.iter(indeks).T(ts,:).^2;
                    end
                    temp_erp_var = sqrt(temp_erp_var);
                    unmx_table.ERP_variance = normalize(temp_erp_var');
                end


                %%% New, Improved ERP Variability: 
                % Norm of ERP Variability across Events 
                % if strcmp(seg.studyType,'ER')
                %     temp_erp_var = iter.T(1,:).^2;
                %     for ts = 2:size(iter.T,1)
                %         temp_erp_var =temp_erp_var + iter.T(ts,:).^2;
                %     end
                %     temp_erp_var = sqrt(temp_erp_var);
                %     unmx_table.ERP_variance = normalize(temp_erp_var');
                % end

                if rc.eeg == 0
                    unmx_table.Eye_variable = normalize(peaks_per_min');
                end

            else  % no EOG available
                % FOR NOW! MD - FIXME!!
                unmx_table = array2table(normalize(elc_signal_corr'),'VariableNames',{'ELC_signal_corr'});
                unmx_table.Elc_spectrum_corr = normalize(elc_spectrum_corr');
                unmx_table.Spect_one_over_f = normalize(spect_one_over_f');
                unmx_table.Spect_flat = normalize(spect_flat');
                unmx_table.Kurtosis = normalize(kurtosis');
                unmx_table.ELC_power_corr = normalize(elc_power_corr');
                unmx_table.Eye_variable =  normalize(peaks_per_min');


                if strcmp(seg.studyType,'ER')
                    temp_erp_var = ICA.iter(indeks).T(1,:).^2;
                    for ts = 2:size(ICA.iter(indeks).T,1)
                        temp_erp_var =temp_erp_var + ICA.iter(indeks).T(ts,:).^2;
                    end
                    temp_erp_var = sqrt(temp_erp_var);
                    unmx_table.ERP_variance = (temp_erp_var');
                end
            
            end

            ims = zeros(263, 350,3,size(ICA.iter(indeks).time_kurtosis,2)); %3,nIC); %Dummy images list

            elc_corr = elc_signal_corr;

            clear  elc_signal_corr elc_power_corr elc_spectrum_corr spect_one_over_f kurtosis spect_flat e eeg_channels elec_channels
            disp('    Extracting image variables...');

            % Getting Images in form of Matrix
            % set(0,'DefaultFigureVisible','off');
            % fprintf('        Processing image no.: 0 ')

            disp(strrep(infile(ds), '_ica2.mat', ''))

            if rc.eeg
                elecloc_ica = sc_trim_sens(data_ica2.label(chan2icanr), data_ica2.hdr.elec);  %provides labels of chans that survived AR and were used for ICA
            else
                elecloc_ica = sc_trim_sens(data_ica2.label(chan2icanr), data_ica2.hdr.grad);  %provides labels of chans that survived AR and were used for ICA
            end
            if exist(saveFolderds)==0; mkdir (saveFolderds); end
            sc_locfile( elecloc_ica, saveFolderds); %trims elc file according to actual no of channels

            parfor k = 1:size(ICA.iter(indeks).mixing,2)  

                fig2 = figure();
                set(fig2,'visible','off', 'position',[0,0,263,250]);
                topoplot(ICA.iter(indeks).mixing(:,k), char(fullfile(trim_folder, 'ch_trimmed.elc')), 'nosedir', chan.ori,'verbose','off','style' , 'map','headrad' ,0, 'electrodes', 'off');
                F =imresize(getframe(fig2).cdata,[263, 350]);

                ims(:,:,:,k) = double(F);
                close(fig2)
                disp(k)

            end

            set(0,'DefaultFigureVisible','on');
            disp('    Data Extracion Complete!\n')

            clear fig z elec_channels eeg_channels k

            %% ------------------------------------------------------
   

            fprintf('Predicting using %s classifier with the %s model...\n',model_type,ica.ICA2_classModel)

            switch model_type

                case 'own' %Custom Architecture

                    sig_corr_test = arrayDatastore(table2array(unmx_table));

                    %splitting image data
                    imdsTest = arrayDatastore(ims,'IterationDimension',4);

                    %Splitting Target Values
                    dsTest = combine(imdsTest,sig_corr_test);
                    mbqTest =  minibatchqueue(dsTest,...
                        'MiniBatchSize',size(ICA.iter(indeks).mixing,2),...
                        'MiniBatchFcn',@sc_preprocessMiniBatch,...
                        'MiniBatchFormat',{'SSCB','CB'});

                    dlNet = OWN_model;

                    % Predirtion
                    [predictions,prob_values] = sc_jointDLPredictions(dlNet,mbqTest,[1 0]);

                    clear d sTest mbqTest dlnet sig_corr_test

                case 'pretrained' %Pretrained Hybrid Model


                    dlNet = Hybrid_model;
                    net = inceptionv3();


                    %getting input layer size for extraction
                    imageSize = net.Layers(1).InputSize;
                    augmentedTestSet = augmentedImageDatastore(imageSize, ims);

                    %feature layer name depends on pretrained architecture anf must be manually
                    featureLayer = net.Layers(end).Name;

                    % Extracting Training and Testing Features ,via Activations
                    testFeatures = activations(net, augmentedTestSet, featureLayer, ...
                        'MiniBatchSize', 22, 'OutputAs', 'columns');

                    clear augmentedTrainingSet augmentedTestSet featureLayer net imageSize imdsTrain imdsTest

                    XTest2 = array2table(testFeatures');
                    XTest2.Elc_sig_corr = unmx_table.ELC_signal_corr;
                    XTest2.Elc_spectrum_corr = unmx_table.Elc_spectrum_corr;
                    XTest2.Spect_one_over_f = unmx_table.Spect_one_over_f;
                    XTest2.Spect_flat = unmx_table.Spect_flat;
                    XTest2.Kurtosis = unmx_table.Kurtosis;
                    XTest2.Elc_power_corr = unmx_table.ELC_power_corr;

                    prob_values_1 = predict(dlNet,  XTest2);
                    prob_values = prob_values_1(:,2);

                    predictions = onehotdecode(prob_values_1,[0 1],2);
                    variable_comparison = XTest2;
                    variable_comparison.predictions = prob_values;

                    clear prob_values_1 XTest2 testFeatures net dlNet

                case 'if' %basic classifier based on two distinct models

              %%     
                   
                if rc.resetGPU && ~isempty(gpuDeviceTable)
                     try
                         reset(gpuDevice(1));
                     catch
                         parallel.gpu.enableCUDAForwardCompatibility(true)
                         reset(gpuDevice(1));
                     end
                end
                %%% Choose proper Image and Feature model ratio in
                %%% weighted probability 
                if ~isfield(rc, 'modelScaling')
                    switch ica.ICA2_classModel

                        %%% Munster and Ctf format
                        case "classification_models_MEG_ER_noelec.mat"
                            augmentedTestSet = augmentedImageDatastore([263 351 3], ims);

                            convScaling = 3;
                            featScaling = 2;

                        %%% Itab or Itab-Sim
                        case "classification_models_MEG_ER_elec.mat"
                            augmentedTestSet = augmentedImageDatastore([332 442 3], ims);

                            convScaling = 5;
                            featScaling = 2;

                        %%% Boston or other "normal" Data
                        case "classification_models_EEG_RS_elec.mat"
                            augmentedTestSet = augmentedImageDatastore([263 350 3], ims);

                            convScaling = 4;
                            featScaling = 2;
                    
                        case "classification_models_EEG_ER_elec.mat"
                            % Features In (5):
                            %  + elc_signal_corr
                            %  + Elc_spectrum_corr
                            %  + Spect_one_over_f
                            %  + Spect_flat
                            %  + Kurtosis
                            %  + Elc_power_corr

                            augmentedTestSet = augmentedImageDatastore([332 442 3], ims);

                            convScaling = 3;
                            featScaling = 2;

                        case "classification_models_EEG_RS_elec_retrained2.mat" 
                            augmentedTestSet = augmentedImageDatastore([263 350 3], ims);

                            convScaling = 3;
                            featScaling = 2;
                        % 
                        % case "classification_models_EEG_ER_mockEOG.mat"
                        %     % Features In (7):
                        %     %  + elc_signal_corr
                        %     %  + Elc_spectrum_corr
                        %     %  + Spect_one_over_f
                        %     %  + Spect_flat
                        %     %  + Kurtosis
                        %     %  + Elc_power_corr
                        %     %  + Eye_variable
                        %     %  + ERP_variance
                        %     augmentedTestSet = augmentedImageDatastore([332 442 3], ims);
                        % 
                        %     convScaling = 2;
                        %     featScaling = 2;
                        otherwise

                            
                            augmentedTestSet = augmentedImageDatastore([263 350 3], ims);

                            convScaling = 2;
                            featScaling = 2;
                    end
                else
                    convScaling = modelScaling(1);
                    featScaling = modelScaling(2);
                end
                %%% Predict Both Images and Features
                predictions_images = predict(IF_conv_model,augmentedTestSet,'MiniBatchSize',24);
                predictions_features = predict(IF_mlp_model,unmx_table);
                prob_values = (predictions_images(:,1)*convScaling+ featScaling*predictions_features(:,2))/(convScaling+featScaling); 
                
                % fix for highly ELC correlated componnets
                %prob_values(find (elc_corr > 0.66)) = 0;

                %%% Create visualization Structure For Classification
                classif_table = table(predictions_images(:,1), predictions_features(:,2),...
                    'VariableNames', {'Images','Features'});

                classif_table.elc_corr = elc_corr'; % Correlation with ELC channels
                classif_table.Is_Problematic = categorical(classif_table.Images - classif_table.Features > 0.3);
                classif_table.Proposed_Classification = double(prob_values > classification_threshold);
                classif_table.Y = Y';
                fieldNameTemp = char(sprintf("Iter_%s",string(indeks)));


                %%% Structure of All Classif Table Used to calculate
                %%% general accuracy

                %%% Predictions Used for all of the Script
                predictions = classif_table.Proposed_Classification;
                classif_table.prob_values = prob_values;
                classif_table_iters{end+1} = classif_table;
                clear  augmentedTestSet
                %%
            end
            disp('    Prediction extraction complete!\n')

        end
        % iter = ICA.iter(indeks);
        elc_signal_corr = ICA.iter(indeks).elc_signal_correlation;%(1:nIC); % correlation with signal of electric electrodes

        %% -----------------------------------------------
        %------ POWER SPECTRAL DENSITY ESTIMATION  -----

        clear mspettro
        win = 1024*ceil(4/(ica_par.dec_fact+1));
        for ix=1:nIC
            spectral_estimator = spectrum.welch('Hamming',win,0);
            Power_Spectrum = psd(spectral_estimator,ICA.iter(indeks).IC_sig(ix,:),'NFFT',win,'Fs',Fsnew);
            mspettro(ix,:) = sqrt((Fsnew/win)*Power_Spectrum.Data)';
            F = Power_Spectrum.Frequencies;
        end
        clear spectral_estimator Power_Spectrum ix
        window = 10;          %expressed in seconds
        window_base = find(ica_par.time_ica >= ica_par.time_ica(1) & ica_par.time_ica < window+ica_par.time_ica(1));     % expressed in number of points


        %% Manual Validation Loop + saving classified results

        % prob_values i s used for setting brain components

        %For now it asks almost every time for user input
        % Future iterations may also have Verbosity argument

        %         disp('    === Manual Validation ===')
        flag = 0;
        flag2 = 0;

        classif_new = [];

        not_sure_ic = [];
        ButtonName ='';

        ICA_after_threshold =[];
        while flag == 0

            iterator_IC = 1:nIC;
            i = 1;

            if ~isempty(not_sure_ic)
                iterator_IC = not_sure_ic;
                not_sure_ic = [];

            elseif Unsure_only==true
                try
                    b_ic = [];
                    g_ic = data_ica2.ICA.iter(indeks).brain_ic;
                    iterator_IC = ICA.iter(indeks).Unsure;
                    not_sure_ic2 = iterator_IC;
                    if isempty(iterator_IC)
                        disp('WARNING! Vector Unsure is empty, No need to consider unsure labels!')
                    end

                catch
                    disp("WARNING! no Unsure vector detected! Proceeding with normal labeling!")
                    b_ic = [];
                    g_ic = [];
                end
            else
                b_ic = [];
                g_ic = [];


            end

            while i <= length(iterator_IC)
                it = iterator_IC(i);

                sig = ICA.iter(indeks).IC_sig(it,:);

                %-------------------------------------------%
                %%% AVERAGING SIGNALS AND MOVING AVERAGES %%%

                if (Verbosity == 1 || strcmp(Mode,'review') || strcmp(Mode,'label')) && ~strcmp(seg.studyType , 'RS')

                    %Averaging Facotr
                    avg_factor = 12;

                    %Variables for average 2
                    averages_table = table();

                    for ind = 1:ncond
                        ave = average(sig,event(ind).peaks_pos,...
                            round(event(ind).t_pre*Fsnew/1000),...
                            round(event(ind).t_post*Fsnew/1000),...
                            round(event(ind).t_baseline_on*Fsnew/1000),...
                            round(event(ind).t_baseline_off*Fsnew/1000));
                        name = string(ind);
                        averages_table.(name) = ave;
                    end

                    %MOVING AVERAGE
                    %Adding Table columns together accoding to every entry on ave
                    Moving_Average_2 = zeros(length(averages_table.(string(1))),1);

                    for ins = 1:ncond
                        variable = averages_table.(ins);

                        for avs = 1:length(Moving_Average_2)
                            Moving_Average_2(avs) = Moving_Average_2(avs)+ variable(avs);
                        end
                    end

                    clear variable averages_table ave name
                end
                if Verbosity == 1 || strcmp(Mode,'review') || strcmp(Mode,'label')

                    % PLOTTING INFORMATION
                    f = figure(i);
                    set(f,'visible','off','Units','pixels','Position',[100 100 1200 800],'WindowState','normal');
                    f.Name = ['   ( ' int2str(ds) ' / ' int2str(length(infile)) ' )  Name: '  char(dsname) ' Iteration: ', num2str(it)];
                    movegui(f)

                    subplot(2,2,4),  %ERP/ERF
                    %%%% UNCORCK IT LATER!!! MD 
                    if ~strcmp(seg.studyType , 'RS')
                        plot(linspace(event(1).t_pre/1000,event(1).t_post/1000,length(Moving_Average_2/ncond)),movmean(Moving_Average_2/ncond,avg_factor));
                        xline(0,'color','#D95319','LineWidth', 2)
                        axis( [event(1).t_pre/1000 event(1).t_post/1000 1.1*min(min(Moving_Average_2/ncond)) 1.1*max(max(Moving_Average_2/ncond)) ]);
                        title('ERP/ERF','FontSize',10); xlabel('ms');
                    end
                    subplot(2,5,1:3),   %timecourse
                    plot(ica_par.time_ica,sig);
                    axis tight; title('time course','FontSize',10); xlabel('s');

                    subplot(2,2,3),   %power
                    plot(F,mspettro(it,:),'r');
                    try
                        axis([ 0  min(100,ica_par.filter.freq_high) 0 1.1*max(mspettro(it,:)) ]);  %FIXME WARNING!!! Axes are not compatible - commenting for now %MD
                    catch
                        axis([ 0 100 0 1.1*max(mspettro(it,:)) ]);  %FIXME WARNING!!! Axes are not compatible - commenting for now %MD

                    end
                    title('power spectrum','FontSize',10); xlabel('Hz');
                    if exist(saveFolderds)==0; mkdir (saveFolderds); end

                    sc_locfile( elecloc_ica, saveFolderds); 

                    if  strcmp(Mode,'label')
                        prob_values(i) = 0;
                    end
                    subplot(2,5,4:5),
                    text(-0.5,-0.6,strcat('ELC correlation:  ',{'     '}, num2str(elc_signal_corr(it),'%4.2f')),'FontSize', 10,'fontweight','bold'),
                    text(-0.5,-0.7,strcat('Brain Sign prob:  ',{'     '}, num2str(prob_values(i),'%4.2f')),'FontSize', 10,'fontweight','bold'),
                    text(-0.5,-0.8,strcat('Peaks per minute: ',{'  '}, num2str(peaks_per_min(it),'%4.2f')),'FontSize', 10,'fontweight','bold'),

                    topoplot(ICA.iter(indeks).mixing(:,it), fullfile(rc.expFolder, rc.analysisName, '6_ICA2', dsname,'ch_trimmed.elc'), 'nosedir', chan.ori, 'electrodes', 'off');
                    title(strcat({'IC'},num2str(it),' [/',num2str(nIC),']'),'FontSize',10);


                    ButtonName = '';
                    while flag2 ~= 1
                        if strcmp(Mode,'review')
                            k = waitforbuttonpress;
                            key_value = double(get(gcf,'CurrentCharacter'));
                            if double(key_value) == 29     % Pressing Right Arrow gives next component
                                ButtonName = 'Next';
                                flag2 =1;
                            elseif double(key_value) == 28 %  Pressing Left Arrow gives previous component
                                ButtonName = 'Back';
                                flag2 =1;
                            elseif double(key_value) == 13 % Pressing Enter gives next iteration
                                i = length(iterator_IC);
                                break;
                            elseif double(key_value) == 27 % Pressing Esc ends viewing
                                ButtonName = 'ex';
                                flag2 =1;
                            end
                        else
                            inps = input('Press b/a/s to label accordingly: Brain/Artifact/Notsure; type x to close','s');

                            if strcmp(inps,'x')
                                flag2 =1;
                                ButtonName ='ex';
                            elseif strcmp(inps,'b')
                                flag2 =1;

                                ButtonName ='Brain signal';
                            elseif strcmp(inps,'a')
                                ButtonName ='Artifact';
                                flag2 =1;
                            elseif strcmp(inps,'s')
                                ButtonName ='Not Sure';
                                flag2 =1;
                            elseif strcmp(inps,'<')
                                ButtonName ='Back';
                                flag2 =1;
                            elseif strcmp(inps,'') && strcmp(Mode,'review')
                                ButtonName = 'Next';
                                flag2 =1;
                            else
                                disp("Press viable button!");
                            end
                        end
                    end

                    flag2 =0;

                    disp(strcat('  User Action:      ' , ButtonName ))
                    disp(strcat('  Component Number: ' ,num2str(i)))
                    if strcmp(ButtonName,'Brain signal')
                        g_ic = [g_ic it];
                        classif_new = [classif_new 1];
                    elseif strcmp(ButtonName,'Artifact')
                        b_ic = [b_ic it];
                        classif_new = [classif_new 0];
                    elseif strcmp(ButtonName,'Not Sure')
                        not_sure_ic = [not_sure_ic it];
                    elseif strcmp(ButtonName,'Back')
                        if i > 1
                            i = i-2;
                        else
                            i = i-1;
                        end
                    elseif strcmp(ButtonName,'Next') && strcmp(Mode,'review')
                        if prob_values(i) > classification_threshold
                            g_ic = [g_ic it];
                            classif_new = [classif_new 1];
                        else
                            b_ic = [b_ic it];
                            classif_new = [classif_new 0];
                        end
                    else
                        %flag = 1;
                        if Unsure_only == true && isempty(not_sure_ic)
                            not_sure_ic = not_sure_ic2;
                        end
                        break;
                    end
                    close all;

                    if i >= nIC
                        for k= i:nIC
                            ICA_after_threshold = [ICA_after_threshold k];
                        end
                        %flag = 1;
                        break;
                    else
                        ICA_after_threshold = [];
                    end
                    i = i+1;

                else
                    if prob_values(i) > classification_threshold
                        g_ic = [g_ic it];
                        classif_new = [classif_new 1];
                    else
                        b_ic = [b_ic it];
                        classif_new = [classif_new 0];
                    end

                    if i == nIC
                        for k= i:nIC
                            ICA_after_threshold = [ICA_after_threshold k];
                        end
                        %flag = 1;
                        break;
                    end
                    i = i+1;
                end
            end

            clear sig mat;
            if Verbosity == 1
                ButtonName=questdlg('Are you sure that the classification is correct?', ...
                    'Independent Components classification','Yes','No','Yes');

                if strcmp(ButtonName,'Yes')
                    flag=1;
                end
            else
                ButtonName = 'Yes'; % WARNING
                flag=1;
            end

            % Validation
            if Validation == true && strcmp(Mode,'class')
                %Checking model accuracy
                class_rep_1 = classification_report(predictions,double(Y));
                confMat_1 = confusionmat(double(Y),predictions);
                figure('Name','Confusions')
                set(0,'DefaultFigureVisible','off');

                subplot(2,2,1)
                heatmap(confMat_1)

                subplot(2,2,2)
                %[Xs,Ys] = perfcurve(double(Y'),prob_values,1);
                %plot(Xs,Ys)

                text(.75,.275,sprintf("Accuracy: %0.2f",class_rep_1.Accuracy))
                text(.75,.2,sprintf("Precision: %0.2f",class_rep_1.Accuracy))
                text(.75,.125,sprintf("Recall: %0.2f",class_rep_1.Accuracy))
                text(.75,.05,sprintf("F1 score: %0.2f",class_rep_1.Accuracy))

                %uit = uitable("Data",struct2table(class_rep_1));

                fig_handle = gcf;
                print(gcf,'-dpng', strcat(currentFolder, '/',dsname,'/',strcat('validation_map',num2str(indeks),'.png')));
                close(fig_handle)

            end

            if ~Validation == true
                if strcmp(ButtonName,'Yes')
                    file_name = char(infile(ds));
                    disp(not_sure_ic)
                    data_ica2.ICA.iter(indeks).Unsure = not_sure_ic;
                    data_ica2.ICA.iter(indeks).brain_ic_number = length(g_ic);
                    data_ica2.ICA.iter(indeks).brain_ic = g_ic;
                    data_ica2.ICA.iter(indeks).ICA_after_threshold = ICA_after_threshold;
                    if n_elec > 3
                        art_con = mean(data_ica2.ICA.iter(indeks).elc_power_correlation(g_ic));
                    else
                        art_con = 0;
                    end
                    data_ica2.ICA.iter(indeks).artifact_contamination = art_con;
                end
            end
        end

        if ~isempty(not_sure_ic) && strcmp(Mode,'class')
            for i= not_sure_ic
                if prob_values(i) > classification_threshold
                    g_ic = [g_ic i ];
                else
                    b_ic = [b_ic i ];
                end
            end
            not_sure_ic = [];
        end

        close all
    end
    
    %% Calculating Artifact Contamination
    for it = 1: num_iters
        corresp_g_ic = data_ica2.ICA.iter(it).brain_ic;

        if n_elec > 2
            art_con = mean(data_ica2.ICA.iter(it).elc_power_correlation(corresp_g_ic));
        else
            art_con = 0;
        end

        goodness_ic(it,1) = length(corresp_g_ic);   % no of brain components
        goodness_ic(it,2) = 1-art_con;           % artifact contamination

    end
    [~,bestiter_new] = max( goodness_ic(:,1).*goodness_ic(:,2) );

    data_ica2.ICA.best_iter.index = bestiter_new;
    data_ica2.ICA.best_iter.IC_sig =  data_ica2.ICA.iter(bestiter_new).IC_sig;
    data_ica2.ICA.best_iter.nonlinearity = ica.ICA2_nonlinearity;

       

    if (strcmp(rc.dataFormat, 'itab-sim') || strcmp(rc.dataFormat, 'syn-sim')) && strcmp(Mode,'class')
        if isfield(rc.simul, 'setallICasbrain')
            if rc.simul.setallICasbrain   %settinga all brain as IC - for sim purposes only!!
                for it =  1:length(data_ica2.ICA.iter)
                    nofic(it)= data_ica2.ICA.iter(it).total_ic_number;
                end
                [icnobestiter,bestindex] = max(nofic);
                data_ica2.ICA.best_iter.index = bestindex;
                data_ica2.ICA.iter(bestindex).brain_ic = 1:icnobestiter;
                data_ica2.ICA.iter(bestindex).brain_ic_number = data_ica2.ICA.iter(bestindex).total_ic_number;
            end
        end
    end

    %prepare IC statistics

    datastat = sc_icavarstat(dsname, infile{ds}, data_ica2);

    ica2stat{1} = dsname;
    ica2stat{2} = bestiter_new;
    ica2stat{3} = data_ica2.ICA.iter(bestiter_new).total_ic_number;
    ica2stat{4} = data_ica2.ICA.iter(bestiter_new).brain_ic_number;
    ica2stat(5) = {datastat.var_ave(1)};
    ica2stat(6) = {datastat.var_ave(2)};
    ica2stat(7) = {datastat.var_ave(3)};
    ica2stat(8) = {datastat.var_ave(2) - datastat.var_ave(3)};
    ica2stat(9) = {round(datastat.ratio_brain)};
    ica2stat(10) = {round(datastat.ratio_resid)}; 

    %%%%% SAVING

    file_name = infile{ds};
    fprintf('final_saving at %s \n',file_name)
    save(file_name,'data_ica2','-v7.3');  % '-v6' v6 added for speedup

    clear goodness_ic
    fid = fopen(progressfile, 'w');
    fprintf(fid, 'recently completed ds: ');
    fprintf(fid, dsname);
    fprintf(fid, '\n');
    fprintf(fid, mat2str(ds));
    fclose(fid);

    %%%% Making best_iter TOPOPLOT %%%%%

    disp('')
    currentFolder = fullfile(rc.analysisFolder, '6_ICA2' );
    it = data_ica2.ICA.best_iter.index;
    % iter = data_ica2.ICA.iter(it);
    g_ic = data_ica2.ICA.iter(it).brain_ic;
    fprintf('Best iter: %d',it)
    disp('')
    all_ic = [1:size(data_ica2.ICA.iter(it).IC_sig,1)];
    b_ic = all_ic(~ismember(all_ic,g_ic));

    %%% Ploting Brain ICs

    chan2ica = [find(strcmp(data_ica2.hdr.actualchantype, 'eeg'))' find(strcmp(data_ica2.hdr.actualchantype, 'meg'))' find(strcmp(data_ica2.hdr.actualchantype, 'mag'))' find(strcmp(data_ica2.hdr.actualchantype, 'meggrad'))' find(strcmp(data_ica2.hdr.actualchantype, 'megmag'))'];
    if rc.eeg
        elecloc_ica = sc_trim_sens(data_ica2.label(chan2ica), data_ica2.hdr.elec);
    else
        elecloc_ica = sc_trim_sens(data_ica2.label(chan2ica), data_ica2.hdr.grad);  %provides labels of chans that survived AR and were used for ICA
    end
    sc_locfile( elecloc_ica, saveFolderds); %trims elc file according to actual no of channels
    topodir = '+X';

   % figure ('visible','off','Position',get(0,'Screensize'));
   figure ('visible','off',   'position',[0,0,4000,3000]);
                
    %%% Adding Failsafe if no ICA Is Good 
    if ~isempty(g_ic)
        lat_good =  ceil( sqrt(length(g_ic)) );
    else
        lat_good = 1;
        fileEXIT = fopen(fullfile(currentFolder,dsname,"NO_GOOD_ICA_WHATSOEVER.txt"),"w");
        fprintf(fileEXIT,"In All Iterations No single Good ICA Had Been Found! Remove Channel IMMEDIATELY")
    end


    tiledlayout(lat_good,lat_good,'TileSpacing','tight','Padding','tight');   %original subplot layout
    for ixrev=1:length(g_ic) % NO PARFOR for now because of No generation
        ix=length(g_ic)-ixrev+1; %par draws ICs in reversed order
        temp_ix = g_ic(ix);
        nexttile
        if rc.eeg
        topoplot(data_ica2.ICA.iter(it).mixing(:,temp_ix), char(fullfile(trim_folder,'ch_trimmed.elc')), 'nosedir', chan.ori, 'electrodes', 'labels');
        else
        topoplot(data_ica2.ICA.iter(it).mixing(:,temp_ix), char(fullfile(trim_folder,'ch_trimmed.elc')), 'nosedir', chan.ori, 'electrodes', 'on');
        end
        title(strcat({'IC'},num2str(temp_ix),' /',num2str(length(g_ic)),''));
    end
    try
        delete(strtrim(ls(fullfile(currentFolder, dsname, '_BRAIN_BEST_ITER*.jpg'))));    % remove old files
    catch
    end
    fig_handle = gcf();
    print(fig_handle,  strcat(currentFolder, '/',dsname,'/', strcat('_BRAIN_BEST_ITER_',num2str(it), '_MAPS.jpg')),'-djpeg','-r180');
    close(fig_handle)
      %%% Ploting nonBrain ICs
    
   figure ('visible','off',   'position',[0,0,4000,3000]);

    if ~isempty(b_ic)
        lat_bad =  ceil( sqrt(length(b_ic)) );
    else
        lat_bad = 1;
    end
        
    tiledlayout(lat_bad,lat_bad,'TileSpacing','tight','Padding','tight');   %original subplot layout
    for ixrev=1:length(b_ic)
        ix=length(b_ic)-ixrev+1; %par draws ICs in reversed order
        temp_ix = b_ic(ix);
        nexttile
        if rc.eeg
            topoplot(data_ica2.ICA.iter(it).mixing(:,temp_ix), fullfile(trim_folder, 'ch_trimmed.elc'), 'nosedir', chan.ori, 'electrodes', 'labels');
        else
            topoplot(data_ica2.ICA.iter(it).mixing(:,temp_ix), fullfile(trim_folder, 'ch_trimmed.elc'), 'nosedir', chan.ori, 'electrodes', 'on');
        end
        title(strcat({'IC'},num2str(temp_ix),' /',num2str(length(g_ic)),''));
    end

    try
        delete(strtrim(ls(fullfile(currentFolder, dsname, '_NONBRAIN_BEST_ITER*.jpg'))));  % remove old files
    catch
    end

    fig_handle = gcf();
    print(fig_handle, strcat(currentFolder, '/',dsname,'/',strcat('_NONBRAIN_BEST_ITER_',num2str(it), '_MAPS.jpg')),'-djpeg','-r180');
    close(fig_handle)

    % clearvars -except Cx_map DC Data_volume IF_conv_model IF_mlp_model Mode Unsure_only Validation Verbosity ar begId chan classif_table_all classification_threshold con  currentFolder  design  event_ctftrig  ica  ica2stat  infile  loc  model_type  predictions_all  prob_values_all  rc  seg  sim ica2stat 
 
end
