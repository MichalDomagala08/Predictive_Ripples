function [subject_level_repr_rp, subject_level_repr_sac, subject_level_repr_fix] = func_create_representational_vectors(ripple_table,saccades_table,ripple_data,dataFolder,saveFolder,params)
    subject_level_repr_rp = cell(1,19);
    subject_level_repr_sac = cell(1,19);
    subject_level_repr_fix = cell(1,19);
    files = {dir(fullfile(dataFolder, '*.mat')).name};
    % Creates representational vectors for ripples, saccades and fixations
    for sb = 1:length(ripple_data.all_subjNames )
    
        
        fprintf("\n SUBJECT: %s\n\n",ripple_data.all_subjNames{sb})
    
        %%% load and unpack the data structu res
        raw_data = load(fullfile(dataFolder,files{sb}),'data').data.data_eeg;
    
        subject_level_repr_rp{sb} = cell(1,params.trialNum);
        subject_level_repr_sac{sb} = cell(1,params.trialNum);
        subject_level_repr_fix{sb} = cell(1,params.trialNum);
        %%% Get Artifact Estimation across all channels - 

        if params.data_preprocessing == 2 % Kasia Processing
            marks =get_artifact_marks(ripple_data.all_subjNames{sb}, raw_data, saccades_table(saccades_table.subject_new == ripple_data.all_subjNames{sb},:));
            [data_artifact, ~] = apply_artifact_marks(raw_data, marks);


        elseif params.data_preprocessing == 1 % My, OLD Processing
            [~,~,~,~,data_artifact,~] =  func_artifact_rejection_alternative(raw_data,params);% OLD PROCESSING
        else                                  % Read Processed data 
            curpath = strsplit(dataFolder,"\"); if curpath(5) == "preprocessing"; else error("Path should be preprocessing"); end
            data_artifact = raw_data.data_artifact;
            data_bha      = raw_data.data_bha;
        end

        if params.bha; [data_bha,~] = func_rippleband_filtering(raw_data,params);end
    
        for it = 1:length(raw_data.trial)
            
    
            if params.bha % If we are addng BHA data to the mix
                data_bha.trial{it}(isnan(data_artifact.trial{it})) = nan;
                currTrialData = [data_artifact.trial{it} - mean(data_artifact.trial{it},2,"omitnan");...
                                 data_bha.trial{it} -  mean(data_bha.trial{it},2,"omitnan")];
            else
                currTrialData = data_artifact.trial{it} - mean(data_artifact.trial{it},2,"omitnan");
            end


            %%% (1) Ripple Locked Representation (Average) 
            currentRipples = table2array(ripple_table(ripple_table.subject ==ripple_data.all_subjNames{sb} & ripple_table.trial == it,[5,6] ));
            subject_level_repr_rp{sb}{it} = zeros(size(currTrialData,1),height(currentRipples));
            for rip = 1:size(currentRipples,1)
                subject_level_repr_rp{sb}{it}(:,rip) = mean(currTrialData(:,currentRipples(rip,1):currentRipples(rip,2)),2,"omitnan");
            end
    
    
            %%% (2) Saccade Loced Representation (Average) 
            curr_sac = saccades_table(strcmp(saccades_table.subject_new,ripple_data.all_subjNames{sb}) & saccades_table.trial_number == it,{'saccade_onset_time','saccade_offset_time'});
            curr_sac  = curr_sac(~isnan(curr_sac.saccade_onset_time) & ~isnan(curr_sac.saccade_offset_time),:);
            subject_level_repr_sac{sb}{it} = zeros(size(currTrialData,1),length(height(curr_sac)));
    
            for sac = 1:height(curr_sac)
                currsactime = round((curr_sac(sac,:).saccade_onset_time*500+2000):(curr_sac(sac,:).saccade_offset_time*500+2000));
                subject_level_repr_sac{sb}{it}(:,sac) = mean(currTrialData(:,currsactime),2,"omitnan");
            end


            %%% (3) Fixations Locked Representation (Overlapping windows of 25 ms (or N) 
            curr_fix = saccades_table(strcmp(saccades_table.subject_new,ripple_data.all_subjNames{sb}) & saccades_table.trial_number == it,{'fixation_beg','fixation_end','fixation_len'});
            curr_fix  = curr_fix(~isnan(curr_fix.fixation_beg) & ~isnan(curr_fix.fixation_end),:);
            curr_fix  = curr_fix(1:end-1,:);

            subject_level_repr_fix{sb}{it} = zeros(size(currTrialData,1),length(height(curr_fix)));

            nWin = ceil(max(curr_fix.fixation_len)*500/params.fix_box_length);

            subject_level_repr_fix{sb}{it} = nan(size(currTrialData,1), height(curr_fix),max(nWin));
            for fix = 1:height(curr_fix)
                fixationSamples =  round((curr_fix(fix,:).fixation_beg*500+2000):(curr_fix(fix,:).fixation_end*500+2000));
                disp(ceil(length(fixationSamples)/params.fix_box_length))

                for win = 1:ceil(length(fixationSamples)/params.fix_box_length)
                    idx = (win-1)*params.fix_box_length+1 : min(win*params.fix_box_length, numel(fixationSamples));
                    subject_level_repr_fix{sb}{it}(:,fix,win) = mean(currTrialData(:,fixationSamples(idx)),2,"omitnan");
                end
            end
        end

        if params.data_preprocessing == 1 && params.data_save
            data.data_artifact = data_artifact;
            data.data_bha = data_bha;
            save(fullfile('D:\Documents_Dell\Predictive_Ripples_2025\data\preprocessed',files{sb}),"data",'-v7.3');
         elseif params.data_preprocessing == 2 && params.data_save
             data.data_artifact = data_artifact;
            data.data_bha = data_bha;
            save(fullfile('D:\Documents_Dell\Predictive_Ripples_2025\data\preprocessed_kasia',files{sb}),"data",'-v7.3');
        end
    
    end
    if params.bha; reprTitle = "representational_data_bha.mat"; else; reprTitle = "representational_data.mat"; end

end