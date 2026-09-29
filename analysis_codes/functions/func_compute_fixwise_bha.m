function [fixation_bha_table_ret,fixation_bha_long] = func_compute_fixwise_bha(fixation_bha_table,ripple_table,ripple_data,dataFolder)

    % This function computes bha (by either loading data or filtering in-place)
    % and then computes:
    % - mean BHA per fixation
    % - Before vs After BHA
    % - replicates with and without ripples 
    % - replicates with channels


    files = {dir(fullfile(dataFolder, '*.mat')).name};


     % Electrode info
    T = readtable('D:\Documents_Dell\Predictive_Ripples_2025\data\Electrodes.xlsx','ReadRowNames',true, 'Sheet', 'responsive');
    electrodeFields = {'antHP',	'postHP'};

    % Select Rippleband channels only:
    clear goodChannels
    for sb = 1:length(T.Properties.RowNames)
        channel_strings = cellfun(@(f) T.(f){sb}, electrodeFields, 'UniformOutput', false);
        channel_strings = channel_strings(~cellfun('isempty', channel_strings));  % usuń puste
        splits          = cellfun(@(s) split(s,','), channel_strings, 'UniformOutput', false);
        goodChannels{sb} = vertcat(splits{:});
    end
    goodChannels(cellfun(@isempty, goodChannels)) = [];

    gc = cellfun(@(c) string(c(:)), goodChannels, 'UniformOutput', false);   % każdy kanał jako string
    uniqueChans = unique(vertcat(gc{:}), 'stable');
    for i = 1:length(uniqueChans)
        fixation_bha_table.(uniqueChans(i)+"_mean") = nan(height(fixation_bha_table),1);
        fixation_bha_table.(uniqueChans(i)+"_mean_norip") = nan(height(fixation_bha_table),1);
        fixation_bha_table.(uniqueChans(i)+"_diff") = nan(height(fixation_bha_table),1);
        fixation_bha_table.(uniqueChans(i)+"_diff_norip") = nan(height(fixation_bha_table),1);
    end



    for sb = 1:length(ripple_data.all_subjNames)

        raw_data = load(fullfile(dataFolder,files{sb}));

        % Perform Artifact Rejection for Data
        if params.data_preprocessing == 2 % Kasia Processing
            marks =get_artifact_marks(ripple_data.all_subjNames{sb}, raw_data, fixation_bha_table(fixation_bha_table.subject_new == ripple_data.all_subjNames{sb},:));
            [data_artifact, ~] = apply_artifact_marks(raw_data, marks);
            [data_bha,~] = func_rippleband_filtering(raw_data,params);
            data_bha.trial{it}(isnan(data_artifact.trial{it})) = nan;
        elseif params.data_preprocessing == 1 % My, OLD Processing
            [~,~,~,~,data_artifact,~] =  func_artifact_rejection_alternative(raw_data,params);% OLD PROCESSING
            [data_bha,~] = func_rippleband_filtering(raw_data,params);
            data_bha.trial{it}(isnan(data_artifact.trial{it})) = nan;
        else                                  % Read Processed data 
            curpath = strsplit(dataFolder,"\"); if curpath(5) == "preprocessing"; else error("Path should be preprocessing"); end
            data_bha      = raw_data.data_bha;
        end
        chanHipIdx = ismember(data_bha.label,goodChannels{sb});
    
        % Wyniki trzymam w macierzy zamiast wpisywać do tabeli w najgłębszej pętli.
        bha_out = nan(height(fixation_bha_table), numel(goodChannels{sb}), 4);
    
        % strcmp po całej tabeli raz na osobę, a nie dla każdej fiksacji x kanału.
        subjMask = strcmp(fixation_bha_table.subject_new,ripple_data.all_subjNames{sb});
    
        for it = 1:length(data_bha.trial)
            fprintf("Subject: %s Trial: %d \n",ripple_data.all_subjNames{sb},it)
            data_bha.trial{it} = data_bha.trial{it}.^2; % Squared for power
            data_bha.trial{it}(isnan(raw_data.data_artifact.trial{it})) = nan;
    
            % Get ripple and fixation data
            currentRipples = ripple_table(ripple_table.subject ==ripple_data.all_subjNames{sb} & ripple_table.trial == it,[2,5,6] );
            curr_fix = fixation_bha_table(subjMask & fixation_bha_table.trial_number == it,8:12); % [ZMIANA 2] subjMask zamiast strcmp
            %  Zapamiętuję numery wierszy curr_fix w fixation_bha_table (użyte niżej zamiast indx_sac == fix-1)
            curr_fix_rows = find(subjMask & fixation_bha_table.trial_number == it);
    
            % Fallback when there are empty fixation rows
            if sum(curr_fix.fixation_beg) == 0 & ~isempty(curr_fix)
                curr_fix.fixation_beg = curr_fix.saccade_offset_time;
                curr_fix.fixation_end = ([curr_fix.saccade_onset_time(2:end); 6]);
                curr_fix.fixation_len = curr_fix.fixation_end - curr_fix.fixation_beg;
                %     fixation_bha_table(strcmp(fixation_bha_table.subject_new,ripple_data.all_subjNames{sb}) &...
                %                               fixation_bha_table.trial_number == it,:).fixation_beg = curr_fix.fixation_beg;
                %     fixation_bha_table(strcmp(fixation_bha_table.subject_new,ripple_data.all_subjNames{sb}) &...
                %                               fixation_bha_table.trial_number == it,:).fixation_end = curr_fix.fixation_end;
                %     fixation_bha_table(strcmp(fixation_bha_table.subject_new,ripple_data.all_subjNames{sb}) &...
                %                               fixation_bha_table.trial_number == it,:).fixation_len = curr_fix.fixation_len;
            end
            curr_fix.next_saccade_onset  = curr_fix.fixation_end;
            curr_fix.next_saccade_offset = [curr_fix.saccade_offset_time(2:end); 6 + nanmean(curr_fix.saccade_offset_time - curr_fix.saccade_onset_time)];
    
            %  Ta sama maska isnan
            keepFix = ~isnan(curr_fix.fixation_beg) & ~isnan(curr_fix.fixation_end);
            curr_fix  = curr_fix(keepFix,:);
            curr_fix_rows = curr_fix_rows(keepFix);
            % curr_fix  = curr_fix(1:end-1,:);
            % Iterate through fixation events
            for fix = 1:height(curr_fix)
                % Get timings of fixations and saccades before and after
                fixationSamples  =  round((curr_fix(fix,:).fixation_beg*500+2000):(curr_fix(fix,:).fixation_end*500+2000));
                sacBeforeSamples =  round((curr_fix(fix,:).saccade_onset_time*500+2000):(curr_fix(fix,:).saccade_offset_time*500+2000));
                sacAfterSamples  =  round((curr_fix(fix,:).next_saccade_onset*500+2000):(curr_fix(fix,:).next_saccade_offset*500+2000));
    
                % [ZMIANA 3] Wiersz docelowy = wiersz, z którego wzięto czasy fiksacji.
                % Uzasadnienie: indx_sac == fix-1 jest poprawne tylko wtedy, gdy żaden wiersz nie wypadł
                % przy filtrze ~isnan powyżej. Jeśli wypadł wiersz indx_sac = 0, wszystko przesuwa się o 1.
                currentFixMain = curr_fix_rows(fix);
                % currentFixMain = strcmp(fixation_bha_table.subject_new,ripple_data.all_subjNames{sb}) &...
                %                  fixation_bha_table.trial_number == it & fixation_bha_table.indx_sac == fix-1 ;
                
                % (2) Per_Channel_BHA:
                channelwise_bha = {}; channelwise_bha_norip = {}; channelwise_bha_diff = {};  channelwise_bha_norip_diff = {};
                for i = 1:numel(goodChannels{sb})
                    currentRipples_chan = currentRipples(currentRipples.channel == goodChannels{sb}{i},[2,3]);
                    ranges = arrayfun(@(b, e) b:e, currentRipples_chan.beg, currentRipples_chan.end, 'UniformOutput', false);
                    rippleSamples = [ranges{:}];
                    channelslice = ismember(data_bha.label,goodChannels{sb}{i});
                    if ~ isempty(data_bha.trial{it}(channelslice))
                        channelwise_bha{i} = data_bha.trial{it}(channelslice,fixationSamples);
                        channelwise_bha_norip{i} = data_bha.trial{it}(channelslice, setdiff(fixationSamples,rippleSamples));
                        channelwise_bha_diff{i} = nanmean(data_bha.trial{it}(channelslice,sacAfterSamples)) - nanmean(data_bha.trial{it}(channelslice,sacBeforeSamples));
                        channelwise_bha_norip_diff{i} = nanmean(data_bha.trial{it}(channelslice, setdiff(sacAfterSamples,rippleSamples))) -...
                                                        nanmean(data_bha.trial{it}(channelslice, setdiff(sacBeforeSamples,rippleSamples)));
                        % [ZMIANA 1] zapis do macierzy zamiast do tabeli nanmean zamiast mean. Uzasadnienie: jedna próbka artefaktu (NaN)
                        % dawała NaN dla całej fiksacji; w _diff i tak używałeś nanmean.
                        bha_out(currentFixMain, i, 1) = nanmean(channelwise_bha{i});
                        bha_out(currentFixMain, i, 2) = nanmean(channelwise_bha_norip{i});
                        bha_out(currentFixMain, i, 3) = channelwise_bha_diff{i};
                        bha_out(currentFixMain, i, 4) = channelwise_bha_norip_diff{i};
                        % fixation_bha_table( currentFixMain,:).(goodChannels{sb}{i}+"_mean")       = mean(channelwise_bha{i});
                        % fixation_bha_table( currentFixMain,:).(goodChannels{sb}{i}+"_mean_norip") = mean(channelwise_bha_norip{i});
                        % fixation_bha_table( currentFixMain,:).(goodChannels{sb}{i}+"_diff")       = channelwise_bha_diff{i};
                        % fixation_bha_table( currentFixMain,:).(goodChannels{sb}{i}+"_diff_norip") = channelwise_bha_norip_diff{i};
                    end
                end
            end
        end
    
        % [ZMIANA 1] Zapis do tabeli raz na osobę, tylko wiersze tej osoby.
        % Jeśli kolumny nie ma, tworzę ją jako NaN, żeby wiersze bez danych nie dostały domyślnej wartości.
        suffixes = ["_mean", "_mean_norip", "_diff", "_diff_norip"];
        for i = 1:numel(goodChannels{sb})
            for m = 1:4
                colName = goodChannels{sb}{i} + suffixes(m);
                if ~ismember(colName, fixation_bha_table.Properties.VariableNames)
                    fixation_bha_table.(colName) = nan(height(fixation_bha_table),1);
                end
                fixation_bha_table.(colName)(subjMask) = bha_out(subjMask, i, m);
            end
        end
    end

    % Compute means
    % Compute means
    chanCols = uniqueChans + suffixes;                          % [nKanałów x 4], np. "RDa1-RDa2_mean"

    for m = 1:numel(suffixes)
        fixation_bha_table.("chanAvg" + suffixes(m)) = mean(fixation_bha_table{:, cellstr(uniqueChans + suffixes(m))}, 2, 'omitnan');
    end
    fixation_bha_table.chanAvg_nChan = sum(~isnan(fixation_bha_table{:, cellstr(chanCols(:, 1))}), 2);
    fixation_bha_table = renamevars(fixation_bha_table,["subject_new","trial_number"],["subject","trial"]);
    fixation_bha_table_ret = fixation_bha_table(:,[1:12,177:181]);


    % (2) Tall: wiersz = fiksacja x kanał
    baseCols = setdiff(string(fixation_bha_table.Properties.VariableNames), ...
                       [chanCols(:); "chanAvg" + suffixes(:); "chanAvg_nChan"], 'stable');
    blocks = cell(numel(uniqueChans), 1);
    for c = 1:numel(uniqueChans)
        V    = fixation_bha_table{:, cellstr(chanCols(c, :))};
        keep = any(~isnan(V), 2);
        if ~any(keep), continue; end
        Tl = fixation_bha_table(keep, cellstr(baseCols));
        Tl.channel    = repmat(uniqueChans(c), sum(keep), 1);
        Tl.channel_id = string(Tl.subject) + "_" + Tl.channel;
        for m = 1:numel(suffixes)
            Tl.("BHA" + suffixes(m)) = V(keep, m);
        end
        blocks{c} = Tl;
    end
    fixation_bha_long = vertcat(blocks{:});
    

end