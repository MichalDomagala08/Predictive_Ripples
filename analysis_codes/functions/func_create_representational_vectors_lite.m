function R_unified = func_create_representational_vectors_lite(R_unified,ripple_data,ripple_table,dataFolder,params)
    % This function uses similar preprocessing but to gather exact times
    % from ripples that are provided

    R_unified(R_unified.is_after_sac,:).event_rip_latency  = -R_unified(R_unified.is_after_sac,:).rip_event_distance_edge;
    R_unified(R_unified.is_before_sac,:).event_rip_latency = R_unified(R_unified.is_before_sac,:).rip_event_distance_edge;
    R_unified(R_unified.is_after_fix,:).event_rip_latency  = -R_unified(R_unified.is_after_fix,:).rip_event_distance_edge;
    R_unified(R_unified.is_before_fix,:).event_rip_latency = R_unified(R_unified.is_before_fix,:).rip_event_distance_edge;


    R_unified  = R_unified(R_unified.is_abs_closest,:);
    R_unified.R_reverse = nan(height(R_unified),1);


    files = {dir(fullfile(dataFolder, '*.mat')).name};

    for sb = 1:length(ripple_data.all_subjNames )

        raw_data = load(fullfile(dataFolder,files{sb}),'data').data.data_eeg;

        fprintf("\n SUBJECT: %s\n\n",ripple_data.all_subjNames{sb})

        if params.data_preprocessing

            if exist(fullfile("D:\Documents_Dell\Predictive_Ripples_2025\data\preprocessed",strrep(files{sb},"reref","preprocessed")))
               data = load(fullfile("D:\Documents_Dell\Predictive_Ripples_2025\data\preprocessed",strrep(files{sb},"reref","preprocessed")));
                data_artifact = data.data_artifact;
                if params.bha; data_bha = data.data_bha;end

            else

                [~,~,~,~,data_artifact,~]   = func_artifact_rejection_alternative(raw_data,params);
                if params.bha; [data_bha,~] = func_rippleband_filtering(raw_data,params);end
                
                save(fullfile("D:\Documents_Dell\Predictive_Ripples_2025\data\preprocessed",strrep(files{sb},"reref","preprocessed")),"data_artifact","data_bha");

            end
        end


        %%
        currentRipples_sac_R_T = [];


        %%
        for it = 1:length(data_artifact.trial)
            
    
            if params.bha % If we are addng BHA data to the mix
                data_bha.trial{it}(isnan(data_artifact.trial{it})) = nan;
                currTrialData = [data_artifact.trial{it} - mean(data_artifact.trial{it},2,"omitnan");...
                                 data_bha.trial{it} -  mean(data_bha.trial{it},2,"omitnan")];
            else
                currTrialData = data_artifact.trial{it} - mean(data_artifact.trial{it},2,"omitnan");
            end


            %%% (1) Ripple Locked Representation (Average) 
            currentRipples = table2array(ripple_table(ripple_table.subject ==ripple_data.all_subjNames{sb} & ripple_table.trial == it,[5,6] ));
            currentRipples_sac_R = nan(size(currentRipples,1),4);

            currentRipples_sac_R(:,3:4) = currentRipples;

         
            for rip = 1:size(currentRipples,1)
                % Get ripple representation:
                curr_rip_repr = mean(currTrialData(:,currentRipples(rip,1):currentRipples(rip,2)),2,"omitnan");
            
                curr_cl_events = R_unified(R_unified.subject ==ripple_data.all_subjNames{sb} & R_unified.trial == it & R_unified.rip_beg == currentRipples_sac_R(rip,3) ,:);

        
                try

%%
                    % Reverse Saccade  and Ripple repr:
                    for i = [0,1]
                        curr_sac_ev = curr_cl_events(curr_cl_events.is_fix == i,:);
    
                        if height(curr_sac_ev) >1; curr_sac_ev = curr_sac_ev(1,:); end
                        % If Saccade is AFTER ripple (then we take exactly same distance before)
                        if curr_sac_ev.rank_cat == 'a'
                            sac_time = curr_sac_ev.event_rip_latency + curr_sac_ev.rip_beg;
                            sac_len  = curr_sac_ev.event_end - curr_sac_ev.event_beg;
        
                            if  sac_time-sac_len <= 0 ; continue; else
                                curr_sac_repr = mean(currTrialData(:,sac_time-sac_len:sac_time),2,"omitnan");
                                ok = ~isnan(curr_sac_repr) & ~isnan(curr_rip_repr);
                
                                if sum(ok) >= 5  % minimalna liczba punktów, ustaw wg. siebie
                                     currentRipples_sac_R(rip,i+1) = corr(curr_sac_repr(ok), curr_rip_repr(ok), 'type', 'Pearson');
                                end
                            end
                        else
                            sac_time = curr_sac_ev.event_rip_latency + curr_sac_ev.rip_end;
                            sac_len  = curr_sac_ev.event_end - curr_sac_ev.event_beg;
        
                             if  sac_time+sac_len > 7001 ; continue; else
                                curr_sac_repr = mean(currTrialData(:,sac_time:sac_time+sac_len),2,"omitnan");
                                ok = ~isnan(curr_sac_repr) & ~isnan(curr_rip_repr);
                
                                if sum(ok) >= 5  % minimalna liczba punktów, ustaw wg. siebie
                                     currentRipples_sac_R(rip,i+1) = corr(curr_sac_repr(ok), curr_rip_repr(ok), 'type', 'Pearson');
                                end
                            end
                        end
    
                        % Saving currently done ripple to an Orignal Table
                        slice = R_unified.subject ==ripple_data.all_subjNames{sb} & R_unified.trial == it & R_unified.is_fix==i & R_unified.rip_beg == currentRipples_sac_R(rip,3);
                        
                            if height(R_unified(slice,:)) > 1
                                R_unified(find(slice,1),:).R_reverse = currentRipples_sac_R(rip,i+1);
                            elseif height(R_unified(slice,:))
                                R_unified(slice,:).R_reverse= currentRipples_sac_R(rip,i+1);
                            end
                      
                    end
%%
                catch
                    sprintf("%s Trial: %d Ripple: %d",ripple_data.all_subjNames{sb},it,rip)
                    disp("error")
                end
            end
          

    
        end

%%
    end
    
end


