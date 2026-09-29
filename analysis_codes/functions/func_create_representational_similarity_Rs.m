function [T_all,T_rip,T_sac,T_fix,T_fix2] = func_create_representational_similarity_Rs(subject_level_repr_rp,subject_level_repr_sac,subject_level_repr_fix,saccades_table,ripple_table,ripple_data,saveFolder,params)
    % This function computes for every trial representationa lsimilariy of
    % every:
    % - saccade x ripple
    % - saccade x saccade
    % - ripple  x ripple
    % - fixation x ripple
    % - fixation x fixation

    T_all = nan; T_rip = nan; T_sac = nan; T_fix = nan; T_fix2 = nan;

    
    for sb = 1:length(ripple_data.all_subjNames )
        for it = 1:params.trialNum
    
            current_sac_rep = subject_level_repr_sac{sb}{it};
            current_rip_rep = subject_level_repr_rp{sb}{it};
            current_fix_rep = subject_level_repr_fix{sb}{it};
    
            curr_sac = saccades_table(strcmp(saccades_table.subject_new,ripple_data.all_subjNames{sb}) & saccades_table.trial_number == it,:);
            curr_sac  = table2array(round(curr_sac(~isnan(curr_sac.saccade_onset_time) & ~isnan(curr_sac.saccade_offset_time),{'saccade_onset_time','saccade_offset_time'}).*500 + 2000));
            curr_fix = saccades_table(strcmp(saccades_table.subject_new,ripple_data.all_subjNames{sb}) & saccades_table.trial_number == it,:);
            curr_fix  = table2array(round(curr_fix(~isnan(curr_fix.fixation_beg) & ~isnan(curr_fix.fixation_end),{'fixation_beg','fixation_end'}).*500 + 2000));
            curr_fix = curr_fix(1:end-1,:);
            curr_rip =  table2array(ripple_table(strcmp(ripple_table.subject,ripple_data.all_subjNames{sb}) & ripple_table.trial == it,1:7));
    
            %%% (1) Saccade vs Ripple Representational Similarity 
            if ~isempty(curr_sac) & ~isempty(curr_rip)
    
                R = nan(size(current_sac_rep,2), size(current_rip_rep,2));
                for i = 1:size(current_sac_rep,2)
                    for j = 1:size(current_rip_rep,2)
                        ok = ~isnan(current_sac_rep(:, i)) & ~isnan(current_rip_rep(:, j));
        
                        if sum(ok) >= 5  % minimalna liczba punktów, ustaw wg. siebie
                            R(i,j) = corr(current_sac_rep(ok, i), current_rip_rep(ok, j), 'type', 'Pearson');
                        else; R(i,j) = NaN; end
                    end
                end
        
                [I,J] = ndgrid(1:size(current_sac_rep,2), 1:size(current_rip_rep,2));
                SacIdx = I(:);   RipIdx = J(:);
                SacOn = curr_sac(SacIdx,1); SacOff = curr_sac(SacIdx,2);
                Subject = curr_rip(RipIdx,1);  Trial = curr_rip(RipIdx,3);  Channel = curr_rip(RipIdx,2); 
                
                RipBeg = curr_rip(RipIdx,5); RipMid = curr_rip(RipIdx,7); RipEnd = curr_rip(RipIdx,6);
                Rvec = R(:);
                T = table(Subject, Trial, Channel,SacIdx,RipIdx, SacOn, SacOff, RipBeg, RipEnd,RipMid, Rvec, ...
                    'VariableNames',{'subject','trial','channel','sac_idx','rip_idx','sac_onset','sac_offset','rip_beg','rip_end','rip_middle','R'});
                T = T(~isnan(T.R), :); T.rip_beg = str2double( T.rip_beg); T.rip_end = str2double( T.rip_end); T.rip_middle = str2double( T.rip_middle); 
                T.R_fisher = atanh(min(max(T.R, -0.9999), 0.9999));  T.trial = str2double( T.trial);
                T.rip_sac_timing = T.sac_onset - T.rip_beg;
        
                if ~istable(T_all); T_all = T([],:); end
                T_all = [T_all; T];
            end
            
            %%% (2) Ripple vs Ripple Representational Similarity 
            if ~isempty(curr_rip)
    
                R_rip  = nan(size(current_rip_rep,2), size(current_rip_rep,2));
                for i = 1:size(current_rip_rep,2)-1
                    for j = i+1:size(current_rip_rep,2)
                        ok = ~isnan(current_rip_rep(:, i)) & ~isnan(current_rip_rep(:, j));
                        if sum(ok) >= 5  % minimalna liczba punktów, ustaw wg. siebie (Nie posiadających NANów
                            R_rip(i,j) = corr(current_rip_rep(ok, i), current_rip_rep(ok, j), 'type', 'Pearson');
                        else; R_rip(i,j) = NaN; end
        
                    end
                end
        
        
                [I,J] = ndgrid(1:size(current_rip_rep,2), 1:size(current_rip_rep,2));
                Rip2Idx = I(:); RipIdx = J(:);  Rvec = R_rip(:);
                Subject = curr_rip(RipIdx,1);  Trial = curr_rip(RipIdx,3);  Channel = curr_rip(RipIdx,2); 
                Rip1Beg =str2double(curr_rip(RipIdx,5)); Rip1Mid = str2double(curr_rip(RipIdx,7)); Rip1End = str2double(curr_rip(RipIdx,6));
                Rip2Beg = str2double(curr_rip(Rip2Idx,5)); Rip2Mid = str2double(curr_rip(Rip2Idx,7)); Rip2End =str2double( curr_rip(Rip2Idx,6));
                
                T = table(Subject, Trial, Channel,RipIdx,Rip2Idx, Rip1Beg, Rip1End, Rip1Mid, Rip2Beg, Rip2End, Rip2Mid, Rvec, ...
                    'VariableNames',{'subject','trial','channel','rip_idx_1','rip_idx_2','rip_beg_1','rip_end_1','rip_middle_1','rip_beg_2','rip_end_2','rip_middle_2','R'});
                        T = T(~isnan(T.R),:);
        
                T.R_fisher = atanh(min(max(T.R, -0.9999), 0.9999));  T.trial = str2double( T.trial);
                T.rip_sac_timing = T.rip_beg_1 - T.rip_beg_2;
        
                if ~istable(T_rip); T_rip = T([],:); end
                T_rip = [T_rip; T];
            end
    
            %%% (3) Saccade vs Saccade Representational Similarity 
            if ~isempty(curr_sac) 
        
                curr_sac_text = saccades_table(strcmp(saccades_table.subject_new,ripple_data.all_subjNames{sb}) & saccades_table.trial_number == it,:);
                curr_sac_text  = table2array(curr_sac_text(~isnan(curr_sac_text.saccade_onset_time) & ~isnan(curr_sac_text.saccade_offset_time),["subject_new","trial_number"]));
        
                R_sac = nan(size(current_sac_rep,2), size(current_sac_rep,2));
                for i = 1:size(current_sac_rep,2)-1
                    for j = i+1:size(current_sac_rep,2)
                        ok = ~isnan(current_sac_rep(:, i)) & ~isnan(current_sac_rep(:, j));
                        if sum(ok) >= 5  % minimalna liczba punktów, ustaw wg. siebie
                            R_sac(i,j) = corr(current_sac_rep(ok, i), current_sac_rep(ok, j), 'type', 'Pearson');
                        else; R_rip(i,j) = NaN; end
                    end
                end
        
        
                [I,J] = ndgrid(1:size(current_sac_rep,2), 1:size(current_sac_rep,2));
                SacIdx1 = I(:); SacIdx2 = J(:);  Rvec = R_sac(:);
                Subject = curr_sac_text(SacIdx1,1);  Trial = curr_sac_text(SacIdx1,2);
                SacOn1 = curr_sac(SacIdx1,1); SacOff1 = curr_sac(SacIdx1,2);
                SacOn2 = curr_sac(SacIdx2,1); SacOff2 = curr_sac(SacIdx2,2);    
        
                T = table(Subject, Trial, SacIdx1,SacIdx2, SacOn1,SacOff1, SacOn2, SacOff2, Rvec, ...
                    'VariableNames',{'subject','trial','sac_idx_1','sac_idx_2','sac_onset_1','sac_offset_1','sac_onset_2','sac_offset_2','R'});
                T = T(~isnan(T.R),:);
        
                T.R_fisher = atanh(min(max(T.R, -0.9999), 0.9999));  T.trial = str2double( T.trial);
                T.rip_sac_timing = T.sac_onset_1 - T.sac_onset_2;
        
                if ~istable(T_sac); T_sac = T([],:); end
                T_sac = [T_sac; T];
            end
    
    
    
            %%% (4) RIpple Fixation Repr Similarity
    
            R = nan(size(current_fix_rep,2), size(current_rip_rep,2),size(current_fix_rep,3));
            for i = 1:size(current_fix_rep,2)
                for j = 1:size(current_rip_rep,2)
                    for m = 1:size(current_fix_rep,3)
    
                        ok = ~isnan(current_fix_rep(:, i,m)) & ~isnan(current_rip_rep(:, j));
        
                        if sum(ok) >= 5  % minimalna liczba punktów, ustaw wg. siebie
                            R(i,j,m) = corr(current_fix_rep(ok, i,m), current_rip_rep(ok, j), 'type', 'Pearson');
                        else; R(i,j,m) = NaN; end
                    end
                end
            end
    
            Rmat = reshape(R, [], size(R,3));     % (fix*rip) rows x windows cols
            Rmat75 = nan(size(Rmat,1), params.fix_box_max);
            Rmat75(:,1:min(size(Rmat,2), params.fix_box_max)) = Rmat(:,1:min(size(Rmat,2), params.fix_box_max));
                    
            [M,J] = ndgrid(1:size(current_fix_rep,2), 1:size(current_rip_rep,2));
            FixIdx = M(:); RipIdx = J(:);
            
            FixBeg = curr_fix(FixIdx,1); FixEnd = curr_fix(FixIdx,2);
            Subject = curr_rip(RipIdx,1); Trial = str2double(curr_rip(RipIdx,3)); Channel = curr_rip(RipIdx,2);
            RipBeg = str2double(curr_rip(RipIdx,5)); RipEnd = str2double(curr_rip(RipIdx,6));
            
            Rnames = "R_" + string(1:params.fix_box_max);
            T = table(Subject, Trial, Channel, FixIdx, RipIdx, FixBeg, FixEnd, RipBeg, RipEnd, ...
                'VariableNames', {'subject','trial','channel','fix_idx','rip_idx','fix_beg','fix_end','rip_beg','rip_end'});
            T.is_excessive = sum(~isnan(Rmat), 2) > params.fix_box_max;
            T.rip_fix_timing_onsets = T.fix_beg - T.rip_end;
            T.rip_fix_timing_middle = (T.fix_beg+ (T.fix_end - T.fix_beg)/2) - (T.rip_end+ (T.rip_end - T.rip_beg)/2);
    
            T = [T, array2table(Rmat75, 'VariableNames', Rnames)];
            T = T(~any(~isnan(Rmat),2)==0, :);   % drop rows all-NaN across windows
            if ~istable(T_fix); T_fix = T([],:); end
            T_fix = [T_fix; T];



            %%% (5) Fixation-Fixation Similarity: (No segments!)

            curr_fix_text = saccades_table(strcmp(saccades_table.subject_new,ripple_data.all_subjNames{sb}) & saccades_table.trial_number == it,:);
            curr_fix_text  = table2array(curr_fix_text(~isnan(curr_fix_text.fixation_beg) & ~isnan(curr_fix_text.fixation_end),["subject_new","trial_number"]));
            curr_fix_text = curr_fix_text(1:end-1,:);

            current_fix_rep2 = squeeze(nanmean(current_fix_rep,3));
            R = nan(size(current_fix_rep2,2), size(current_fix_rep2,2),size(current_fix_rep2,3));
            for i = 1:size(current_fix_rep2,2)
                for j = 1:size(current_fix_rep2,2)
                    ok = ~isnan(current_fix_rep2(:, i)) & ~isnan(current_fix_rep2(:, j));
        
                    if sum(ok) >= 5  % minimalna liczba punktów, ustaw wg. siebie
                        R(i,j) = corr(current_fix_rep2(ok, i), current_fix_rep2(ok, j), 'type', 'Pearson');
                    else; R(i,j) = NaN; end
                end
            end
    
            Rmat = reshape(R, [], size(R,3));     % (fix*rip) rows x windows cols


            
            [I,J] = ndgrid(1:size(current_fix_rep2,2), 1:size(current_fix_rep2,2));
            fix2Idx = I(:); fixIdx = J(:);  Rvec = R(:);
            Subject = curr_fix_text(fixIdx,1);  Trial = str2double(curr_fix_text(fixIdx,2));
            Fix1Beg = curr_fix(fixIdx,1);  Fix1End = curr_fix(fixIdx,2);
            Fix2Beg = curr_fix(fix2Idx,1); Fix2End = curr_fix(fix2Idx,2);
            
            T = table(Subject, Trial,fixIdx,fix2Idx, Fix1Beg, Fix1End, Fix2Beg, Fix2End, Rvec, ...
                'VariableNames',{'subject','trial','fix_idx_1','fix_idx_2','fix_beg_1','fix_end_1','fix_beg_2','fix_end_2','R'});
            T = T(~isnan(T.R),:);

            T.rip_fix_len_1 = T.fix_beg_1 - T.fix_end_1;
            T.rip_fix_len_2 = T.fix_beg_2 - T.fix_end_2;
            T.R_fisher = atanh(min(max(T.R, -0.9999), 0.9999));
    
            if ~istable(T_fix2); T_fix2 = T([],:); end
            T_fix2 = [T_fix2; T];


        end
    end
    if params.bha; reprTitle = "RSA_tables_bha.mat"; else; reprTitle = "RSA_tables.mat"; end

    save(fullfile(saveFolder,params.exact_analysisName,reprTitle),"T_all","T_rip","T_sac","T_fix","T_fix2");
end