function R_unified = func_RSA_preprocessing_nearest_sac_fix(T_fix_org,T_sac,ripple_table,params)

    T_sac.R_abs = abs(T_sac.R);
    T_sac.R_fisher_abs = abs(T_sac.R_fisher);
    
    % Get Markers of Timing for Ripple - Before or after 
    T_sac.is_after  = T_sac.sac_onset   >= T_sac.rip_end + params.before_cond;    
    T_sac.is_before = T_sac.sac_offset  <= T_sac.rip_beg - params.before_cond;                                
    
    
    if ismember("is_refix", T_fix_org.Properties.VariableNames)
        T_fix = T_fix_org(:,[1:12 88:98]);
    else
        T_fix = T_fix_org(:,1:12 );
    end
    T_fix.R = table2array(mean(T_fix_org(:,13:87),2,"omitnan"));
    T_fix.R_abs = abs(T_fix.R);
    T_fix.R_fisher = atanh(min(max(T_fix.R, -0.9999), 0.9999)); 
    T_fix.R_fisher_abs = abs(T_fix.R_fisher);
    
    % Get markers of Before/After:
    T_fix.is_before = T_fix.fix_end <= T_fix.rip_beg -  params.before_cond;
    T_fix.is_after  = T_fix.fix_beg >= T_fix.rip_end +  params.before_cond;                                % if the Rip to Sac timing is negative and less than a boudnary of 5 samples
      
    
    T_sacfix = outerjoin(T_sac,T_fix,"LeftKeys",["subject","trial","channel","rip_idx","sac_offset"],"RightKeys",["subject","trial","channel","rip_idx","fix_beg"]);
    T_sacfix(:, [18 19 20 22 25 26]) = [];
    T_sacfix  = T_sacfix(~isnan(T_sacfix.fix_idx),:);
    T_sacfix = renamevars(T_sacfix, {'subject_T_sac','trial_T_sac','channel_T_sac','rip_idx_T_sac','rip_beg_T_sac','rip_end_T_sac'}, {'subject','trial','channel','rip_idx' 'rip_beg','rip_end'});
    % Overlaping mask
    overlap = (T_sacfix.rip_beg <= T_sacfix.sac_offset) & (T_sacfix.rip_end >= T_sacfix.sac_onset) & ...
              (T_sacfix.rip_beg <= T_sacfix.fix_end)    & (T_sacfix.rip_end >= T_sacfix.fix_beg);       
    
    % Compute absolute distance
    T_sacfix.rip_sac_distance_edge = nan(height(T_sacfix),1);
    T_sacfix.rip_sac_distance_edge(T_sacfix.is_before_T_sac) = T_sacfix.rip_beg(T_sacfix.is_before_T_sac)  - T_sacfix.sac_offset(T_sacfix.is_before_T_sac);
    T_sacfix.rip_sac_distance_edge(T_sacfix.is_after_T_sac)  = T_sacfix.sac_onset(T_sacfix.is_after_T_sac) - T_sacfix.rip_end(T_sacfix.is_after_T_sac);
    
    T_sacfix.rip_fix_distance_edge = nan(height(T_sacfix),1);
    T_sacfix.rip_fix_distance_edge(T_sacfix.is_before_T_fix) = T_sacfix.rip_beg(T_sacfix.is_before_T_fix)  - T_sacfix.fix_end(T_sacfix.is_before_T_fix);
    T_sacfix.rip_fix_distance_edge(T_sacfix.is_after_T_fix)  = T_sacfix.fix_beg(T_sacfix.is_after_T_fix) - T_sacfix.rip_end(T_sacfix.is_after_T_fix);
    
    
    T_NotTheSame = T_sacfix(~overlap,:); % Get ONLY NON-Coincidental ripples (based on params.maxtime)
    
    % Get closest saccade:
    %%% group by Subject x Trial x Ripple Index getting 3 befoe and after
    g = findgroups(T_NotTheSame.subject, T_NotTheSame.trial,T_NotTheSame.channel, T_NotTheSame.rip_idx) ;G = max(g);K = 1;
    nearest_before_sac = nan(G,K);nearest_after_sac  = nan(G,K);
    nearest_before_fix = nan(G,K);nearest_after_fix  = nan(G,K);
    
    % Iterating through group index and getting Before andafer timiing
    for gg = 1:G
        idx = find(g==gg); if isempty(idx), continue; end
        bs = idx(T_NotTheSame.is_before_T_sac(idx)); [~,o]=sort(abs(T_NotTheSame.rip_sac_distance_edge(bs)),'ascend'); mb=min(K,numel(o)); nearest_before_sac(gg,1:mb)=bs(o(1:mb));
        asl= idx(T_NotTheSame.is_after_T_sac(idx));  [~,o]=sort(abs(T_NotTheSame.rip_sac_distance_edge(asl)),'ascend'); ma=min(K,numel(o)); nearest_after_sac(gg,1:ma)=asl(o(1:ma));
        bf = idx(T_NotTheSame.is_before_T_fix(idx)); [~,o]=sort(abs(T_NotTheSame.rip_fix_distance_edge(bf)),'ascend'); mb=min(K,numel(o)); nearest_before_fix(gg,1:mb)=bf(o(1:mb));
        af = idx(T_NotTheSame.is_after_T_fix(idx));  [~,o]=sort(abs(T_NotTheSame.rip_fix_distance_edge(af)),'ascend'); ma=min(K,numel(o)); nearest_after_fix(gg,1:ma)=af(o(1:ma));
    end
        
    
    %%
    %%% (5) Unified Closest-Event Table: 1 row per Ripple × closest Event (saccade OR fixation)
    
    % Recover rip_num via physical key (rip_idx is NOT the true ID after is_repeated filter)
    keyRip = string(ripple_table.subject) + "|" + string(ripple_table.channel) + "|" + string(ripple_table.trial) + "|" + string(ripple_table.beg);
    
    nearests = {nearest_before_sac, nearest_after_sac, nearest_before_fix, nearest_after_fix};
    sides    = ["b" "a" "b" "a"];
    is_fix_v = [false false true true];
    parts    = cell(1,4);
    
    for ci = 1:4
        idxs = nearests{ci};  good = ~isnan(idxs);
        aaa  = T_NotTheSame(idxs(good),:);                                % closest rows
    
        keyA = string(aaa.subject) + "|" + string(aaa.channel) + "|" + string(aaa.trial) + "|" + string(aaa.rip_beg);
        [~, loc] = ismember(keyA, keyRip);  aaa.rip_num = ripple_table.rip_num(loc);
    
        base = table(aaa.subject, aaa.channel, aaa.trial, aaa.rip_num, aaa.rip_beg, aaa.rip_end, ...
                     'VariableNames', {'subject','channel','trial','rip_num','rip_beg','rip_end'});
    
        if is_fix_v(ci)
            ev = table(aaa.fix_idx, aaa.fix_beg, aaa.fix_end, ...
                       'VariableNames', {'event_idx','event_beg','event_end'});
            rs = table(aaa.R_T_fix, aaa.R_abs_T_fix, aaa.R_fisher_T_fix, aaa.R_fisher_abs_T_fix, ...
                       'VariableNames', {'R','R_abs','R_fisher','R_fisher_abs'});
            ds = aaa.rip_fix_distance_edge;
            bf = aaa.is_before_T_fix;   af = aaa.is_after_T_fix;


            if ismember("is_refix", T_fix_org.Properties.VariableNames)
                rf = table(aaa.is_refix, aaa.precursor_idx, aaa.precursor_distance, aaa.precursor_lag, ...
                       'VariableNames', {'is_refix','precursor_idx','precursor_distance','precursor_lag'});
                new = [base ev rs rf];
            else
                new = [base ev rs];
            end

            new.rip_event_distance_edge = ds;
            new.is_before_fix = bf;   new.is_after_fix  = af;   
            new.is_before_sac = false(height(new), 1);
            new.is_after_sac  = false(height(new), 1);
        else
            ev = table(aaa.sac_idx, aaa.sac_onset, aaa.sac_offset, ...
                       'VariableNames', {'event_idx','event_beg','event_end'});
            rs = table(aaa.R_T_sac, aaa.R_abs_T_sac, aaa.R_fisher_T_sac, aaa.R_fisher_abs_T_sac, ...
                       'VariableNames', {'R','R_abs','R_fisher','R_fisher_abs'});

            if ismember("is_refix", T_fix_org.Properties.VariableNames)
                rf = table(nan(height(aaa),1), nan(height(aaa),1),nan(height(aaa),1), nan(height(aaa),1), ...
                       'VariableNames', {'is_refix','precursor_idx','precursor_distance','precursor_lag'});
                new = [base ev rs rf];
            else
                new = [base ev rs];
            end


            ds = aaa.rip_sac_distance_edge;
            bs = aaa.is_before_T_sac;  as = aaa.is_after_T_sac;
            new.rip_event_distance_edge = ds;
            new.is_before_sac = bs;   new.is_after_sac  = as;   
            new.is_before_fix = false(height(new), 1);
            new.is_after_fix  = false(height(new), 1);
        end
        new.is_fix = repmat(is_fix_v(ci), height(new), 1);
        new.side   = repmat(sides(ci),    height(new), 1);
        parts{ci}  = new;
    end
    
    R_unified = vertcat(parts{:});
    
    % Extra RSA / latency metrics (mirror of Tall_cleaned in reference)
    [~,~,ic] = unique(R_unified.subject, 'stable');   R_unified.subject_num = ic;
    R_unified.rank_cat           = categorical(R_unified.side);
    R_unified.event_rip_latency  = R_unified.rip_beg - R_unified.event_beg;
    R_unified = sortrows(R_unified, {'subject','trial','rip_num','is_fix','side'});
    
    % Safety: drop rows where distance wasn't defined (event coincident with ripple, no is_before/is_after)
    R_unified = R_unified(~isnan(R_unified.rip_event_distance_edge),:);
    %%
    absd = abs(R_unified.rip_event_distance_edge);
    grp = findgroups(table(R_unified.subject, R_unified.channel, R_unified.trial, R_unified.rip_num, R_unified.is_fix));
    mins_per_grp = accumarray(grp, absd, [], @min);            % vector length = nGroups
    mins_for_row = mins_per_grp(grp);                         % expand do length(absd)
    R_unified.is_abs_closest = absd == mins_for_row;
    
    R_unified.restricted = abs(R_unified.event_rip_latency)  < params.max_time;
    R_unified.event_rip_latency_abs = abs(R_unified.event_rip_latency);


end