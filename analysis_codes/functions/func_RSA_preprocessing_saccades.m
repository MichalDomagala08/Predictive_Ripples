function [overlap,T_all_cleaned,Tall_cleaned,R_table]=func_RSA_preprocessing_saccades(T_all,ripple_table,params)

    T_all_cleaned = T_all;



    %%% (2) Clean Saccade x Ripple table
    
    % Compute FIsher and other Representational Similarity metrics
    T_all_cleaned.R_abs = abs(T_all_cleaned.R);
    T_all_cleaned.R_fisher_abs = abs(T_all_cleaned.R_fisher);
    T_all_cleaned.rip_sac_timing_abs = abs(T_all_cleaned.rip_sac_timing);
    T_all_cleaned.R_abs_log = log(abs(T_all_cleaned.R));
    T_all_cleaned.R_fisher_abs_log = log(abs(T_all_cleaned.R_fisher));
    T_all_cleaned.rip_sac_timing_abs_log = log(abs(T_all_cleaned.rip_sac_timing));
    
    % Get Markers of Timing for Ripple - Before or after 
    T_all_cleaned.is_after  = T_all_cleaned.sac_onset   >= T_all_cleaned.rip_end + params.before_cond;    
    T_all_cleaned.is_before = T_all_cleaned.sac_offset  <= T_all_cleaned.rip_beg - params.before_cond;                                
    T_all_cleaned.is_after_restricted =  T_all_cleaned.is_after  & (T_all_cleaned.sac_onset - T_all_cleaned.rip_end) <= params.max_time;
    T_all_cleaned.is_before_restricted=  T_all_cleaned.is_before & (T_all_cleaned.rip_beg - T_all_cleaned.sac_offset) <= params.max_time;
    
    % Compute absolute distance
    T_all_cleaned.rip_sac_distance_edge = nan(height(T_all_cleaned),1);
    T_all_cleaned.rip_sac_distance_edge(T_all_cleaned.is_before) = T_all_cleaned.rip_beg(T_all_cleaned.is_before)  - T_all_cleaned.sac_offset(T_all_cleaned.is_before);
    T_all_cleaned.rip_sac_distance_edge(T_all_cleaned.is_after)  = T_all_cleaned.sac_onset(T_all_cleaned.is_after) - T_all_cleaned.rip_end(T_all_cleaned.is_after);

    %%% (1) Masking
    
    % get mask for fitlering ripple & saccade coincidence
    if isnumeric(params.sac_rip_coincidence)
        overlap =  (T_all_cleaned.rip_sac_distance_edge < params.sac_rip_coincidence & T_all_cleaned.rip_sac_distance_edge >0);
    else
        overlap = (T_all_cleaned.rip_beg <= T_all_cleaned.sac_offset) & (T_all_cleaned.rip_end >= T_all_cleaned.sac_onset);
    end
    

    %%% (3) Finding closest saccades

    T_NotTheSame = T_all_cleaned(~overlap,:); % Get ONLY NON-Coincidental ripples (based on params.maxtime)
    
    %%% group by Subject x Trial x Ripple Index getting 3 befoe and after
    g = findgroups(T_NotTheSame.subject, T_NotTheSame.trial,T_NotTheSame.channel, T_NotTheSame.rip_idx) ;G = max(g);K = 3;
    nearest_before = nan(G,K);nearest_after  = nan(G,K);
    
    % Iterating through group index and getting Before andafer timiing
    for gg = 1:G
        idx = find(g==gg); if isempty(idx), continue; end
        d = T_NotTheSame.rip_sac_distance_edge(idx);                        % signed timing
        bidx = T_NotTheSame(idx,:).is_before; 
        if ~isempty(bidx), [~,ob] = sort(abs(d(bidx)),'ascend'); mb=min(K,numel(ob)); nearest_before(gg,1:mb)= idx(ob(1:mb)); end
        aidx = T_NotTheSame(idx,:).is_after; 
        if ~isempty(aidx), [~,oa] = sort(abs(d(aidx)),'ascend'); ma=min(K,numel(oa)); nearest_after(gg,1:ma)  = idx(oa(1:ma)); end
    end
    
    
    % Getting separate columns for before up to K saccades 
    R_table = ripple_table(:,1:6);
    for it = 1:K
        aaa = T_NotTheSame(nearest_before(~isnan(nearest_before(:,it)),it),:);
        
        % Repair rip_num using physical timestamp (rip_idx is NOT the true ripple ID after is_repeated filter)
        keyA = string(aaa.subject) + "|" + string(aaa.channel) + "|" + string(aaa.trial) + "|" + string(aaa.rip_beg);
        keyB = string(ripple_table.subject) + "|" + string(ripple_table.channel) + "|" + string(ripple_table.trial) + "|" + string(ripple_table.beg);
        [~, loc] = ismember(keyA, keyB);
        aaa.rip_num = ripple_table.rip_num(loc);
        
        % Select columns and rename to b1/b2/b3 BEFORE join
        aaa = aaa(:,{'subject','channel','trial','rip_num','sac_idx','sac_onset','R','R_fisher'});
        idx = ~ismember(aaa.Properties.VariableNames, {'subject','channel','trial','rip_num'});
        aaa.Properties.VariableNames(idx) = strcat(aaa.Properties.VariableNames(idx), sprintf('_b%d',it));
        
        R_table = outerjoin(R_table,aaa, ...
                  'LeftKeys',  {'subject','channel','trial','rip_num'}, ...
                  'RightKeys', {'subject','channel','trial','rip_num'}, ...
                  'Type', 'left', 'MergeKeys', true);
    end
    
    % Getting separate columns for after up to K saccades
    for it = 1:K
        aaa = T_NotTheSame(nearest_after(~isnan(nearest_after(:,it)),it),:);
        
        keyA = string(aaa.subject) + "|" + string(aaa.channel) + "|" + string(aaa.trial) + "|" + string(aaa.rip_beg);
        keyB = string(ripple_table.subject) + "|" + string(ripple_table.channel) + "|" + string(ripple_table.trial) + "|" + string(ripple_table.beg);
        [~, loc] = ismember(keyA, keyB);
        aaa.rip_num = ripple_table.rip_num(loc);
        
        aaa = aaa(:,{'subject','channel','trial','rip_num','sac_idx','sac_onset','R','R_fisher'});
        idx = ~ismember(aaa.Properties.VariableNames, {'subject','channel','trial','rip_num'});
        aaa.Properties.VariableNames(idx) = strcat(aaa.Properties.VariableNames(idx), sprintf('_a%d',it));
        
        R_table = outerjoin(R_table,aaa, ...
                  'LeftKeys',  {'subject','channel','trial','rip_num'}, ...
                  'RightKeys', {'subject','channel','trial','rip_num'}, ...
                  'Type', 'left', 'MergeKeys', true);
    end
    
    for i = 1:6
        R_table.(strcat(R_table.Properties.VariableNames{9+4*(i-1)},'_abs'))      =  table2array(abs(R_table(:,9+4*(i-1))));
        R_table.(strcat(R_table.Properties.VariableNames{9+4*(i-1)},'_abs_log'))  =  table2array(log(abs(R_table(:,9+4*(i-1)))));
        R_table.(strcat(R_table.Properties.VariableNames{10+4*(i-1)},'_abs'))     =  table2array(abs(R_table(:,10+4*(i-1))));
        R_table.(strcat(R_table.Properties.VariableNames{10+4*(i-1)},'_abs_log')) =  table2array(log(abs(R_table(:,10+4*(i-1)))));
    end
    R_table = sortrows(R_table, {'subject','trial','rip_num'});
    
    
    %%% (4) Create Tall Table of closest saccades  
    % Create Tall structure of Before/After and Distance (n. of saccades) 
    base = R_table(:,{'subject','channel','trial','rip_num','beg','end'}); vars = ["sac_idx","sac_onset","R","R_fisher"];
    combos = ["b1","b2","b3","a1","a2","a3"]; sides = ["b","b","b","a","a","a"]; ranks = [1 2 3 1 2 3];
    parts = cell(1,numel(combos));
    for ci=1:numel(combos)
        cols = vars + "_" + combos(ci); new = base;
        for v=1:numel(vars)
            if ismember(cols(v), R_table.Properties.VariableNames), new.(vars(v)) = R_table{:,cols(v)};
            else new.(vars(v)) = nan(height(R_table),1); end
        end
        new.side = repmat(sides(ci), height(R_table), 1); new.rank = repmat(ranks(ci), height(R_table), 1);
        parts{ci} = new;
    end
    Tall = vertcat(parts{:});
    
    
    %%% Create cleaned (without NANs) table with additional columns
    Tall_cleaned = Tall(~isnan(Tall.R),:);
    [~, ~, ic] = unique(Tall_cleaned.subject, 'stable'); Tall_cleaned.subject_num = ic;
    Tall_cleaned.R_abs = abs(Tall_cleaned.R);
    Tall_cleaned.R_abs_log = log(abs(Tall_cleaned.R));
    Tall_cleaned.R_fisher_abs = abs(Tall_cleaned.R_fisher);
    Tall_cleaned.R_fisher_abs_log = log(abs(Tall_cleaned.R_fisher));
    Tall_cleaned.rank_cat = categorical(Tall_cleaned.rank);
    Tall_cleaned.sac_rip_latency = Tall_cleaned.beg -   Tall_cleaned.sac_onset;
    
    if params.sac_rip_beginning
        Tall_cleaned = Tall_cleaned(Tall_cleaned.beg > 2000+params.sac_rip_beginning/2,:);
    end





end