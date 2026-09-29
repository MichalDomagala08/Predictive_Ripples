function [T_fix_cleaned, T_fix_average, T_fix_tall] = func_RSA_preprocessng_fixation(T_fix,params,masking)

    % Thsi script preprocesses the Fixation RSA tables to separate -
    % analysis usefull tables


    %%% Filtering out Data
    T_fix_cleaned = T_fix;


    % (1) Mask: Get mask for filtering Coincidental Ripples and fixations: 
    segStart = T_fix.fix_beg + (0:params.fix_box_max-1)* params.fix_box_length;           % N x 75, start of each segment
    segEnd   = min(segStart +  params.fix_box_length - 1, T_fix.fix_end); % N x 75, clipped end
    mask_no_overlap_segments = segEnd < T_fix.rip_beg | segStart > T_fix.rip_end;
    
    Rvars = startsWith(T_fix.Properties.VariableNames,"R_");
    Rmat = T_fix_cleaned{:,Rvars};
    Rmat(~mask_no_overlap_segments) = NaN;
    T_fix_cleaned{:,Rvars} = Rmat;
    
    if masking

    
        
        % (2) Mask:  Whether Ripple happened during that trial
        if params.fix_rip_trial 
             mask_no_overlap_fixations = ~(T_fix_cleaned.rip_beg <= T_fix_cleaned.fix_end) & (T_fix_cleaned.rip_end >= T_fix_cleaned.fix_beg);
        else mask_no_overlap_fixations = ones(height(T_fix_cleaned),1); end
    
        % (3) Mask:  for too short fixations
        mask_short_fixations = T_fix_cleaned.fix_end-T_fix_cleaned.fix_beg > params.fix_box_min*params.fix_box_length; % at least 5 segments of length)
        
        T_fix_cleaned = T_fix_cleaned(mask_short_fixations & mask_no_overlap_fixations & ~T_fix_cleaned.is_excessive ,:);
    end

    %%% Compute average R fix
    T_fix_average = T_fix_cleaned(:,1:12);
    T_fix_average.rip_fix_timing_middle_abs = abs(T_fix_average.rip_fix_timing_middle);
    T_fix_average.R = table2array(mean(T_fix_cleaned(:,13:86),2,"omitnan"));
    T_fix_average.R_abs = abs(T_fix_average.R);
    T_fix_average.R_fisher = atanh(min(max(T_fix_average.R, -0.9999), 0.9999)); 
    T_fix_average.R_fisher_abs = abs(T_fix_average.R_fisher);
    
    % Get markers of Before/After:
    T_fix_average.is_before = T_fix_average.fix_end <= T_fix_average.rip_beg -  params.before_cond;
    T_fix_average.is_after  = T_fix_average.fix_beg >= T_fix_average.rip_end +  params.before_cond;                                % if the Rip to Sac timing is negative and less than a boudnary of 5 samples
    T_fix_average.is_before_restricted =  T_fix_average.is_before &  (T_fix_average.rip_beg - T_fix_average.fix_end) <= params.max_time;
    T_fix_average.is_after_restricted = T_fix_average.is_after & (T_fix_average.fix_beg - T_fix_average.rip_end) <= params.max_time;
    
    % Compute absolute distance
    T_fix_average.rip_fix_distance_edge = nan(height(T_fix_average),1);
    T_fix_average.rip_fix_distance_edge( T_fix_average.is_before) = T_fix_average.rip_beg(T_fix_average.is_before) - T_fix_average.fix_end(T_fix_average.is_before);
    T_fix_average.rip_fix_distance_edge(T_fix_average.is_after)   = T_fix_average.fix_beg(T_fix_average.is_after)  - T_fix_average.rip_end(T_fix_average.is_after);
    

    if ismember("is_refix", T_fix.Properties.VariableNames)
        extra_cols = ["Xpx_original","Ypx_original","fixation_len","is_refix", ...
              "precursor_idx","precursor_distance","precursor_lag",...
              'local_precursor_idx','local_lag','local_precursor_distance','is_return_saccade','is_valid_refix'];
        T_fix_average = horzcat(T_fix_average, T_fix_cleaned(:,extra_cols));
    end
    
    %%% Compute Segment-wise Tall Table
    Rvars = T_fix_cleaned.Properties.VariableNames(startsWith(T_fix_cleaned.Properties.VariableNames,'R_'));
    metaVars = T_fix_cleaned.Properties.VariableNames( ~startsWith(T_fix_cleaned.Properties.VariableNames,'R_'));
    
    rowIdx = repelem((1:height(T_fix_cleaned))',numel(Rvars));
    T_fix_tall = T_fix_cleaned(rowIdx,metaVars);
    Rmat = T_fix_cleaned{:,Rvars}; phaseMat = zeros(size(Rmat)); 
    
    % Ger R per segments
    T_fix_tall.segment = repmat((1:numel(Rvars))',height(T_fix_cleaned),1);
    T_fix_tall.R       = reshape(Rmat.',[],1);
    T_fix_tall.R_abs = abs(T_fix_tall.R);
    T_fix_tall.R_fisher = atanh(min(max(T_fix_tall.R, -0.9999), 0.9999)); 
    T_fix_tall.R_fisher_abs = abs(T_fix_tall.R_fisher);
    
    % Ascertainign early,middle,late segments for a given fixation
    for r = 1:size(Rmat, 1)
        validSeg = find(~isnan(Rmat(r, :)));
        nValid   = numel(validSeg);
        if nValid >= params.fix_box_min
            phaseMat(r, validSeg) = ceil((1:nValid) * params.fix_box_min / nValid);
        end
    end
    T_fix_tall.fix_phase = reshape(phaseMat.', [], 1); 
    
    % removing R-segments that do not have a valid R (i.e. are artefact
    T_fix_tall = T_fix_tall(~isnan(T_fix_tall.R), :);
    
    % Get segment tming and distnaces
    T_fix_tall.seg_beg = T_fix_tall.fix_beg + (T_fix_tall.segment-1)*params.fix_box_length;
    T_fix_tall.seg_end = min(T_fix_tall.seg_beg + params.fix_box_length - 1, T_fix_tall.fix_end);
    seg_mid = (T_fix_tall.seg_beg + T_fix_tall.seg_end)/2;
    rip_mid = (T_fix_tall.rip_beg + T_fix_tall.rip_end)/2;
    T_fix_tall.seg_to_ripple_distance_beg = T_fix_tall.seg_beg - T_fix_tall.rip_beg;
    T_fix_tall.seg_to_ripple_distance_mid = seg_mid - rip_mid;
    T_fix_tall.seg_to_ripple_distance_mid_abs = abs(T_fix_tall.seg_to_ripple_distance_mid);

    % Compute segment vs ripple timings
    T_fix_tall.is_before = T_fix_tall.seg_end < T_fix_tall.rip_beg - params.before_cond;
    T_fix_tall.is_after  = T_fix_tall.seg_beg > T_fix_tall.rip_end + params.after_cond;
    T_fix_tall.is_during = ~T_fix_tall.is_before & ~T_fix_tall.is_after;
    T_fix_tall.is_before_restricted = T_fix_tall.is_before & T_fix_tall.rip_beg -T_fix_tall.seg_end <= params.max_time;
    T_fix_tall.is_after_restricted  = T_fix_tall.is_after  & T_fix_tall.seg_beg -T_fix_tall.rip_end <= params.max_time;
    % Removing missng values 
    T_fix_tall = T_fix_tall(~isnan(T_fix_tall.R) & T_fix_tall.fix_phase ~= 0, :);
    
end