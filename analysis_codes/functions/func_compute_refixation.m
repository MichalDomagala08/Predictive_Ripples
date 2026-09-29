
function saccades_table_refix = func_compute_refixation(saccades_table,params)

    % This funciton computes re-fixation using simple thresholding
    % parameter. 
    %
    % It computes euclidean distance between X and Y eye position after the
    % saccades for every saccade pair in a given trial, and then thresholds it
    %
    % The saccades are written up incrementaly, thus allowing for first
    % stablishing initial saccades to location and then focusing only on
    % pairs to said locations.
    %
    % The code also inputs the Lag distance between the precursor and a
    % current fixation, as well as a distance
    %
    % Main return is is_refix - a binary variable tracking whether the
    % Current fixation is a refixation 

    saccades_table_refix = saccades_table;
    
    saccades_table_refix.is_refix = zeros(height(saccades_table_refix),1);
    saccades_table_refix.precursor_idx = nan(height(saccades_table_refix),1);
    saccades_table_refix.precursor_distance =  nan(height(saccades_table_refix),1);
    [G,~,~] = findgroups(saccades_table_refix.subject_new, saccades_table_refix.trial_number);
    for g = 1:max(G)
        idxs = find(G==g);
        X = saccades_table_refix.Xpx_square(idxs);
        Y = saccades_table_refix.Ypx_square(idxs);
        % saccades_table_refix.fixation_mean_X(idxs) = 
        % saccades_table_refix.fixation_mean_Y(idxs) = 

        valid = ~isnan(X) & ~isnan(Y) & ~isnan(saccades_table_refix.saccade_onset_time(idxs));
        for k = 2:numel(idxs)
            if ~valid(k); continue; end
            for m = 1:k-1% Iterating thoruh all th other 
                if valid(m) && ~saccades_table_refix.is_refix(idxs(m))  % tylko initial jako kotwica
                    if hypot(X(k)-X(m), Y(k)-Y(m)) <= params.p_px
                        saccades_table_refix.is_refix(idxs(k)) = 1;
                        saccades_table_refix.precursor_idx(idxs(k)) = m-1;
                        saccades_table_refix.precursor_distance(idxs(k)) =  hypot(X(k)-X(m), Y(k)-Y(m));
                        break
                    end
                end
            end
        end
    end
    
    saccades_table_refix.precursor_lag = saccades_table_refix.indx_sac - saccades_table_refix.precursor_idx; 
end
