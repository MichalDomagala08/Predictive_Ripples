function saccades_table_refix = func_compute_refixation2(saccades_table,params)

    
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
    n = height(saccades_table_refix);
    saccades_table_refix.is_refix           = zeros(n,1);
    saccades_table_refix.precursor_idx      = nan(n,1);
    saccades_table_refix.precursor_distance = nan(n,1);
    saccades_table_refix.is_valid_refix = false(n,1);   % pierwsza refiksacja danej kotwicy, niezależnie od lagu
    [G,~,~] = findgroups(saccades_table_refix.subject_new, saccades_table_refix.trial_number);


    saccades_table_refix.local_precursor_idx = nan(n,1);
    saccades_table_refix.local_lag           = nan(n,1);

    for g = 1:max(G)

        %% First run - to have a primary Fixation lists and all
        idxs = find(G==g);
        X = saccades_table_refix.Xpx_square(idxs);
        Y = saccades_table_refix.Ypx_square(idxs);
        valid = ~isnan(X) & ~isnan(Y) & ~isnan(saccades_table_refix.saccade_onset_time(idxs));

        for k = 2:numel(idxs)
            if ~valid(k); continue; end
            for m = 1:k-1                                            % <- bez zmian: od najstarszej
                if valid(m) && ~saccades_table_refix.is_refix(idxs(m))  % <- bez zmian: tylko initial jako kotwica
                    d = hypot(X(k)-X(m), Y(k)-Y(m));
                    if d <= params.p_px
                        saccades_table_refix.is_refix(idxs(k))           = 1;
                        saccades_table_refix.precursor_idx(idxs(k))      = m-1;
                        saccades_table_refix.precursor_distance(idxs(k)) = d;
                        break
                    end
                end
            end
        end


        %% lag "lokalny" liczony w obrębie klastra tej samej kotwicy, a nie względem samej kotwic
        idxs = find(G==g);
        indx_sac_g = saccades_table_refix.indx_sac(idxs);
        is_refix_g = saccades_table_refix.is_refix(idxs) == 1;
        anchor_g   = saccades_table_refix.precursor_idx(idxs);

        % id klastra: kotwica -> jej własny indx_sac; refiksacja -> indx_sac jej kotwicy
        cluster_id = indx_sac_g;
        cluster_id(is_refix_g) = anchor_g(is_refix_g);

        u = unique(cluster_id(~isnan(cluster_id)));
        for c = 1:numel(u)
            members = find(cluster_id == u(c) & ~isnan(indx_sac_g));
            [~,ord] = sort(indx_sac_g(members));
            members = members(ord);              % posortowane chronologicznie w obrębie klastra

            if numel(members) >= 2
                saccades_table_refix.is_valid_refix(idxs(members(2))) = true;   % <-- NOWE: pierwsza refiksacja w łańcuchu
            end

            for j = 2:numel(members)
                cur  = members(j);
                prev = members(j-1);
                saccades_table_refix.local_precursor_idx(idxs(cur)) = indx_sac_g(prev);
                saccades_table_refix.local_precursor_distance(idxs(k)) = hypot(X(k)-X(m), Y(k)-Y(m));
                saccades_table_refix.local_lag(idxs(cur))           = indx_sac_g(cur) - indx_sac_g(prev);
            end
        end
    end

    %%% Lag between the refix and the precursor 
    saccades_table_refix.precursor_lag = saccades_table_refix.indx_sac - saccades_table_refix.precursor_idx;

    %%% Refix Saccades = refiksacje o krótkim LOKALNYM lagu
    saccades_table_refix.is_return_saccade = saccades_table_refix.is_refix & ...
        (saccades_table_refix.local_lag > 1);
end