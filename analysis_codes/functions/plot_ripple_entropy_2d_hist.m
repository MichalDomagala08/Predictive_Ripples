function plot_ripple_entropy_2d_hist(currentEntropy,ripplewise_ent)
    
    
    % time bins: 1500:100:5500 (edges)
    time_edges = 1500:100:5500;          % -> 40 bins between 1500-1600...5400-5500
    time_centers = time_edges(1:end-1) + diff(time_edges)/2;
    
    % get vectors
    tvec = ripplewise_ent.median;    % assumed numeric
    awsvec = ripplewise_ent.(currentEntropy); % assumed numeric (if different col name, change)
    
    
    % entropy (AWS) bins: 0:0.1:2
    aws_edges = floor(min(awsvec)):0.1:ceil(max(awsvec));                 % -> 20 bins
    aws_centers = aws_edges(1:end-1) + diff(aws_edges)/2;
    
    
    % remove NaNs in either
    valid = ~isnan(tvec) & ~isnan(awsvec);
    tvec = tvec(valid);
    awsvec = awsvec(valid);
    
    % Digitize into bins (returns bin index: 1..Nbins, or NaN if out of range)
    tbin = discretize(tvec, time_edges);
    awsbin = discretize(awsvec, aws_edges);
    
    % Build 2D count matrix
    nb_time = length(time_edges)-1;
    nb_aws = length(aws_edges)-1;
    counts2d = zeros(nb_aws, nb_time);
    
    % accumulate counts
    for k = 1:length(tbin)
        ti = tbin(k);
        ai = awsbin(k);
        if ~isnan(ti) && ~isnan(ai)
            counts2d(ai, ti) = counts2d(ai, ti) + 1;
        end
    end
    
    
    % Optionally normalize per column or overall:
    counts2d_norm = counts2d; % keep raw
    % normalize across time bins (per column) example:
    counts2d_colnorm = counts2d ./ (sum(counts2d,1) + eps);

    imagesc((time_centers-2000) *2, aws_centers, counts2d_colnorm);   % x=czasy, y=aws
    axis xy; colormap(parula); colorbar; hold on;
    xline(0, 'w', 'LineWidth', 3);    % teraz x=2000 jest w tej samej skali co time_centers
    xline(6000, 'w', 'LineWidth', 3);

end