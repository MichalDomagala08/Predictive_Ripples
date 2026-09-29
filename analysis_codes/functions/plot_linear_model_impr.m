function plot_linear_model_impr(current_data, lme, currentDependant, currentPredictor, varNames, plot_errors, nBins, useJitter)
    if nargin < 8, useJitter = true; end
    nBins = max(2, round(nBins));

    x = current_data.(currentPredictor);
    y = current_data.(currentDependant);
    ok = ~isnan(x) & ~isnan(y);
    x = x(ok); y = y(ok);

    if plot_errors
        edges = unique(quantile(x, linspace(0,1,nBins+1)));   % equal-count bins, robust to skew
        edges(1) = edges(1) - eps; edges(end) = edges(end) + eps;
        ic = discretize(x, edges);
        v = ~isnan(ic);
        meanY = accumarray(ic(v), y(v), [], @(z) mean(z,'omitnan'));
        semY  = accumarray(ic(v), y(v), [], @(z) std(z,'omitnan')./sqrt(numel(z)));
        xu    = accumarray(ic(v), x(v), [], @mean);   % bin center = actual mean x in that bin
    end

    x_pred = linspace(min(x), max(x), 100)';
    holdNames = string(varNames(2:end));
    newVars = cell(1, numel(holdNames));
    for k = 1:numel(holdNames)
        v = current_data.(holdNames(k));
        if isnumeric(v)
            newVars{k} = repmat(mean(v,'omitnan'), numel(x_pred), 1);
        else
            [u,~,ic2] = unique(v);
            [~,mi] = max(accumarray(ic2,1));
            newVars{k} = repmat(u(mi), numel(x_pred), 1);   % keeps original class (cell/string/categorical)
        end
    end
    newT = table(x_pred, newVars{:});
    newT.Properties.VariableNames = cellstr(string(varNames));

    [yhat, yci] = predict(lme, newT, 'Conditional', false);   % population-average curve
    if all(isnan(yhat))
        warning('predict() returned all-NaN — check newT for NaNs or unseen categorical levels.');
    end

    hold on;
    if useJitter
        yj = y + (rand(size(y))-0.5)*0.15;   % display-only jitter, doesn't touch the binned stats
    else
        yj = y;
    end
    scatter(x, yj, 10, [0.4 0.4 0.9], 'filled', 'MarkerFaceAlpha', 0.10);
    if plot_errors
        errorbar(xu, meanY, semY, 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor',[0.4 0.4 0.9]);
    end
    plot(x_pred, yhat, 'r', 'LineWidth', 2);
    fill([x_pred; flipud(x_pred)], [yci(:,1); flipud(yci(:,2))], 'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    xlabel(currentPredictor); ylabel(currentDependant);
    idx = find(string(lme.CoefficientNames) == string(varNames{1}), 1);
    title(sprintf('%s: t=%.2f, p=%.3g', string(varNames{1}), lme.Coefficients.tStat(idx), lme.Coefficients.pValue(idx)), 'Interpreter','none');
    legend({'Raw (jittered)','Binned mean ± SEM','Model fit (population avg.)'}, 'Location','best');
    grid on; box on; hold off;
end