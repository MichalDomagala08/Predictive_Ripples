function plot_linear_model(current_data,lme,currentDependant,currentPredictor,varNames,plot_errors,binWidth)
    % Plots linear model 
    x = current_data.(currentPredictor);
    y = current_data.(currentDependant); 
    if plot_errors
        if binWidth
            edges = min(x)-eps : binWidth : max(x)+eps;
            ic = discretize(x, edges); valid = ~isnan(ic);
            xu = edges(1:end-1); %+ binWidth/2;
            meanY = accumarray(ic(valid), y(valid), [numel(xu),1], @(v) mean(v,'omitnan'), NaN);
            semY  = accumarray(ic(valid), y(valid), [numel(xu),1], @(v) std(v,'omitnan'), NaN);
            valid = ~isnan(meanY);
            xu = xu(valid); meanY = meanY(valid); semY = semY(valid);
        else         
            [xu, ~, ic] = unique(x);                      % xu posortowane rosnąco
            meanY = accumarray(ic, y, [], @nanmean);         % średnie dla każdej unikalnej wartości x
            semY  = accumarray(ic, y, [], @(v) nanstd(v)./sqrt(numel(v))); % SEM (opcjonalnie)
        end
    end
    
    % Predictions
    x_pred = linspace(min(x), max(x), 100)';
    newVars = arrayfun(@(vn) repmat(current_data.(vn)(1), numel(x_pred), 1), string(varNames(2:end)), 'UniformOutput', false);
    newT = table(x_pred, newVars{:});     newT.Properties.VariableNames = cellstr(string(varNames));
    [yhat, yci] = predict(lme, newT);

    coeffs = lme.Coefficients; names = string(coeffs.Name);
    i = find(names==string(varNames{1}),1);

    % rysunek: surowe punkty + uśrednione per-unique x
    hold on;
    scatter(x, y, 14, [0.4 0.4 0.9], 'filled', 'MarkerFaceAlpha', 0.3);     % raw
    if plot_errors
        errorbar(xu, meanY, semY, 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor',[0.4 0.4 0.9]); % binned mean ± SEM
    end
    plot(x_pred, yhat, 'r', 'LineWidth', 2);
    fill([x_pred; flipud(x_pred)], [yci(:,1); flipud(yci(:,2))], 'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none')
    xlabel(currentPredictor);
    ylabel(currentDependant);
    idx = find(string(lme.CoefficientNames) == string(varNames{1}), 1);
    title(sprintf('%s: t=%.2f, p=%.3g', string(varNames{1}), lme.Coefficients.tStat(idx), lme.Coefficients.pValue(idx)), 'Interpreter', 'none');
    
    legend({'Raw observations','Mean per unique TotalRippleNumber'}, 'Location','best');
    grid on; box on; hold off;
end