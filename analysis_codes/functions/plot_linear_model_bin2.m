function plot_linear_model_bin2(current_data,lme,currentDependant,currentPredictor,varNames,plot_errors,plot_alls,binWidth)

    % In this version of this plotting function, we have more 
    x = current_data.(varNames{1});
    y = current_data.(currentDependant); 
    
    % Predictions 1:
    x_pred = linspace(0, max(x),length(x)/binWidth)';
    newVars = arrayfun(@(vn) repmat(current_data.(vn)(find(~current_data.(varNames{2}),1)), numel(x_pred), 1), string(varNames(2:end)), 'UniformOutput', false);
    newT = table(x_pred, newVars{:});     newT.Properties.VariableNames = cellstr(string(varNames));
    [yhat, yci] = predict(lme, newT);

    % Predictions 2:
    newVars = arrayfun(@(vn) repmat(current_data.(vn)(find(current_data.(varNames{2}),1)), numel(x_pred), 1), string(varNames(2:end)), 'UniformOutput', false);
    newT = table(x_pred, newVars{:});     newT.Properties.VariableNames = cellstr(string(varNames));
    [yhat2, yci2] = predict(lme, newT);

    coeffs = lme.Coefficients; names = string(coeffs.Name);

    hold on;
    if plot_alls
        scatter(x(current_data.(currentPredictor) == 0), y(current_data.(currentPredictor) == 0), 14, [0 0.45 0.74], 'filled', 'MarkerFaceAlpha', 0.1);
        scatter(x(current_data.(currentPredictor) > 0), y(current_data.(currentPredictor) > 0), 14, [0.85 0.33 0.10], 'filled', 'MarkerFaceAlpha', 0.1);
    end

    if plot_errors
        bw = 2; if exist('binWidth','var') && ~isempty(binWidth), bw = binWidth; end
        tmin = floor(min(x)); tmax = ceil(max(x)); edges = tmin:bw:tmax; centers = edges + bw/2;
        m0 = nan(numel(edges),1); s0 = m0; m1 = m0; s1 = m0;
        for i=1:numel(edges)
            inbin = x>=edges(i) & x<edges(i)+bw;
            y0 = y(inbin & current_data.(currentPredictor) == 0); y1 = y(inbin & current_data.(currentPredictor) > 0);
            m0(i)=mean(y0,'omitnan'); s0(i)=std(y0,'omitnan')/sqrt(max(1,numel(y0)));
            m1(i)=mean(y1,'omitnan'); s1(i)=std(y1,'omitnan')/sqrt(max(1,numel(y1)));
        end
        valid = ~isnan(m0) | ~isnan(m1);
        errorbar(centers(valid), m0(valid), s0(valid), 'o-', 'Color',[0 0.45 0.74], 'LineWidth',1.2,'MarkerFaceColor',[0 0.45 0.74]);
        errorbar(centers(valid), m1(valid), s1(valid), 's-', 'Color',[0.85 0.33 0.10], 'LineWidth',1.2,'MarkerFaceColor',[0.85 0.33 0.10]);
    end



    % draw model predictions as before
    plot(x_pred, yhat, 'b', 'LineWidth', 2);
    plot(x_pred, yhat2, 'r', 'LineWidth', 2);
    fill([x_pred; flipud(x_pred)], [yci(:,1); flipud(yci(:,2))], 'b', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    fill([x_pred; flipud(x_pred)], [yci2(:,1); flipud(yci2(:,2))], 'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none');

    % Labels and all:
    idx = find(contains(string(lme.CoefficientNames), string(varNames{2})), 1);
    title(sprintf('%s ~ %s: t=%.2f, p=%.3g',currentDependant, string(varNames{2}), lme.Coefficients.tStat(idx), lme.Coefficients.pValue(idx)), 'Interpreter', 'none');
    xlabel(strrep(varNames{1},'_','-'));ylabel(strrep(currentDependant,'_','-'));
    legend({sprintf(' %s False',strrep(currentPredictor,'_','-')),sprintf('%s True',strrep(currentPredictor,'_','-')),}, 'Location','best');
    grid on; box on; hold off;
end