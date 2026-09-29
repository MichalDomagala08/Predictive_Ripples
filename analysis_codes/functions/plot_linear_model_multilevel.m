function plot_linear_model_multilevel(current_data, lme, currentDependant, currentPredictor, varNames, plot_errors, plot_alls, binWidth)

    x = current_data.(varNames{1});
    y = current_data.(currentDependant);
    pred_cat = current_data.(currentPredictor);

    % Unikalne poziomy predyktora kategorycznego
    categories = unique(pred_cat(~isnan(pred_cat)));
    nCats = numel(categories);

    % Dynamiczna paleta kolorów (dla 5 poziomów np. lines, turbo lub parula)
    if nCats <= 7
        cmap = lines(nCats);
    else
        cmap = jet(nCats);
    end

    % Punkty siatki predykcji X
    bw = 2; 
    if exist('binWidth','var') && ~isempty(binWidth), bw = binWidth; end
    nPoints = max(50, ceil(range(x(~isnan(x))) / bw));
    x_pred = linspace(min(x(~isnan(x))), max(x(~isnan(x))), nPoints)';

    % Obliczenie binów do error barów
    if plot_errors
        tmin = floor(min(x(~isnan(x)))); 
        tmax = ceil(max(x(~isnan(x)))); 
        edges = tmin:bw:tmax; 
        centers = edges + bw/2;
    end

    hold on;

    % 1. Scatter (surowe punkty)
    if plot_alls
        for c = 1:nCats
            cat_val = categories(c);
            idx_c = (pred_cat == cat_val) & ~isnan(x) & ~isnan(y);
            scatter(x(idx_c), y(idx_c), 14, cmap(c,:), 'filled', 'MarkerFaceAlpha', 0.1, 'HandleVisibility', 'off');
        end
    end

    % 2. Binned Error Bars (Średnia +/- SEM)
    if plot_errors
        for c = 1:nCats
            cat_val = categories(c);
            m = nan(numel(edges), 1);
            s = nan(numel(edges), 1);
            for b = 1:numel(edges)
                inbin = (x >= edges(b)) & (x < edges(b) + bw) & (pred_cat == cat_val);
                y_bin = y(inbin);
                if ~isempty(y_bin) && any(~isnan(y_bin))
                    m(b) = mean(y_bin, 'omitnan');
                    s(b) = std(y_bin, 'omitnan') / sqrt(sum(~isnan(y_bin)));
                end
            end
            valid = ~isnan(m);
            errorbar(centers(valid), m(valid), s(valid), 'o-', ...
                'Color', cmap(c,:), 'LineWidth', 1.2, ...
                'MarkerFaceColor', cmap(c,:), 'MarkerSize', 4, ...
                'HandleVisibility', 'off');
        end
    end

    % 3. Predykcje LME i Wstęgi Przedziałów Ufności (95% CI)
    legHandles = gobjects(nCats, 1);
    legLabels  = cell(nCats, 1);

    for c = 1:nCats
        cat_val = categories(c);
        first_row_idx = find(pred_cat == cat_val, 1);

        % Budowanie tabeli do predykcji
        newVars = cell(1, numel(varNames) - 1);
        for v = 2:numel(varNames)
            var_name = varNames{v};
            if strcmp(var_name, currentPredictor)
                newVars{v-1} = repmat(cat_val, numel(x_pred), 1);
            else
                % Utrzymanie stałej wartości dla pozostałych kowariantów (np. subject)
                newVars{v-1} = repmat(current_data.(var_name)(first_row_idx), numel(x_pred), 1);
            end
        end

        newT = table(x_pred, newVars{:});
        newT.Properties.VariableNames = cellstr(string(varNames));

        % Predykcja z modelu z wyłączeniem efektów losowych per-subject
        try
            [yhat, yci] = predict(lme, newT, 'Conditional', false);
        catch
            [yhat, yci] = predict(lme, newT);
        end

        % Wstęga CI
        fill([x_pred; flipud(x_pred)], [yci(:,1); flipud(yci(:,2))], cmap(c,:), ...
            'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');

        % Linia dopasowania
        legHandles(c) = plot(x_pred, yhat, 'Color', cmap(c,:), 'LineWidth', 2);
        legLabels{c} = sprintf('%s = %s', strrep(currentPredictor, '_', '-'), string(cat_val));
    end

    % 4. Formatowanie wykresu
    xlabel(strrep(varNames{1}, '_', '-'));
    ylabel(strrep(currentDependant, '_', '-'));
    title(sprintf('Model: %s ~ %s', currentDependant, currentPredictor), 'Interpreter', 'none');
    legend(legHandles, legLabels, 'Location', 'best');
    grid on; box on; 
    hold off;
end