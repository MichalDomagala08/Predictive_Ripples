

function plot_diagnostic_lm(Out, model, yname,xname)
% diagnosticLMEPlot(Out, model, yname)
% Out - tabela z danymi (z kolumną xname)
% model - dopasowany model (fitlme / fitlm / fitglme)
% yname - nazwa zmiennej zależnej (string/char)

% ---- extract residuals & fitted robustnie ----
if isa(model, 'LinearMixedModel') || isa(model, 'GeneralizedLinearMixedModel')
    r = residuals(model);
    f = fitted(model);
    modelVars = model.Variables;
elseif isa(model, 'LinearModel')
    if isprop(model, 'Residuals')
        r = model.Residuals.Raw;
    else
        r = residuals(model);
    end
    if isprop(model, 'Fitted')
        f = model.Fitted;
    else
        f = predict(model, model.Variables);
    end
    modelVars = model.Variables;
else
    error('Unsupported model class: %s', class(model));
end

r_std = (r - mean(r)) ./ std(r);
x = Out.(xname);

% tests / summaries
try
    [h_lillie, p_lillie] = lillietest(r);
catch
    h_lillie = NaN; p_lillie = NaN;
end


% plotting

% 1 Residuals vs Fitted (LOESS on sorted f)
subplot(3,3,1);
scatter(f, r, 12, [0.2 0.6 1], 'filled'); hold on;
[f_sort, ifs] = sort(f);
r_smooth_f = smooth(f_sort, r(ifs), 0.2, 'loess');        % same length as f_sort
plot(f_sort, r_smooth_f, 'r-', 'LineWidth', 1.2);
yline(0,'k--'); xlabel('Fitted'); ylabel('Residuals'); title('Residuals vs Fitted'); box on; hold off;

% 2 Residuals vs Predictor (LOESS on sorted x)
subplot(3,3,2);
scatter(x, r, 12, [0.2 0.6 1], 'filled'); hold on;
[x_sort, ix] = sort(x);
r_smooth_x = smooth(x_sort, r(ix), 0.2, 'loess');
plot(x_sort, r_smooth_x, 'r-', 'LineWidth', 1.2);
yline(0,'k--'); xlabel(xname); ylabel('Residuals'); title('Residuals vs Predictor'); box on; hold off;

% 3 QQ
subplot(3,3,3);
qqplot(r); title(sprintf('QQ plot (Lillie p=%.3g)', p_lillie));

% 4 Hist
subplot(3,3,4);
histogram(r, 40, 'FaceColor',[0.2 0.6 1],'EdgeColor','none');
xlabel('Residuals'); title(sprintf('Hist (skew=%.2f kurt=%.2f)', skewness(r), kurtosis(r))); box on;

% 5 Box by quartile
% subplot(3,3,5);
% try
% boxplot(r, grp, 'Labels', {'Q1','Q2','Q3','Q4'});
% xlabel(sprintf('%s quartiles',xname)); ylabel('Residuals'); title('Residuals by X quartile');
% catch; end

subplot(3,3,5);
outlier_idx = find(abs(r_std) > 3);
dispersion = var(Out.(yname)) / mean(Out.(yname));

try
    grp = discretize(x, 4);
    if numel(r) < 20
        axis off; text(0.1,0.5,'Too few obs for quartiles','FontSize',11);
    else
        grpq = discretize(x,4);
        counts = accumarray(grpq(~isnan(grpq)),1,[4,1]);
        if any(counts <10)
            % brak dostatecznej liczby w którymś kwartylu -> pokaż surowe punkty i informację
            scatter(x, r, 10, [0.2 0.6 1], 'filled'); xlabel(xname); ylabel('Residuals');
            title(sprintf('Residuals (not enough obs for quartiles, minPerBin=10)'));
        else
            boxplot(r, grpq, 'Labels', {'Q1','Q2','Q3','Q4'});
            xlabel(sprintf('%s quartiles', xname)); ylabel('Residuals'); title('Residuals by X quartile');
        end
    end
catch;end;

% 6 Autocorrelation of residuals (toolbox-free)
subplot(3,3,6);
maxlag = min(50, numel(r)-1);
N = numel(r);
mu = mean(r);
% oblicz ACF dla lag = 0..maxlag
acf = zeros(maxlag+1,1);
den = sum((r - mu).^2);
for k = 0:maxlag
    acf(k+1) = sum((r(1:N-k)-mu) .* (r(1+k:N)-mu)) / den;
end
stem(0:maxlag, acf, 'filled');
xlabel('Lag'); ylabel('ACF'); title('Autocorrelation (ACF)');
grid on;
% przyblizone 95% CI dla ACF to +/- 1.96/sqrt(N)
ci = 1.96 / sqrt(N);
hold on;
plot([0 maxlag], [ci ci], 'r--', 'LineWidth', 1);
plot([0 maxlag], [-ci -ci], 'r--', 'LineWidth', 1);
hold off;

% 7 Top residuals table (text fallback)
subplot(3,3,7);
tblOut = table((1:numel(r))', r, r_std, 'VariableNames', {'Idx','Residual','Res_z'});
[~, ord] = sort(abs(tblOut.Res_z), 'descend');
tblShow = tblOut(ord(1:min(15,height(tblOut))), :);
txt = evalc('disp(removevars(tblShow,{}))'); 
txt = regexprep(txt,'<[^>]*>','');   % usuń znaczniki HTML/strong jeśli się pojawią
text(0,0.95,txt,'FontName','Courier','FontSize',9,'Interpreter','none','VerticalAlignment','top');
axis off;
title(sprintf('Top %d residuals (|z|>3 => %d)', min(15,height(tblOut)), numel(outlier_idx)));

% 8 Model summary & checks
try
    subplot(3,3,8);
    axis off;
    txt = sprintf('Model class: %s\n', class(model));
    if isprop(model, 'Coefficients')
        coeffs = model.Coefficients;
        for k = 1:height(coeffs)
            txt = [txt, sprintf('%s: Est=%.4g  SE=%.4g  p=%.3g\n', coeffs.Name{k}, coeffs.Estimate(k), coeffs.SE(k), coeffs.pValue(k))];
        end
    end
    try
        [p_levene, ~] = vartestn(r, grp, 'TestType','LeveneAbsolute','Display','off');
        txt = [txt, sprintf('\nLillie p=%.3g\nLevene p (quartiles): %.3g\n', p_lillie, p_levene)];
    catch
        txt = [txt, sprintf('\nLillie p=%.3g\nLevene p (quartiles): n/a\n', p_lillie)];
    end
    txt = [txt, sprintf('Outliers (|z|>3): %d\nDispersion (var/mean) of Y: %.2f\n', numel(outlier_idx), dispersion)];
    text(0,0.95, txt, 'FontName','Courier','FontSize',10, 'VerticalAlignment','top');
    title('Model summary & checks');
catch;end;
% 9 Local trend + fit overlays (LOESS on x)
subplot(3,3,9);
scatter(x, Out.(yname), 18, [0.6 0.6 0.6], 'filled'); hold on;
[xs_sort, idxs] = sort(x);
y_smooth = smooth(xs_sort, Out.(yname)(idxs), 0.2, 'loess');
plot(xs_sort, y_smooth, ':', 'Color',[0.5 0 0.7], 'LineWidth', 1.6);

% fitlm line (if possible)
try
    lm = fitlm(x, Out.(yname));
    y_lm = predict(lm, xs_sort');
    plot(xs_sort, y_lm, 'k-', 'LineWidth', 1.6);
    legendEntries = {'Data','LOESS','fitlm'};
catch
    legendEntries = {'Data','LOESS'};
end

% model fixed-effects population prediction (robust)
nGrid = numel(xs_sort);
try
    baseRow = model.Variables(1,:);
    newT_model = repmat(baseRow, nGrid, 1);
    newT_model.(xname) = xs_sort;

    try
        y_model = predict(model, newT_model, 'Conditional', false);
    catch
        y_model = predict(model, newT_model);
    end
    plot(xs_sort, y_model, 'r-', 'LineWidth', 1.6);
    legendEntries{end+1} = 'model fixed';
catch
    % skip
end

xlabel(xname); ylabel(yname); title('Data + LOESS + fitlm + model-fixed');
legend(legendEntries, 'Location','best'); box on; hold off;

end