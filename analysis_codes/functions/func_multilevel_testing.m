function [curr_est, curr_tStat, curr_p, nam, pred, anname] = func_multilevel_testing( ...
    data, formula_tpl, gr_names, coef_range, anname_val, ...
    group_var, group_labels, pair_testinf, ...
    plot_var, plot_covars, bin_params, ...
    x_label, fig_title_prefix, fig_subtitle_tpl, closeness_label_prefix, ...
    categorical_var, save_ps_path, data_lme)
%FUNC_MULTILEVEL_TESTING  Multilevel LME testing + (annotated multilevel
%boxplot + closeness LM) visualisation.
%
%Consolidates the repeated multilevel block used in PR_5_2 (fix_phase,
%5-level ordinal predictor). For All vs Restricted just CALL TWICE on
%row-filtered data.
%
%ZACHOWANIE vs oryginał: W "Restricted" bloku boxplot (subplot 1) i
%y_max parwise bracketów były rysowane na FULL data, podczas gdy LME,
%closeness i post-plot LME były fitowane na <250ms subset. To idiosynkrazja
%oryginału; tu zachowane przez opcjonalny argument data_lme:
%   data      -> boxplot i y_max bracketów
%   data_lme  -> LME, closeness plot, pairwise aaaa, post-plot LME
%Jeśli data_lme=[] (lub pominięty) → data_lme=data (przypadek "All").
%
%INPUTS
%   data              - tabela używana do boxplota.
%   data_lme          - opcjonalna tabela dla LME/closeness/post-plot LME.
%                       [] = traktuj jako data (domyślnie).
%   formula_tpl       - string formuły z jednym %s na zmienną zależną, np.
%                       "%s ~ fix_phase + seg_to_ripple_distance_mid_abs + (1|subject)"
%   gr_names          - 1xN cell nazw zmiennych zależnych.
%   coef_range        - indeksy współczynników do wyciągnięcia (np. 2:3).
%                       sgtitle i tytuł closeness pokazują coef_range(1).
%   anname_val        - string wpisany w kolumnę anName ("RipProx",
%                       "RipProx_Restr").
%   group_var         - kolumna wielopoziomowa dla boxplota ("fix_phase").
%   group_labels      - 1xN cell etykiet boxplota (N = liczba poziomów,
%                       np. N=5: {'early2','early1','middle','late1','late2'}).
%   pair_testinf      - 1xP cell par indeksów poziomów do post-hoc
%                       overlay na boxplocie, np. {{1,2},{2,3},{3,4},{4,5}}.
%                       [] = wyłącz overlay.
%   plot_var          - predyktor do plot_linear_model_multilevel.
%   plot_covars       - 1xM cell kowariatów do plot_linear_model_multilevel.
%   bin_params        - 1x3 wektor [cov_flag plot_flag n_bins].
%   x_label           - xlabel subplota closeness.
%   fig_title_prefix  - prefiks sgtitle ("All", "Restricted", ...).
%   fig_subtitle_tpl  - sprintf-template z jednym %s (dep_disp) dla title
%                       boxplota, np. "After - Before Boxplot %s".
%   closeness_label_prefix
%                     - prefiks w title closeness, "" lub "Restricted ".
%   categorical_var   - kolumna-categorical-wersja group_var (np.
%                       "fix_phase_cat"); użyta w post-plot LME i do
%                       wstrzykiwania p-values w legendę. [] = wyłącz.
%   save_ps_path      - pełna ścieżka pliku .ps (append).
%
%OUTPUTS
%   curr_est, curr_tStat, curr_p, nam, pred, anname — kolumny do tabeli T.
%
%SEE ALSO: func_binary_testing, fitlme, plot_linear_model_multilevel

    if nargin < 17 || isempty(data_lme)
        data_lme = data;
    end

    curr_est = []; curr_tStat = []; curr_p = []; nam = []; pred = []; anname = [];
    n_coef   = numel(coef_range);
    n_levels = numel(group_labels);
    n_dummies = n_levels - 1;
    disp_idx = coef_range(1);

    for i = 1:numel(gr_names)
        dep      = string(gr_names{i});
        dep_disp = strrep(dep, '_', '-');

        % --- Fit LME (na data_lme) ---------------------------------------
        lme = fitlme(data_lme, sprintf(formula_tpl, dep));

        % --- Wyciągnięcie współczynników ----------------------------------
        curr_est   = [curr_est;   lme.Coefficients.Estimate(coef_range)];
        curr_tStat = [curr_tStat; lme.Coefficients.tStat(coef_range)];
        curr_p     = [curr_p;     lme.Coefficients.pValue(coef_range)];
        pred       = [pred;       string(lme.Coefficients.Name(coef_range))];
        nam        = [nam,        repmat(dep, 1, n_coef)];
        anname     = [anname;     repmat(string(anname_val), n_coef, 1)];

        % --- Figura -------------------------------------------------------
        figure("Visible","Off","PaperOrientation","landscape", ...
               "Units","normalized","Position",[0 0 1 1]);
        sgtitle(sprintf("%s: %s: P: %2.3f T: %2.3f", ...
                        fig_title_prefix, dep_disp, ...
                        lme.Coefficients.pValue(disp_idx), ...
                        lme.Coefficients.tStat(disp_idx)));

        % --- Boxplot + pairwise brackets (na data, brackets na data_lme) --
        subplot(1,2,1); hold on;
        boxplot(data.(dep), data.(group_var), 'Labels', group_labels);

        if ~isempty(pair_testinf)
            y_max = max(data.(dep), [], 'all');
            for j = 1:numel(pair_testinf)
                lv1 = pair_testinf{j}{1};
                lv2 = pair_testinf{j}{2};
                aaaa = data_lme(data_lme.(group_var) == lv1 | ...
                                data_lme.(group_var) == lv2, :);
                lme_t = fitlme(aaaa, sprintf(formula_tpl, dep));
                if lme_t.Coefficients.pValue(disp_idx) < 0.05
                    star = ' *';
                else
                    star = '';
                end
                y = y_max - 0.05;
                plot([j j j+1 j+1], [y y+0.02 y+0.02 y], 'k', 'LineWidth', 1.2);
                text(mean([j j+1]), y+0.05, ...
                     sprintf('p=%.3g%s', lme_t.Coefficients.pValue(disp_idx), star), ...
                     'HorizontalAlignment', 'center');
            end
        end

        ylabel(dep_disp);
        title(sprintf(fig_subtitle_tpl, dep_disp));

        % Kolorowanie boxów
        cmap = lines(n_levels);
        hBox = flipud(findobj(gca,'Tag','Box'));
        hMed = flipud(findobj(gca,'Tag','Median'));
        hOut = flipud(findobj(gca,'Tag','Outliers'));
        for k = 1:min(n_levels, numel(hBox))
            set(hBox(k), 'Color', cmap(k,:));
            set(hMed(k), 'Color', max(0, cmap(k,:) - 0.15));
            if k <= numel(hOut)
                set(hOut(k), 'MarkerEdgeColor', cmap(k,:));
            end
        end

        % --- Closeness multilevel LM (na data_lme) -----------------------
        subplot(1,2,2);
        plot_linear_model_multilevel(data_lme, lme, dep, plot_var, plot_covars, ...
                                     bin_params(1), bin_params(2), bin_params(3));
        xticks(sort([xticks 25]));
        xticklabels(string(str2double(xticklabels) * 2));
        xlabel(x_label);
        title(sprintf("%s%s Closeness:  P: %2.3f T: %2.3f", ...
                      closeness_label_prefix, dep_disp, ...
                      lme.Coefficients.pValue(disp_idx), ...
                      lme.Coefficients.tStat(disp_idx)));

        % --- Post-plot: wstrzyknięcie pairwise p-values do legendy --------
        if ~isempty(categorical_var)
            cat_formula = strrep(formula_tpl, group_var, categorical_var);
            lme_c = fitlme(data_lme, sprintf(cat_formula, dep));
            c_idx = find(contains(lme_c.Coefficients.Name, categorical_var));
            pvals = lme_c.Coefficients.pValue(c_idx);
            lg = legend(gca);
            s = cellstr(lg.String);
            n = numel(s);
            for k = 1:min(numel(pvals), n_dummies)
                pv = pvals(k);
                if ~isnan(pv)
                    s{n-n_dummies+k} = sprintf('%s (p=%.3g)', s{n-n_dummies+k}, pv);
                end
            end
            lg.String = s;
        end

        print(gcf, save_ps_path, '-dpsc', '-append', '-fillpage');
        close(gcf);
    end
end