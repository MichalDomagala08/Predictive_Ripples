function [curr_est, curr_tStat, curr_p, nam, pred, anname] = func_binary_testing2( ...
    data, formula_tpl, gr_names, coef_range, anname_val, lme_par,...
    group_var, group_labels, plot_covars, bin_params, ...
    x_label, fig_title_prefix, fig_subtitle, save_ps_path)
    %FUNC_BINARY_TESTING  Binary LME testing + (boxplot + closeness LM) visualisation.
    %
    %Consolidates the repeated binary-test block used across PR_5_2:
    %reversed vs not, is_fix vs not, is_after vs not, etc.
    %For Restricted vs All just CALL TWICE — restricted is simply the same
    %call on a row-filtered table with different bin_params / x_label /
    %anname_val / fig_title_prefix.
    %
    %INPUTS
    %   data              - tabela na której fitujemy LME (już przefiltrowana,
    %                       np. R_unified_reverse_fix z maską
    %                       is_abs_closest & rip_event_distance_edge>0).
    %   formula_tpl       - string formuły z jednym %s na zmienną zależną, np.
    %                       "%s ~ reversed + event_rip_latency_abs + (1|subject) + (1|rip_num)"
    %   gr_names          - 1xN cell nazw zmiennych zależnych, np.
    %                       {'R','R_abs','R_fisher','R_fisher_abs'}.
    %   coef_range        - indeksy współczynników do wyciągnięcia, np. 2:3.
    %                       Każdy gr_name dokłada tyle wierszy ile tu jest.
    %   anname_val        - string wpisany w kolumnę anName dla każdego wiersza
    %                       z tego wywołania, np. "RipProx" lub "RipProx_Restr".
    %   lme_par           - parametr określający typ rozkładu: 
    %                        ""         oznacza zwyczajny
    %                        "Poisson"  oznacza rozkład Poissona
    %                        "Gamma"    oznacza Gamma log link
    %   group_var         - nazwa kolumny binarnej dla boxplota, np. "reversed",
    %                       "is_fix", "is_after_restricted".
    %   group_labels      - 1x2 cell etykiet boxplota, np.
    %                       {'Fixation','Not Fixation'} lub {'Saccade','Fixation'}.
    %   plot_var          - nazwa predyktora przekazywana do plot_linear_model_bin2
    %                       (4. arg pozycyjny), np. "reversed".
    %   plot_covars       - 1xM cell kowariatów do plot_linear_model_bin2 (5. arg).
    %   bin_params        - 1x3 wektor [cov_flag plot_flag n_bins] — końcówka
    %                       pozycyjna plot_linear_model_bin2, np. [1 1 250] lub
    %                       [1 1 20].
    %   x_label           - xlabel dla subplota closeness, np.
    %                       "Absolute Saccade to Ripple distance (ms)".
    %   fig_title_prefix  - prefiks w sgtitle, np. "All", "Restricted",
    %                       sprintf("Restricted (+/- %d ms)", params.max_time/2).
    %   fig_subtitle      - napis używany jako title boxplota, np.
    %                       "Fixation vs Reverse to Ripple RSA".
    %   save_ps_path      - pełna ścieżka pliku .ps do którego appendujemy figury.
    %
    %OUTPUTS
    %   curr_est, curr_tStat, curr_p, nam, pred, anname
    %                     - kolumny gotowe do złożenia w tabelę T i zapisu do
    %                       Excela site-owo (store + writetable zostaje u Ciebie).
    %
    %SEE ALSO: fitlme, plot_linear_model_bin2

    curr_est = []; curr_tStat = []; curr_p = []; nam = []; pred = []; anname = [];
    n_coef = numel(coef_range);

    for i = 1:numel(gr_names)
        dep      = string(gr_names{i});
        dep_disp = strrep(dep, '_', '-');

        % --- Fit LME ------------------------------------------------------
        switch lme_par{i}
            case "";        lme = fitlme(data, sprintf(formula_tpl, dep));
            case "poisson"
                lme = fitglme(data, sprintf(formula_tpl, dep), 'Distribution','Poisson','Link', 'log','Offset', log(data.(plot_covars{1})));
           case "binom"
                lme = fitglme(data, sprintf(formula_tpl, dep), "Distribution", "binomial", "Link", "logit");
            case "gamma"  
                lme = fitglme(data, sprintf(formula_tpl, dep), 'Distribution','Gamma', 'Link','log');
        end
            

        % --- Wyciągnięcie współczynników ----------------------------------
        curr_est   = [curr_est;   lme.Coefficients.Estimate(coef_range)];
        curr_tStat = [curr_tStat; lme.Coefficients.tStat(coef_range)];
        curr_p     = [curr_p;     lme.Coefficients.pValue(coef_range)];
        pred       = [pred;       string(lme.Coefficients.Name(coef_range))];
        nam        = [nam,        repmat(dep, 1, n_coef)];
        anname     = [anname;     repmat(string(anname_val), n_coef, 1)];
        main_idx = find(contains(lme.Coefficients.Name, group_var), 1);

        % --- Figura -------------------------------------------------------
        figure("Visible","Off","PaperOrientation","landscape", ...
               "Units","normalized","Position",[0 0 1 1]);
        sgtitle(sprintf("%s: %s: P: %2.3f T: %2.3f", fig_title_prefix, dep_disp,  lme.Coefficients.pValue(main_idx),  lme.Coefficients.tStat(main_idx)));

        % Boxplot
        subplot(1,2,1);
        try
            boxplot(data.(dep), data.(group_var), 'Labels', group_labels);
        catch; end
        ylabel(dep_disp);
        title(fig_subtitle);

        % Closeness linear model
        subplot(1,2,2);
        plot_linear_model_bin2(data, lme, dep, group_var, plot_covars,  bin_params(1), bin_params(2), bin_params(3));

        xticklabels(string(str2double(xticklabels) * 2));
        xlabel(x_label);
        title("Closeness linear Model: ");

        print(gcf, save_ps_path, '-dpsc', '-append', '-fillpage');
        close(gcf);
    end
end