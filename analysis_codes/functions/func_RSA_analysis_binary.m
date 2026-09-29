function func_RSA_analysis_binary(T_data,overlap,gr_names,dependentName,predictorName,pars,pars_restricted,an_title)


% TO DO!!!
% This function s performng RSA statstcal analyss on bnary data (i.e is before)
pars = {'rip_fix_distance_edge',"is_after","subject"};

curr_est = []; curr_tStat = []; curr_p = []; nam = []; pred = [];anname = [];
for i = 1:length(gr_names)

    %%% For all Before vs After 
    lme = fitlme(T_data((T_data.is_before | T_data.is_after) & overlap,:),sprintf("%s ~ is_after +  rip_fix_distance_edge +(1|subject)", string(gr_names{i})));
    curr_est = [curr_est;lme.Coefficients.Estimate(2:3)]; curr_tStat = [curr_tStat;lme.Coefficients.tStat(2:3)]; curr_p = [curr_p;lme.Coefficients.pValue(2:3)];
    pred = [pred; string(lme.Coefficients.Name(2:3))]; nam = [nam, string(gr_names{i}),string(gr_names{i})]; anname = [anname; "RipProx_Restr"; "RipProx_Restr"];


    figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
    sgtitle(sprintf("All: %s: P: %2.3f T: %2.3f",string(strrep(gr_names{i},'_','-')),lme.Coefficients.pValue(2),lme.Coefficients.tStat(2)))
    subplot(1,2,1)
    boxplot(T_data((T_data.is_before | T_data.is_after)  & overlap,:).(string(gr_names{i})), T_data((T_data.is_before | T_data.is_after)  & overlap,:).is_after, 'Labels', {'after','before'}); % upewnij się kolejność etykiet
    ylabel(string(strrep(gr_names{i},'_','-'))); title(sprintf("Restricted Before vs After %s",string(strrep(gr_names{i},'_','-'))));
    title("After - Before Boxplot");

    subplot(1,2,2)
    plot_linear_model_bin2(T_data( T_data.is_after | T_data.is_before  & overlap,:),lme, string(gr_names{i}),'is_after',{'rip_fix_distance_edge',"is_after","subject"},1,1,250)
    xticks ([sort([xticks 25])])
    xticklabels(string(str2double(xticklabels) * 2));
    xlabel("Absolute Saccade to Ripple distance (ms)");    
    title("Closeness linear Model: ");

    print(gcf,fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_aft_vs_bef.ps',params.exact_analysisName)),'-dpsc','-append','-fillpage');
    close(gcf);

    % For restricted (-/+ params.max ms around Onset)
    lme = fitlme(T_data((T_data.is_before_restricted | T_data.is_after_restricted )  & overlap ,:),sprintf("%s ~ is_after_restricted +  rip_fix_distance_edge +(1|subject) ", string(gr_names{i})));
    curr_est = [curr_est;lme.Coefficients.Estimate(2:3)]; curr_tStat = [curr_tStat;lme.Coefficients.tStat(2:3)]; curr_p = [curr_p;lme.Coefficients.pValue(2:3)];
    pred = [pred; string(lme.Coefficients.Name(2:3))]; nam = [nam, string(gr_names{i}),string(gr_names{i})]; anname = [anname; "RipProx_Restr"; "RipProx_Restr"];

    figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
    sgtitle(sprintf("Restricted (+/- %d ms): %s: P: %2.3f T: %2.3f",params.max_time/2,string(strrep(gr_names{i},'_','-')),lme.Coefficients.pValue(2),lme.Coefficients.tStat(2)))
    subplot(1,2,1)
    boxplot(T_data((T_data.is_before_restricted | T_data.is_after_restricted )  & overlap ,:).(string(gr_names{i})), T_data((T_data.is_before_restricted | T_data.is_after_restricted)  & overlap ,:).is_after_restricted, 'Labels', {'after','before'}); % upewnij się kolejność etykiet
    ylabel(string(strrep(gr_names{i},'_','-'))); title(sprintf("Restricted Before vs After %s",string(strrep(gr_names{i},'_','-'))));
    title("After - Before Boxplot");

    subplot(1,2,2)
    plot_linear_model_bin2(T_data((T_data.is_after_restricted | T_data.is_before_restricted)  & overlap,:),lme, string(gr_names{i}),'is_after_restricted',{'rip_fix_distance_edge',"is_after_restricted","subject"},1,1,20)
    xticks ([sort([xticks 25])])
    xticklabels(string(str2double(xticklabels) * 2));
    xlabel("Absolute Saccade to Ripple distance (ms)");    
    title("Closeness linear Model: ");

    print(gcf,fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_aft_vs_bef.ps',params.exact_analysisName)),'-dpsc','-append','-fillpage');
    close(gcf);

end


psDir =fullfile(saveFolder,params.exact_analysisName,sprintf('%s_ripples_fix_aft_vs_bef.ps',params.exact_analysisName));
system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 


T = table(anname,nam',pred,curr_est,curr_tStat,curr_p, ...
'VariableNames', {'anName','Dependent','Predictor','Coef','tStat','pValue'});


fn = fullfile(saveFolder,params.exact_analysisName,char(sprintf('LME_%s_Restr_Fixations_%s.xlsx',params.exact_analysisName,string(params.max_time)))); sheet = matlab.lang.makeValidName('saccades_aft_vs_bef');
writetable(T,fn,'Sheet',sheet,'WriteVariableNames',true);
ex = actxserver('Excel.Application'); wb = ex.Workbooks.Open(fn); sht = wb.Sheets.Item(sheet);
sigRows = find(T.pValue < 0.05); for r = sigRows'; sht.Range(sprintf('A%d:G%d',r+1,r+1)).Font.Bold = true; end
wb.Save(); wb.Close(false); ex.Quit(); delete(ex);



end