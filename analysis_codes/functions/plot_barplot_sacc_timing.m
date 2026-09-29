function plot_barplot_sacc_timing_med_split(mspl_sacc,saveFolder,currentDependent)

    figure("Visible","On","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);

    subplot(1,2,1);
    bar([sum(table2array(mspl_sacc(mspl_sacc.(currentDependent+"_mspl") == 0,64:69)));...
         sum(table2array(mspl_sacc(mspl_sacc.(currentDependent+"_mspl") == 1,64:69)));]')
    legend("Lower","Upper"); title("Ripple numbers in saccade Distance")
    xticks([1,2,3,4,5,6]);xticklabels(["before",'Sac','Fix1','Fix2','Fix3','Fix4'])


    subplot(1,2,2);
    bar([sum(table2array(mspl_sacc(mspl_sacc.(currentDependent+"_mspl") == 0,[58,59,60,61,62,63])));...
         sum(table2array(mspl_sacc(mspl_sacc.(currentDependent+"_mspl") == 1,[58,59,60,61,62,63])));]')
    legend("Lower","Upper"); title("Ripple numbers in mean saccade Distance")
    xticks([1,2,3,4,5,6]);xticklabels(["before",'Sac','Fix1','Fix2','Fix3','Fix4'])

    sgtitle(sprintf("Timing Difference with median split: %s",currentDependent))
    print(gcf,fullfile(saveFolder,"Results","Lower_Upper_Ripple_Sac_timing.ps"),'-dpsc','-append','-fillpage')
