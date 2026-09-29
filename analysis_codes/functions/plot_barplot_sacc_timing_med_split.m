function plot_barplot_sacc_timing_med_split(mspl_sacc,saveFolder,currentDependent,titl,params)

    figure("Visible","On","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
    ripsac_n1 = {'rippleNumBefore_saccade','rippleNumIn_saccade','rippleNumAfter_saccade','rippleNumAfter2_saccade','rippleNumAfter3_saccade','rippleNumAfter4_saccade'};
    ripsac_n2 = {'rippleNumBefore_saccade_M','rippleNumIn_saccade_M','rippleNumAfter_saccade_M','rippleNumAfter2_saccade_M','rippleNumAfter3_saccade_M','rippleNumAfter4_saccade_M'};

    subplot(1,2,1);
    bar([sum(table2array(mspl_sacc(mspl_sacc.(currentDependent+"_mspl") == 0,ripsac_n2)));...
         sum(table2array(mspl_sacc(mspl_sacc.(currentDependent+"_mspl") == 1,ripsac_n2)));]')
    legend("Lower","Upper"); title("In Saccade Time")
    xticks([1,2,3,4,5,6]);xticklabels(["before",'Sac','Fix1','Fix2','Fix3','Fix4'])


    subplot(1,2,2);
    bar([sum(table2array(mspl_sacc(mspl_sacc.(currentDependent+"_mspl") == 0,ripsac_n1)));...
         sum(table2array(mspl_sacc(mspl_sacc.(currentDependent+"_mspl") == 1,ripsac_n1)));]')
    legend("Lower","Upper"); title("In Mean Saccade Time")
    xticks([1,2,3,4,5,6]);xticklabels(["before",'Sac','Fix1','Fix2','Fix3','Fix4'])

    sgtitle(sprintf("BarTimings with median split: %s",currentDependent))
    print(gcf,fullfile(saveFolder,"Results",sprintf("%s.ps",titl)),'-dpsc','-append','-fillpage')
end