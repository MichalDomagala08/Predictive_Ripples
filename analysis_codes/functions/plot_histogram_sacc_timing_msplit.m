function plot_histogram_sacc_timing_msplit(mspl_sacc,saveFolder,currentDependent,currentDistFeature,titl,params)
    %%% Plots histogrms of a given feature (Normaly Saccadic Distance)
    % Based on median split of some other feature!

    
    msl_upper = table2array(mspl_sacc(mspl_sacc.(currentDependent+"_mspl") == 1,currentDistFeature));
    msl_lower= table2array(mspl_sacc(mspl_sacc.(currentDependent+"_mspl") == 0,currentDistFeature));

    ldata = {msl_lower,msl_upper}; titls = {sprintf("Lower %s",currentDependent),sprintf("Upper %s",currentDependent)};
    figure("Visible","On","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);

    for dat = 1:2
        subplot(1,2,dat);  binWidth = 2/params.nbins;
        edges = -params.bounds:binWidth:params.bounds;  ;    
        histogram(ldata{dat}(ldata{dat} >= -params.bounds & ldata{dat} <= params.bounds), 'BinEdges', edges, 'FaceColor', [0.2 0.6 1], 'EdgeColor','none');
        xlim([-params.bounds params.bounds]); ylim([0 250])
        xlabel('Timing difference (s)'); ylabel('Count'); title(titls{dat})
    end
    sgtitle(sprintf("Timing Difference with median split: %s",currentDependent))
    print(gcf,fullfile(saveFolder,"Results",sprintf("%s.ps",titl)),'-dpsc','-append','-fillpage')

end