function plot_ripple_trials(it, channel, hipothetical_ripples_cluster, data_viz, data_eeg, data_ripples_viz, psFile)
    figure("Visible", "Off","PaperOrientation","landscape",'Units','normalized','Position',[0 0 1 1]);
    
    subplot(3,1,1)
    plot(data_eeg.trial{it}(channel,:))
    xticks(1:500:length(data_viz.trial{it}(channel,:)))
    xticklabels([-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10])
    drawnow;
    yl = ylim;
    if isempty(hipothetical_ripples_cluster{it})
        art = [];
    else
        art = sort([hipothetical_ripples_cluster{it}{:}]);
    end
    if ~isempty(art)
        gaps = find(diff(art) > 1);
        s = art([1, gaps+1]); e = art([gaps, end]);
        for k = 1:numel(s)
            patch([s(k) e(k) e(k) s(k)], [yl(1) yl(1) yl(2) yl(2)], ...
                'g', 'FaceAlpha', .3, 'EdgeColor', 'none')
        end
    end
    plot(data_eeg.trial{it}(channel,:))

    subplot(3,1,2)
    plot(data_viz.trial{it}(channel,:))
    xticks(1:500:length(data_viz.trial{it}(channel,:)))
    xticklabels([-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10])
    drawnow;
    yl = ylim;
    if ~isempty(art)
        for k = 1:numel(s)
            patch([s(k) e(k) e(k) s(k)], [yl(1) yl(1) yl(2) yl(2)], ...
                'g', 'FaceAlpha', .3, 'EdgeColor', 'none')
        end
    end
    plot(data_viz.trial{it}(channel,:))

    subplot(3,1,3)
    hold on;
    yline(2)
    plot(data_ripples_viz.trial{it}(channel,:))
    xticks(1:500:length(data_viz.trial{it}(channel,:)))
    xticklabels([-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10])
    drawnow;
    yl = ylim;
    if ~isempty(art)
        for k = 1:numel(s)
            patch([s(k) e(k) e(k) s(k)], [yl(1) yl(1) yl(2) yl(2)], ...
                'g', 'FaceAlpha', .3, 'EdgeColor', 'none')
        end
    end
    plot(data_ripples_viz.trial{it}(channel,:))

    print(gcf,psFile,'-dpsc','-append','-bestfit');
    close(gcf)
end