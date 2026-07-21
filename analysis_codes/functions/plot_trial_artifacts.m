function plot_trial_artifacts(raw_data,iqr_byChan,range_byChan,ieds_byChan,psFile,chan)

    %%% Plot Trial by Trial timecourse for a given channel, with coloring of different type of detected artifacts
    for i = 1:numel(raw_data.trial)
        figure("Visible","off"); clf; hold on
        plot(raw_data.time{i}, raw_data.trial{i}(chan,:), 'k')

        %%% IED artefact Coloring
        art = sort(ieds_byChan{chan}{i}(:)');
        if ~isempty(art)
            yl = ylim; gaps = find(diff(art)>1); s = art([1,gaps+1]); e = art([gaps,end]);
            for k = 1:numel(s)
                patch(raw_data.time{i}([s(k) e(k) e(k) s(k)]),[yl(1) yl(1) yl(2) yl(2)],'r','FaceAlpha',.3,'EdgeColor','none')
            end
        end

        %%% IQR artefact Coloring
        art = sort(iqr_byChan{chan}{i}(:)');
        if ~isempty(art)
            yl = ylim; gaps = find(diff(art)>1); s = art([1,gaps+1]); e = art([gaps,end]);
            for k = 1:numel(s)
                patch(raw_data.time{i}([s(k) e(k) e(k) s(k)]),[yl(1) yl(1) yl(2) yl(2)],'g','FaceAlpha',.3,'EdgeColor','none')
            end
        end

        %%% Range artefact Coloring
        art = sort(range_byChan{chan}{i}(:)');
        if ~isempty(art)
            yl = ylim; gaps = find(diff(art)>1); s = art([1,gaps+1]); e = art([gaps,end]);
            for k = 1:numel(s)
                patch(raw_data.time{i}([s(k) e(k) e(k) s(k)]),[yl(1) yl(1) yl(2) yl(2)],'b','FaceAlpha',.3,'EdgeColor','none')
            end
        end
        plot(raw_data.time{i}, raw_data.trial{i}(chan,:), 'k');

    end

    title(sprintf('%s | Trial %d',raw_data.label{chan},i))
    if i==1; print(gcf,psFile,'-dpsc'); else; print(gcf,psFile,'-dpsc','-append'); end
    close(gcf);
end