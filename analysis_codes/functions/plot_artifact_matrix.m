function plot_artifact_matrix(data,saveFolder,currSubjName)

    figure('Visible','Off','Units','normalized','Position',[0 0 0.7 1]);
    nChan = length(data.label); nCols = ceil(sqrt(nChan)); nRows = ceil(nChan / nCols);

    tl =tiledlayout(nRows, nCols, 'TileSpacing','compact','Padding','compact');
    for chanNum = 1:length(data.label)
        tmp      = cellfun(@(x) x(chanNum,:), data.trial, 'UniformOutput', false);
        trialMat = vertcat(tmp{:});  % [nTrials x nSamples]
    
        nexttile
        imagesc(trialMat); colorbar();
        ylabel(data.label{chanNum});
    end
    
    xlabel('Sample'); title(tl, 'Trials x Time per Channel');
    exportgraphics(gcf, fullfile(saveFolder, 'ArtifactRejection', currSubjName, sprintf('%s_artifacts_matrix_interp.png', currSubjName)), 'Resolution', 300);
end