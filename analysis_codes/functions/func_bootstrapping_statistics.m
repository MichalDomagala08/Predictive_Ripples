%%% This script computes bootstrapping statistics for a given 

function func_bootstrapping_statistics(saccades_table_cleaned,subj,saccadic_features,contrast,numIter,saveFolder)

    
    rippleSacRate = groupcounts(saccades_table_cleaned(logical(saccades_table_cleaned.(contrast)),:),"subject_new");
    rowNumsBootstrap = nan(numIter,sum(saccades_table_cleaned.(contrast)~=0));
    
    %%% Compute bootsrrapping permutations
    parfor b = 1:numIter
        current_opposites = [];
        for sb  = 1:length(subj)
        
            %Get current candidate rows 
            nDraw   = rippleSacRate(rippleSacRate.subject_new == subj{sb},:).GroupCount;
            if isempty(nDraw); continue; end;
            candIdx = find(strcmp(saccades_table_cleaned.subject_new,subj{sb}) & ~saccades_table_cleaned.(contrast));
            if isempty(candIdx), sampledIdx = []; else;  selLocal = randsample(numel(candIdx),nDraw ,numel(candIdx)<nDraw); sampledIdx = candIdx(selLocal); end
        
            current_opposites = [current_opposites; sampledIdx];
        end

        if length(current_opposites) ~= sum(saccades_table_cleaned.(contrast)~=0); 
            error("The length of both is not ok!"); end
        rowNumsBootstrap(b,:) = current_opposites;
    end
    
    
    
    %%% Compute statistic of saccadic features
    for sac = saccadic_features
    
        feature_name = saccades_table_cleaned.Properties.VariableNames{sac};
    
        true_value = table2array(saccades_table_cleaned(logical(saccades_table_cleaned.(contrast)),sac))';
        bootrstap_values = reshape(table2array(saccades_table_cleaned(rowNumsBootstrap,sac)),size(rowNumsBootstrap));
        stat_true = mean(true_value,"omitnan");  stat_boot = mean(bootrstap_values,2,"omitnan"); 
        center = mean(stat_boot,"omitnan"); p_two_center = (1 + sum(abs(stat_boot - center) >= abs(stat_true - center))) / (1 + numel(stat_boot));
    
    
        figure("Visible","Off","PaperOrientation","landscape","Units","normalized","Position",[0 0 1 1]);
        b = histogram(stat_boot); hold on; xline(stat_true,'r','LineWidth',2);
        center = median(stat_boot);
        text(stat_true, 0.95*max(b.Values), sprintf('True=%.3f\n', stat_true), 'HorizontalAlignment','left', 'BackgroundColor','w');
        title(sprintf("%s P: %.3f, stat: %.3f",feature_name,p_two_center,(stat_true - center)/center))
        xlabel(feature_name); ylabel("Frequency")
    
        print(gcf,fullfile(saveFolder,sprintf("Bootstrap_%s.ps",contrast)),'-dpsc','-append','-fillpage')
        close(gcf)
    end
    %%
    psDir = fullfile(saveFolder,sprintf("Bootstrap_%s.ps",contrast));
    system(strcat("gswin64c ", "-sDEVICE=pdfwrite ","-o ", strrep(psDir,".ps", ".pdf "),psDir));    delete(psDir); 
          
end


