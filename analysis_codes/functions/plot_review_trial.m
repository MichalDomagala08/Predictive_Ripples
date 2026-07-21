
function plot_review_trial(data_artf)
    i = 1;
    while i >= 1 && i <= 80
        figure(1); clf;
        plot(data_artf.time{i}, data_artf.trial{i})
        legend(data_artf.label);
        title(sprintf('Trial %d / 80', i));
    
        [~, ~, key] = ginput(1);  % czeka na kliknięcie/klawisz
    
        % Odczytaj ostatni klawisz
        key = get(gcf, 'CurrentKey');
        if strcmp(key, 'rightarrow')
            i = min(i + 1, 80);
        elseif strcmp(key, 'leftarrow')
            i = max(i - 1, 1);
        elseif strcmp(key, 'escape')
            break;
        end
    end
end