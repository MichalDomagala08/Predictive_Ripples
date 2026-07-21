% Sygnał wejściowy: x (7001 sampli)
% Fs = 500 Hz;


fsample = 500;
% 1. Usuwamy offset
data_concat = cell2mat(raw_data.trial);
x = data_concat(6,:);
xBP = ft_preproc_bandpassfilter(x,fsample, [0.3 150], ...
    [], 'fir', 'twopass');

x_clean = xBP - mean(xBP);

% 2. Bierzemy moduł, aby znaleźć peaki w obu kierunkach
x_abs = abs(zscore(x_clean));

% 3. Parametry detekcji
% - MinPeakHeight: Ustawiamy na ok. 4-5 sigma tła lub sztywny próg (np. 0.00005)
% - MinPeakDistance: 500 ms (Fs * 0.5), aby nie wykrywać tego samego IED wielokrotnie
% - MinPeakProminence: Kluczowy parametr - jak bardzo peak odstaje od otoczenia
%%
[pks, locs_pos] = findpeaks(zscore(x_clean), ...
    'MinPeakHeight', 3);  % Musi wystawać o tyle ponad lokalne tło

[thr, locs_neg] = findpeaks(-zscore(x_clean), ...
    'MinPeakHeight', 3, ...    % Wysokość progu (na podstawie danych z pliku)
    'MinPeakDistance', 300, ...   % Odstęp min 0.5s między wyładowaniami
    'MinPeakProminence', 2.5);  % Musi wystawać o tyle ponad lokalne tło

dist_matrix = abs(locs_pos(:) - locs_neg(:)');
amp_matrix = pks(:) + thr(:)';
set(0, 'DefaultFigureVisible', 'on');

for  i = 1:7001:length(x)-7001
    figure(1)
    subplot(2,1,1)
    plot(i:i+7001,xBP(i:i+7001))
    if ~isempty(locs_pos(locs_pos > i & locs_pos < i+7001))
        xline(locs_pos(locs_pos > i & locs_pos < i+7001))
    end
    if ~isempty(locs_neg(locs_neg > i & locs_neg < i+7001))
        xline(locs_neg(locs_neg > i & locs_neg < i+7001))
    end
    subplot(2,1,2)
    plot(i:i+7001,x_abs(i:i+7001))
    if ~isempty(locs_pos(locs_pos > i & locs_pos < i+7001))
        xline(locs_pos(locs_pos > i & locs_pos < i+7001),'b')
    end
    if ~isempty(locs_neg(locs_neg > i & locs_neg < i+7001))
        xline(locs_neg(locs_neg > i & locs_neg < i+7001),'r')
    end

    w = waitforbuttonpress();
end

%%