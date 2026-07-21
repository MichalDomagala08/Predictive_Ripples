
function [data_ripples,data_viz] = func_rippleband_filtering(dataForPreproc,params)
    % This function is utilised for filtering in ripple range (or any specified) and computing Envelope WITH Hilbert Transform:


    % Copy the data to new structures 
    data_ripples = dataForPreproc;
    data_viz     = dataForPreproc;
    % Mirror the signal in both ways:

    for trial = 1:length(dataForPreproc.trial)

        currentTrial_orig = dataForPreproc.trial{trial};
        trial_t = zeros(size(currentTrial_orig));

        currentTrial = [currentTrial_orig(:, end:-1:1), currentTrial_orig, currentTrial_orig(:, end:-1:1)];
        swr_range_all=zeros(size(trial_t));
        for Currchannel = 1:size(dataForPreproc.label,1)
            % Using BandPass Fieldtrip twopass FIR filter for Every Channel
            swr_range = ft_preproc_bandpassfilter(currentTrial(Currchannel,:),dataForPreproc.fsample ,[params.lowpassfreq,params.highpassfreq],[],'fir','twopass');
            hilbertComp = abs(hilbert(swr_range)); % computing abs.^2 from Hilber transform

            swr_range_all(Currchannel,:) = swr_range(length(currentTrial_orig):2*length(currentTrial_orig)-1);
            trial_t(Currchannel,:) = hilbertComp(length(currentTrial_orig):2*length(currentTrial_orig)-1);
            
        end
        data_ripples.trial{trial} = trial_t;
        data_viz.trial{trial} = swr_range_all;   % Save Non-hilber data for visualization

    end
end


