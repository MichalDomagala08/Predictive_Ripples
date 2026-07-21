function [powerSWR,timecourseSWR,ripplerangeSWR,FFT_All] = func_RippleChar(params,hipothetical_ripples_cluster,data_viz,data_norm,channel,varargin)
    % Used for computing ripple characterstichs: Timecourse, Power TF representation and Power spectrum

    marg   = [params.vizbefT, params.vizaftT];
    sigLen = length(data_viz(channel,:));

    timecourseSWR = nan(length(hipothetical_ripples_cluster), 1+marg(1)+marg(2));
    ripplerangeSWR = nan(length(hipothetical_ripples_cluster), 1+marg(1)+marg(2));

    powerSWR      = nan(length(hipothetical_ripples_cluster), 40, 1+marg(1)+marg(2));

    %%% FFT precomputation
    FreqResolution = (500)/ceil((params.FFT_Length*2+1));
    hzLimit        = ceil(params.fftmaxFreq/FreqResolution);
    FrequenciesX   = 501/(params.FFT_Length*2+1)*(0:ceil((params.FFT_Length*2+1)/2)-1);
    [~,lzLimit]    = min(abs(FrequenciesX-params.fftminFreq));
    fftLen         = 2*params.FFT_Length + 1;
    FFT_All        = nan(length(hipothetical_ripples_cluster), length(FrequenciesX));

    for clust = 1:length(hipothetical_ripples_cluster)
        if hipothetical_ripples_cluster{clust}(end) > sigLen, continue; end

        %%% Identify Maximum Peak
        [~,I] = max(data_viz(channel, hipothetical_ripples_cluster{clust}));
        peak  = hipothetical_ripples_cluster{clust}(1) + I - 1;

        %%% FFT — zero-padding do stałej długości
        f_lo = max(1, peak-params.FFT_Length);
        f_hi = min(sigLen, peak+params.FFT_Length);
        sig_pad = zeros(1, fftLen);
        sig_pad(f_lo-peak+params.FFT_Length+1 : f_hi-peak+params.FFT_Length+1) = data_norm(channel, f_lo:f_hi);
        Y = fft(sig_pad);
        FFT_All(clust,:) = Y(1:length(FrequenciesX)) / fftLen;

        %%% TimeCourse — NaN padding, peak zawsze w centrum
        p_lo = max(1, peak-marg(1));
        p_hi = min(sigLen, peak+marg(2));
        tc   = nan(1, 1+marg(1)+marg(2)); tc_rr = tc;
        tc(p_lo-peak+marg(1)+1 : p_hi-peak+marg(1)+1) = data_norm(channel, p_lo:p_hi);
        tc_rr(p_lo-peak+marg(1)+1 : p_hi-peak+marg(1)+1) = data_viz(channel, p_lo:p_hi);
        timecourseSWR(clust,:) = tc;
        ripplerangeSWR(clust,:) = tc_rr;

        %%% Power — wavelet na dostępnym sygnale, wytnij środek
        t_lo = max(1, peak-marg(1)*4);
        t_hi = min(sigLen, peak+marg(2)*4);

        pov    = abs(squeeze(ft_specest_wavelet(data_norm(channel, t_lo:t_hi), ...
                    [0:(t_hi-t_lo)]/500, 'width',7, 'freqoi',5:5:200, 'verbose',0))).^2;
        
        % indeks środka (peak) w oknie wavelet
        peak_in_tw = peak - t_lo + 1;
        pov_lo     = max(1, peak_in_tw - marg(1));
        pov_hi     = min(size(pov,2), peak_in_tw + marg(2));
        pov_win    = nan(size(pov,1), 1+marg(1)+marg(2));
        pov_win(:, pov_lo-peak_in_tw+marg(1)+1 : pov_hi-peak_in_tw+marg(1)+1) = pov(:, pov_lo:pov_hi);
        pov_db = 10*log10(pov_win + eps);
        for freqs = 1:size(pov_db,1)
            pov_db(freqs,:) = (pov_db(freqs,:) - nanmean(pov_db(freqs,:))) / nanstd(pov_db(freqs,:));
        end
        powerSWR(clust,:,:) = pov_db;
    end
end