function [powerSWR,timecourseSWR,FFT_All] = getSubjRipples_review(params,hipothetical_ripples_cluster,data_viz,data_norm,channel,varargin)
    % Used for computing ripple characterstichs: Timecourse, Power TF representation and Power spectrum
    
    marg = [params.vizbefT,params.vizaftT];
    timecourseSWR = zeros(length(hipothetical_ripples_cluster),1+marg(1)+marg(2));
    powerSWR = zeros(length(hipothetical_ripples_cluster),40,1+marg(1)+marg(2));

    %%% FFT precomputation
    FreqResolution = (500)/ceil((params.FFT_Length*2+1));
    hzLimit = ceil(params.fftmaxFreq/FreqResolution);
    FrequenciesX = 501/(params.FFT_Length*2+1)*(0:ceil((params.FFT_Length*2+1)/2)-1); %Original Frequencies that we have (One sided) 
    [~,lzLimit] = min(abs(FrequenciesX-params.fftminFreq));

    FFT_All = zeros(length(hipothetical_ripples_cluster),length(FrequenciesX));
    for clust = 1:length(hipothetical_ripples_cluster)
        if hipothetical_ripples_cluster{clust}(end) > length(data_viz.trial{1}(channel,:))
            continue;
        end
        %%% Identify Maximum Peak 
        [~,I] = max(data_viz.trial{1}(channel,hipothetical_ripples_cluster{clust}));
        peak = hipothetical_ripples_cluster{clust}(1) +I-1;

        
        

        %%% take 1.5s from each side (750 timesamples around) - Staresina et al., 2016
        if peak > marg(1)*4 && peak+marg(2)*4 < length(data_viz.trial{1}(channel,:))

            %%% FFT of a Ripple Signal, 250ms not 2000 it doesnt make sense to
            %%% take more than this 
            Y = fft(data_norm.trial{1}(channel,peak-params.FFT_Length:peak+params.FFT_Length));
            FFT_All(clust,:)  = Y(1:length(FrequenciesX))/length(FrequenciesX);


            currTrial = data_norm.trial{1}(channel,peak-marg(1):peak+marg(2));

             %%% Get TimeCourse 
            timecourseSWR(clust,:) = currTrial;
    
            %%% Get Power - it is 1.5 larger to exclude Edge Artefacts 
            ft_warning("off")
            % Lo w Width may cause Issues!!! - if your Signal is Jagged
            % then increasing the Width may Help 
            pov = abs(squeeze(ft_specest_wavelet(data_norm.trial{1}(channel,peak-marg(1)*4:peak+marg(2)*4),[0:(marg(1)*4+marg(2)*4)]/500,'width',7,'freqoi',5:5:200,'verbose',0))).^2;
            
            %%% dB Transform: 
            pov_db = 10*log10(pov(:,marg(1)*3:marg(1)*3+marg(1)+marg(2))+eps); % Sometimes to help with near zero-values you should add eps 
            
            %pov_db_cor = pov_db - nanmean(pov_db(:,150:250),2);
            for freqs = 1:size(pov_db,1) % we Are Normalising it Frequency Based... but maybe we should do it differently:
                %pov(freqs,:) = (pov(freqs,:) - nanmean(pov(freqs,:)))/nanstd(pov(freqs,:));

                pov_db(freqs,:) = (pov_db(freqs,:) - nanmean(pov_db(freqs,:)))/nanstd(pov_db(freqs,:));
            end
            ft_warning("on")
            powerSWR(clust,:,:) = pov_db; %pov_db_cor(:,marg(1)*2:marg(1)*2+marg(1)+marg(2));
        end

       


    end

    path1temp      = char(varargin{3});

    if params.vizualization ==  2 || (params.vizualization ==3 && strcmp(path1temp(end-7:end-3),'Undet')) %% Check whether we visualise undetected
    tic

        goodChannels    = varargin{1};
        current_ds_name = varargin{2};
        path1      = varargin{3};
        fprintf("Processing Ripple n.: %s")

        for clust = 1:length(hipothetical_ripples_cluster)
          %%% Long Haul Visualization of EVERY RIPPlE EVENT (Why? XD) TO DO: 
            set(0, 'DefaultFigureVisible', 'off')
    
            figure('units','normalized','outerposition',[0 0 1 1],'paperOrientation','landscape');
            sgtitle(sprintf("Ripple n. %s : %s // %s ",string(clust),current_ds_name,goodChannels{channel}))
            subplot(1,3,1)
            annotation('textbox',[0.01 0.1 0.9 0.9],'FontSize',20,'String',"B",'EdgeColor','none')
            hold on
            plot(squeeze(timecourseSWR(clust,:)))
            xticks([1:(params.vizaftT+params.vizbefT)/10:params.vizaftT+params.vizbefT+1])
            xticklabels([-params.vizbefT*0.002:(params.vizaftT+params.vizbefT)*0.001/5:params.vizaftT*0.002])
            ylabel("Voltage (\muV)")
            xlabel("Time from ripple peak (s)")
            hold off
            title("Time Course")
            
            subplot(1,3,2)
            imagesc(squeeze( powerSWR(clust,:,:)))
            c  = colorbar;
            c.Label.String = 'z-power';
            axis xy
            yticks([1:4:40])
            yticklabels([1:20:200]-1)
            xticks([1:(params.vizaftT+params.vizbefT)/10:params.vizaftT+params.vizbefT+1])
            xticklabels([-params.vizbefT*0.002:(params.vizaftT+params.vizbefT)*0.001/5:params.vizaftT*0.002])
            xlabel("Time (ms)")
            ylabel("Frequency (Hz)")
            title("Power")

            subplot(1,3,3)

            figure('units','normalized','outerposition',[0 0 1 1],'paperOrientation','landscape');
            plot(FrequenciesX(lzLimit:hzLimit),abs(FFT_All(clust,lzLimit:hzLimit)).^2)
            ylabel("Power (\muV2)")
            xlabel("Frequency (Hz)")
            title("FFT Power Spectrum")

            fig_handle = gcf;
            print(fig_handle, '-dpsc','-append','-fillpage' ,  path1)

            fprintf('\b\b\b\b%4d', clust);  % \b is used to remove characters, so the number updates

            close(fig_handle)
        end
    toc
    end
    
end