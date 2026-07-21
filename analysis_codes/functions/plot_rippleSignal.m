
function plot_rippleSignal(timecourseData,powerData,fftData,Ripple_Distance,RippleQuantity,path,params,titl,varargin)
    % Plotting averager Riple (or IED) characteristics:
    % dataRipple - struct containing info about 


    if ~isempty(varargin)
        indivPath = varargin{1};
    end

    %------------------------%
    %%% Timecourse Average Plot %%%
    %------------------------%

 
    set(0, 'DefaultFigureVisible', 'off')

    figure('paperOrientation','landscape');
    annotation('textbox',[0.01 0.1 0.9 0.9],'FontSize',20,'String',"B",'EdgeColor','none');
    hold on
    if size(timecourseData,1) >1
        stdshade(squeeze(timecourseData),0.5,[0.8500 0.3250 0.0980]);
    else
        plot(timecourseData);
    end
    title(sprintf("%s: Average of n = %s ripples",string(titl),string(RippleQuantity)));
    xticks([1:(params.vizaftT+params.vizbefT)/10:params.vizaftT+params.vizbefT+1]);
    xticklabels([-params.vizbefT*0.002:(params.vizaftT+params.vizbefT)*0.001/5:params.vizaftT*0.002]);
    ylabel("Voltage (\muV)");
    xlabel("Time from ripple peak (s)");
    hold off
    fig_handle = gcf;
    if contains(path,'.ps')
        print(fig_handle, '-dpsc','-append','-fillpage' , strrep(path,'XXX',"Ripple"));
    else
        saveas(fig_handle,  strrep(strrep(path,'XXX',"Ripple"),'ps','svg'));
    end
    if ~isempty(varargin); saveas(fig_handle, strrep(indivPath,'XXX',"Ripple")); end

    %------------------------%
    %%% Power Average Plot %%%
    %------------------------%
    figure('paperOrientation','landscape');
    annotation('textbox',[0.01 0.1 0.9 0.9],'FontSize',20,'String',"C",'EdgeColor','none')
    
    imagesc(squeeze(powerData))
    
    c  = colorbar;
    c.Label.String = 'dB';
    axis xy
    yticks([1:4:40]);
    yticklabels([1:20:200]-1);
    xticks([1:(params.vizaftT+params.vizbefT)/10:params.vizaftT+params.vizbefT+1]);
    xticklabels([-params.vizbefT*0.002:(params.vizaftT+params.vizbefT)*0.001/5:params.vizaftT*0.002]);
    xlabel("Time (ms)");
    ylabel("Frequency (Hz)");
    title(sprintf("%s",string(titl)));

    fig_handle = gcf;
    try
    if contains(path,'.ps') 
        print(fig_handle, '-dpsc','-append','-fillpage' ,  strrep(path,'XXX',"Power"));
    else
        saveas(fig_handle,  strrep(strrep(path,'XXX',"Power"),'ps','svg'));  end

    catch
        error(sprintf("There was a problem  with generating: %s",titl))
    end
    if ~isempty(varargin); saveas(fig_handle, strrep(indivPath,'XXX',"Power")); end
    close(fig_handle)
    %----------------------------%
    %%% Inter-ripple Histogram %%%
    %----------------------------%
    if RippleQuantity > 5

        figure('paperOrientation','landscape');
        annotation('textbox',[0.01 0.1 0.9 0.9],'FontSize',20,'String',"D",'EdgeColor','none')
        histogram(Ripple_Distance/500,50,"FaceColor","#D95319")
        xlabel("Inter-ripple interval (ms)");
        ylabel("Probability ")
        title(sprintf("%s:, mean = %0.f ms // iqr = %0.f",string(titl),mean(Ripple_Distance/500),iqr(Ripple_Distance/500)))
        fig_handle = gcf;
        if contains(path,'.ps')
            print(fig_handle, '-dpsc','-append','-fillpage' ,  strrep(path,'XXX',"Hist"));
        else
            saveas(fig_handle,  strrep(strrep(path,'XXX',"Hist"),'ps','svg'));
        end
        if ~isempty(varargin); saveas(fig_handle, strrep(indivPath,'XXX',"Hist")); end
        close(fig_handle);
    end

    %---------------------%
    %%% FFT of a  Signal %%%
    %---------------------%

    if size(fftData,1) >1

        FreqResolution = (500)/ceil((params.FFT_Length*2+1));
        hzLimit = ceil(params.fftmaxFreq/FreqResolution);
        FrequenciesX = 501/(params.FFT_Length*2+1)*(0:ceil((params.FFT_Length*2+1)/2)-1); %Original Frequencies that we have (One sided) 
        [~,lzLimit] = min(abs(FrequenciesX-params.fftminFreq));
        %lzLimit = lzLimit;
    
        %%% Fooof:
        f_range = [lzLimit hzLimit];
        settings = struct('aperiodic_mode','knee', ...
                          'peak_width_limits',[2 20], ...
                          'min_peak_height',0.1, ...
                          'max_n_peaks',6);
        try
            out = fooof(FrequenciesX(lzLimit-5:hzLimit), mean(abs(fftData(:,lzLimit-5:hzLimit)).^2,1), f_range, settings);
            yfit = 10.^( out.aperiodic_params(1) - log10( out.aperiodic_params(2) + FrequenciesX.^ out.aperiodic_params(3)));
            figure('units','normalized','outerposition',[0 0 1 1],'paperOrientation','landscape'); hold on;
            stdshade(abs(fftData(:,lzLimit:hzLimit)).^2,0.5,[0.8500 0.3250 0.0980],FrequenciesX(lzLimit:hzLimit));
            plot(FrequenciesX(lzLimit:hzLimit), yfit(lzLimit:hzLimit), 'r--', 'LineWidth',2);
        
            ylabel("Power (\muV2)");
            xlabel("Frequency (Hz)");
            title(sprintf("%s: FFT Power Spectrum",string(titl)));
        catch
            figure('units','normalized','outerposition',[0 0 1 1],'paperOrientation','landscape'); hold on;
            title(sprintf("%s: FFT Power Spectrum: ERROR! %s",string(titl),string(size(fftData))));
        end
    
    
       
        fig_handle = gcf;
        if contains(path,'.ps')
            print(fig_handle, '-dpsc','-append','-fillpage' ,  strrep(path,'XXX',"FFT"));
        else
            saveas(fig_handle,  strrep(strrep(path,'XXX',sprintf("FFT_%s_%s",string(lzLimit),string(hzLimit))),'ps','svg'));
        end
        if ~isempty(varargin); saveas(fig_handle, strrep(indivPath,'XXX',"FFT")); end
        close(fig_handle);

    end
end 