function plot_single_ripple(powerData,timecourseData,ripplerangeData, currSubjName,channelLabel,params,path)
% this function plots a  single ripple event with Normal, Ripple Range and Power Estimations
     set(0, 'DefaultFigureVisible', 'off')

    figure('paperOrientation','landscape');

    sgtitle(sprintf("%s // %s",currSubjName,channelLabel))
    subplot(2,2,1)
    hold on
    plot(timecourseData);

    xticks([1:(params.vizaftT+params.vizbefT)/10:params.vizaftT+params.vizbefT+1]);
    xticklabels([-params.vizbefT*0.002:(params.vizaftT+params.vizbefT)*0.001/5:params.vizaftT*0.002]);
    ylabel("Voltage (\muV)");
    xlabel("Time from ripple peak (s)");
    hold off
    title("Raw Timecourse")


    subplot(2,2,3)
    hold on
    plot(ripplerangeData);

    xticks([1:(params.vizaftT+params.vizbefT)/10:params.vizaftT+params.vizbefT+1]);
    xticklabels([-params.vizbefT*0.002:(params.vizaftT+params.vizbefT)*0.001/5:params.vizaftT*0.002]);
    ylabel("Voltage (\muV)");
    xlabel("Time from ripple peak (s)");
    hold off
    title("Ripple Range Timecourse")

    %%% Power  Plot %%
    subplot(2,2,[2,4])

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
    title("Time Frequency Power");

    fig_handle = gcf;
    print(fig_handle, '-dpsc','-append','-fillpage' ,  path);
    close(gcf)

end