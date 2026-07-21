function  plot_ripple_contrasts(data,isOn,pathName)
    % Plots Line and violin plots for our main contrasts of data: 
    % - data - sructure containing field responsible for condition and
    %    measure. last field is table of our data in long form
    % - isOn - if plots are to be visible
    % - pathName - savepath

    % Assuming that LongForm is last 
    set(0, 'DefaultFigureVisible', isOn); 
    fieldnam = fieldnames(data);
    longTable = data.LongForm;
    disp(fieldnam)
    for name = 1:length({fieldnam{1:end-1}})
        figure();  
        sgtitle(fieldnam{name})

        %%% Online/Offline difference - Line Plot
        subplot(1,2,1)
        hold on
        plot(ones(size(data.(fieldnam{name})(:,1))),data.(fieldnam{name})(:,1),'.')
        plot(ones(size(data.(fieldnam{name})(:,2)))+1,data.(fieldnam{name})(:,2),'.')
        plot(1,[mean(data.(fieldnam{name})(:,1),"omitnan")],'.',"Color",[0 0.4470 0.7410],"MarkerSize",15)
        plot(2,[mean(data.(fieldnam{name})(:,2),"omitnan")],'.',"Color",[0.8500 0.3250 0.0980],"MarkerSize",15)
        plot([1,2],[mean(data.(fieldnam{name})(:,1),"omitnan"), mean(data.(fieldnam{name})(:,2),"omitnan")],'o',"Color",[0.9290 0.6940 0.1250])

        for i = 1:size(data.(fieldnam{name}),1)
            if data.(fieldnam{name})(i,1) > data.(fieldnam{name})(i,2)
                plot([1,2],[data.(fieldnam{name})(i,1) data.(fieldnam{name})(i,2)],"Color",[0 0.4470 0.7410,0.3])
            else
                plot([1,2],[data.(fieldnam{name})(i,1) data.(fieldnam{name})(i,2)],"Color",[0.8500 0.3250 0.0980,0.3])
            end
        end
        plot([1,2],[mean(data.(fieldnam{name})(:,1),"omitnan") mean(data.(fieldnam{name})(:,2),"omitnan")],"Color",[0.9290 0.6940 0.1250],"LineWidth",1.5)
        xlim([0 3])
        xticks([1 2])
        xticklabels({'Offline'  'Online'})


        %%% Online/Offline difference violin plot
        subplot(1,2,2)
        datVil = table2array(longTable(:,3+name));
        % Offline
        vs = Violin({datVil(strcmp(cellstr(table2cell(longTable(:,3))),'Off'))},1,...
            'QuartileStyle','none',... % boxplot, none
            'DataStyle', 'none',... % scatter, histogram
            'ShowMean',false,...
            "HalfViolin","left","MedianColor",[0.9290 0.6940 0.1250]);
        vs = Violin({datVil(strcmp(cellstr(table2cell(longTable(:,3))),'Off'))},1,...
            'QuartileStyle','none',... % boxplot, none
            'DataStyle', 'histogram',... % scatter, histogram
            'ShowMean',false,...
            "HalfViolin","left","MedianColor",[0.9290 0.6940 0.1250]);
        %Online
        vs = Violin({datVil(strcmp(cellstr(table2cell(longTable(:,3))),'On'))},2,...
            'QuartileStyle','none',... % boxplot, none
            'DataStyle', 'none',... % scatter, histogram
            'ShowMean',false,...
            "HalfViolin","right","MedianColor",[0.9290 0.6940 0.1250]);
        vs = Violin({datVil(strcmp(cellstr(table2cell(longTable(:,3))),'On'))},2,...
            'QuartileStyle','none',... % boxplot, none
            'DataStyle', 'histogram',... % scatter, histogram
            'ShowMean',false,...
            "HalfViolin","right","MedianColor",[0.9290 0.6940 0.1250]);
        xticks([1,2]);
        xticklabels(["Offline","Online"]);
        
        
        fig_handle = gcf;
        saveas(fig_handle,strrep(strrep(pathName,"XXX",fieldnam{name}),"png",'svg')) % CHANGE TO SVG IF NEEDED!

    end

end