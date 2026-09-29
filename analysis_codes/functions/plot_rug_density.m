function plot_rug_density(data,titl,rug_alpha,varargin)
    %%% This function plots Kernel Density Estimation of Trial-wise phenomenon

    %Parameters

    if ~isempty(varargin)
        ks_bandwidth = varargin{1}; 
        if  length(varargin) >1; line_color = varargin{2}; else; line_color = [0 0.447 0.741]; end
        if  length(varargin) >2; rug_color = varargin{3}; else; rug_color = [0.8500 0.3250 0.0980 rug_alpha]; end
        if  length(varargin) >3; rug_height_multip = varargin{4}; else; rug_height_multip = 0.05; end
    else
        ks_bandwidth = 200;
        line_color = [0 0.447 0.741];
        rug_color = [0.8500    0.3250    0.0980 rug_alpha];
        rug_height_multip = 0.05;
    end

    if ~rug_alpha; rug_color = [0.8500    0.3250    0.0980]; end

    hold on;

    % KS Densiy Plot
    [f, xi] = ksdensity(data, 'Bandwidth', ks_bandwidth); % density with kernel estimataion
    plot(xi, f, 'LineWidth', 2, 'Color', line_color); % Niebieski

    % Rug Plot
    rug_height = max(f) * rug_height_multip; 
    line([data; data], [zeros(size(data)); ones(size(data)) * rug_height], ...
         'Color', rug_color, ... % Kolor szary z 20% alfa
         'LineWidth', 1);

    % Trial Boundaries
    xline(2001,'--'); xline(5001,'--');
    
    % Stylizacja osi
    xticks(1:500:7001)
    xticklabels([-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10])
    xlabel('Time (s)')
    ylabel('Density')
    title(titl);

    grid on; hold off;

end