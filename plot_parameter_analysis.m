function plot_parameter_analysis(data_path,s)

%this is for ploting the analysis result based on GenomeScanMutaion.m
%data_path is same as save_path in 'analyse_and_save_parameters.m', i.e.
%optimised parameters output .mat file.

%s-this can be changed based on the parameters you want to plot,but also
%all 'sopt_list' and 's' need to be changed to other parameters.


    % Get a list of all .mat files in the data path
    mat_files = dir(fullfile(data_path, '*.mat'));
    
    % Initialize arrays to store data
    x0_values = [0, 0.001, 0.01, 0.1, 0.5];
    dt_values = [10, 20, 30, 40, 50, 100, 150, 200, 300, 400, 500];
    num_x0 = length(x0_values);
    num_dt = length(dt_values);
    
    abs_diff_N_list = zeros(num_x0, num_dt);
    std_N_list = zeros(num_x0, num_dt);
    
    for file_idx = 1:length(mat_files)
        file_name = mat_files(file_idx).name;
        
        % Extract x0 and dt values from the file name
        x0 = str2double(regexp(file_name, 'x0_([\d.]+)_dt_', 'tokens', 'once'));
        dt = str2double(regexp(file_name, 'dt_(\d+)_analysis.mat', 'tokens', 'once'));
        
        % Load the data from the file
        load(fullfile(data_path, file_name), 'sopt_list');
        
        % Calculate absolute difference and standard deviation values
        abs_diff_N = median(abs(log10(sopt_list) - log10(s)));
        std_N = std(log10(sopt_list) - log10(s));
        
        % Find the corresponding indices in x0_values and dt_values
        x0_idx = find(x0_values == x0);
        dt_idx = find(dt_values == dt);
        
        abs_diff_N_list(x0_idx, dt_idx) = abs_diff_N;
        std_N_list(x0_idx, dt_idx) = std_N;
        std_err_N_list = std_N_list ./ sqrt(500);
    end
    
    % Plot the data
     figure;
    hold on;
    for x0_idx = 1:num_x0
        line_color = [96, 96, 96] / 255; % Default line color
        marker_color = [106/255, 76/255, 36/255]; % Default marker color
        
        if x0_idx == 2
            line_color = 'black'; % Change line color for x0_0.001
            marker_color = [36/255, 169/255, 225/255]; % Change marker color for x0_0.001

        elseif x0_idx == 3
            line_color = 'black'; % Change line color for x0_0.001
            marker_color = [106/255, 225/255, 36/255]; % Change marker color for x0_0.001

        elseif x0_idx == 4
            line_color = 'black'; % Change line color for x0_0.001
            marker_color = [225/255, 213/255, 36/255]; % Change marker color for x0_0.001

        elseif x0_idx == 5
            line_color = 'black'; % Change line color for x0_0.001
            marker_color = [225/255, 36/255, 76/255]; % Change marker color for x0_0.001

        end
        
       
        line_handle = errorbar(dt_values, abs_diff_N_list(x0_idx, :), std_err_N_list(x0_idx, :), 'LineStyle', '-', 'Color', line_color);
        set(line_handle, 'Marker', 'o', 'MarkerFaceColor', marker_color, 'MarkerEdgeColor', 'black', 'MarkerSize', 8);
    end
    
    xlim([0,510])
    xlabel('dt');
    ylabel('$Median\ Abs.\ Difference\ (|\log_{10}(s^*)-\log_{10}(s)|)$','Interpreter','latex');
    title('$Median\ Abs.\ Difference\ (|\log_{10}(s^*)-\log_{10}(s)|)\ Analysis\ for\ fixing\ N\ \mu$',Interpreter='latex');
    legend('x0=0', 'x0=0.001', 'x0=0.01', 'x0=0.1', 'x0=0.5');
    xticks([10, 20, 30, 40, 50, 100, 150, 200, 300, 400, 500]);
    xticklabels({'10','20','30','40','50','100','150','200','300','400','500'});
    grid on
    hold on
end







