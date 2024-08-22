function [data_cell, save_path] = WF_single_timeseries(timesRun, N, mu, s, fr, ns, plotFig)
    % Simulate single site allele frequency time series multiple times and store the results.
    %
    % Parameters:
    % timesRun - Number of simulations to run
    % N - Population size
    % mu - Mutation rate per generation
    % s - Selection coefficient
    % fr - Initial frequency of the allele
    % ns - Sample size for binomial sampling
    % plotFig - Set to 1 to enable frequency time series plot, 0 to disable
    
    dt_list = [400, 500];  % List of sampling time intervals
    T = 10000;             % Total time of generation for each simulation
    x0_list = [0.02, 0.05, 0.1];  % List of initial frequencies

    % Initialize cell array to store data from all configurations
    data_cell = cell(length(x0_list), length(dt_list));
    
    for x0_idx = 1:length(x0_list)
        for dt_idx = 1:length(dt_list)
            % Initialize cells to store simulation results
            frequencyDataCell = cell(1, timesRun);
            timePointsCell = cell(1, timesRun);
            
            for run_idx = 1:timesRun
                [sampledFrequencies, timePoints] = WF_single_simscp_time(N, mu, s, x0_list(x0_idx), T, ns, dt_list(dt_idx), plotFig);
                frequencyDataCell{run_idx} = sampledFrequencies;
                timePointsCell{run_idx} = timePoints;
            end
            
            % Store results for this configuration
            data_cell{x0_idx, dt_idx} = struct('x0', x0_list(x0_idx), 'dt', dt_list(dt_idx), 'frequencies', frequencyDataCell, 'timePoints', timePointsCell);
            
            % Define save path
            save_path = fullfile('C:', 'Users', 'cmcdd', 'Desktop', 'project', 'results', 'single site results', 'negative selection', 'timeseries');
            if ~exist(save_path, 'dir')
                mkdir(save_path);
            end

            % Save the data for this configuration
            save_name = sprintf('x0_%g_dt_%d_data.mat', x0_list(x0_idx), dt_list(dt_idx));
            save(fullfile(save_path, save_name), 'frequencyDataCell', 'timePointsCell');
            disp(['Saved data for x0 = ', num2str(x0_list(x0_idx)), ' and dt = ', num2str(dt_list(dt_idx))]);
        end
    end

    % Save aggregated data from all configurations
    save(fullfile(save_path, 'all_data.mat'), 'data_cell');
    disp('Saved all_data.mat');
end
