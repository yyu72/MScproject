function analyse_and_save_parameters_loop(data_path, save_path,I, ns, Dt, N0, mu0, samplcorrections, initialparams, sinitopt, null, fine, prior, plotfig, suppressoutput, CI, boundary, optimmethod)
  %This script is used for single site analysis and parameters
  %optimisation.

  %data_path-the path that saves your time series data in .mat
  %save_path-the path that you want to save the optimised parameters
  %rest of them are same for setting 'GenomeScanMutation.m'


% Load the all_data.mat file
    load(fullfile(data_path, 'all_data.mat'), 'data_cell');
    
    num_x0 = size(data_cell, 1); % Number of x0 values
    num_dt = size(data_cell, 2); % Number of dt values
    
    for x0i = 1:num_x0
        for dti = 1:num_dt
            cell_data = data_cell{x0i, dti}; % Get the cell data for the current x0 and dt
            
            Nopt_list = [];
            sopt_list = [];
            muopt_list = [];
            
            for run_idx = 1:numel(cell_data)
                xx_2d_pos = cell_data(run_idx).xx_2d_pos; % Get xx_2d_pos from the current run
                time_points = cell_data(run_idx).time_points; % Get time_points from the current run
                
                % Perform parameter analysis here to get Nopt, sopt, muopt
                [Nopt, sopt, muopt] = analyse_parameters(time_points,xx_2d_pos, I, ns, Dt, N0, mu0, samplcorrections, initialparams, sinitopt, null, fine, prior, plotfig, suppressoutput, CI, boundary, optimmethod);
                
                Nopt_list(run_idx) = Nopt;
                sopt_list(run_idx) = sopt;
                muopt_list(run_idx) = muopt;
            end
            
            % Save the lists as a .mat file for each combination of x0 and dt
            output_filename = sprintf('x0_%g_dt_%d_analysis.mat', cell_data(1).x0, cell_data(dti).dt);
            save(fullfile(save_path, output_filename), 'Nopt_list', 'sopt_list', 'muopt_list');
            disp(['Saved parameter analysis for x0 = ', num2str(cell_data(1).x0), ' and dt = ', num2str(cell_data(1).dt)]);
        end
    end
end

function [Nopt, sopt, muopt] = analyse_parameters(time_points,xx_2d_pos, I, ns, Dt, N0, mu0, samplcorrections, initialparams, sinitopt, null, fine, prior, plotfig, suppressoutput, CI, boundary, optimmethod)
    % Perform your parameter analysis here
    % Replace the following lines with your actual analysis
    [Nopt, sopt, muopt, ~, ~, ~, ~, ~, ~] = GenomeScanMutation(time_points, xx_2d_pos, I, ns, Dt, N0, mu0, samplcorrections, initialparams, sinitopt, null, fine, prior, plotfig, suppressoutput, CI, boundary, optimmethod);
end






