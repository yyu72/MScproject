function [xx_2d_pos_cell, time_points_cell]=WF_single_timeseries(timesrun,N,mu,s,fr,ns,plotfig)
%this script is for generating single site frequency based on WF model with
%selection, mutation and drift.

%For initial variables, 
%timesrun is for setting how many time of run based on you wish
%N-population size
%mu-mutation rate
%s-selection coefficient
%ns-binomial sampling sites
%plotfig- set to 1 to enable frequency time series plot, 0 disable

%dt_list  list of sampling time interval
%T_list   list of total time of generation
%x0_list  list of initial frequency


dt_list=[400,500];
T_list=10000;
x0_list=0.02;

num_x0 = length(x0_list);
num_T = length(T_list);
num_dt = length(dt_list);
    
    data_cell = cell(num_x0, num_dt);
    
    for x0i = 1:num_x0
        for dti = 1:num_dt
            xx_2d_pos_cell = cell(1, timesrun);
            time_points_cell = cell(1, timesrun);
            
            for run_idx = 1:timesrun
                [xx_2d_pos, time_points] = WF_single_simscp_time(N, mu, s, x0_list(x0i), T_list, ns, dt_list(dti), plotfig);
                xx_2d_pos_cell{run_idx} = xx_2d_pos;
                time_points_cell{run_idx} = time_points;
            end
            
            data_cell{x0i, dti} = struct('x0', x0_list(x0i), 'dt', dt_list(dti), 'xx_2d_pos', xx_2d_pos_cell, 'time_points', time_points_cell);
            
            %change the save path on you own job path
            save_path='C:\Users\cmcdd\Desktop\project\results\single site results\negative selection\timeseries';
            save(fullfile(save_path, ['x0_', num2str(x0_list(x0i)), '_dt_', num2str(dt_list(dti)), '_data.mat']), 'xx_2d_pos_cell', 'time_points_cell');
            disp(['Saved data for x0 = ', num2str(x0_list(x0i)), ' and dt = ', num2str(dt_list(dti))]);
        end
    end
    
    save(fullfile(save_path, 'all_data.mat'), 'data_cell');
    disp('Saved all_data.mat');
end

