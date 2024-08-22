% Define the patient list and conditions
patients = {'Patient.6', 'Patient.9', 'Patient.31', 'Patient.32', 'Patient.72', 'Patient.78', 'Patient.104', 'Patient.105', 'Patient.106'};
conditions = {'N=1e4', 'N=1e6', 'mu=1e-4', 'mu=1e-6'};
numPatients = length(patients);
numConditions = length(conditions);
data_original = cell(numPatients, numConditions);
means = zeros(numPatients, numConditions);
medians = zeros(numPatients, numConditions);
basePaths = {
    '/home/yi/Videos/RESULT/RESULT',
    '/home/yi/Videos/N0=1e6',
    '/home/yi/Videos/mu0=-1e4',
    '/home/yi/Videos/mu0=-1e6'
};
outputFolder = '/home/yi/Videos/sscatterplot';

% Create output folder if it doesn't exist
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Loop through each patient's data for each condition
for i = 1:numPatients
    for j = 1:numConditions
        patientNumber = strrep(patients{i}, 'Patient.', '');
        ptNumber = ['pt' patientNumber];
        filename = fullfile(basePaths{j}, patients{i}, ptNumber, [ptNumber '_OptimumParameters.mat']);
        
        if exist(filename, 'file')
            load(filename, 'sopt'); % Load the 'sopt' variable
            disp(['Data loaded for ', patients{i}, ' under ', conditions{j}]);
            
            % Remove NaN values 
            s_clean = sopt(~isnan(sopt));
            
            disp(['Number of entries after NaN filter: ', num2str(length(s_clean))]);

            % Store original values
            data_original{i, j} = s_clean;
            
            if ~isempty(s_clean)
                means(i, j) = mean(s_clean);
                medians(i, j) = median(s_clean);
            else
                means(i, j) = NaN;
                medians(i, j) = NaN;
            end
        else
            disp(['Warning: File not found for ', patients{i}, ' under ', conditions{j}]);
            data_original{i, j} = [];
            means(i, j) = NaN;
            medians(i, j) = NaN;
        end
    end
end

% Plot box plots for each patient under different conditions with log Y-axis
for i = 1:numPatients
    figure;
    hold on;
    
    % Prepare data for boxplot
    data_for_boxplot = [];
    group_labels = [];
    
    for j = 1:numConditions
        data_for_boxplot = [data_for_boxplot; data_original{i, j}];
        group_labels = [group_labels; repmat(conditions(j), length(data_original{i, j}), 1)];
    end
    
    % Plot the boxplot
    boxplot(data_for_boxplot, group_labels);
    
    % Set Y-axis to log scale
    set(gca, 'YScale', 'log');
    
    title(['Log-Scaled Selection Coefficients for ', patients{i}]);
    ylabel('Log-Scaled Selection Coefficient');
    xlabel('Condition');
    grid on;
    hold off;
    
    % Save the figure
    saveas(gcf, fullfile(outputFolder, [patients{i} '_selection_box_plot_log.png']));
end

% Save mean and median values to a .mat file
save(fullfile(outputFolder, 'mean_median_sopt_values.mat'), 'means', 'medians', 'patients', 'conditions');

% Display mean and median values for all patients and conditions
disp('Mean and Median Selection Coefficients for all patients and conditions:');
disp(table(patients', means, medians, 'VariableNames', {'Patient', 'Means', 'Medians'}));

disp('Log-scaled box plots and mean/median values have been generated and saved.');
