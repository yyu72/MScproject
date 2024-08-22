% Define the patient list and conditions
patients = {'Patient.6', 'Patient.9', 'Patient.31', 'Patient.32', 'Patient.72', 'Patient.78', 'Patient.104', 'Patient.105', 'Patient.106'};
conditions = {'N=[]', 'N=1e4', 'N=1e6'};
numPatients = length(patients);
numConditions = length(conditions);
data = cell(numPatients, numConditions);
basePaths = {'/home/yi/Videos/RESULT2/RESULT', '/home/yi/Videos/RESULT/RESULT', '/home/yi/Videos/N0=1e6'};
outputFolder = '/home/yi/Videos/muscatterplot';

% Create output folder if it doesn't exist
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Loop through each patient's data for each condition
means = zeros(numPatients, numConditions);
medians = zeros(numPatients, numConditions);

for i = 1:numPatients
    for j = 1:numConditions
        patientNumber = strrep(patients{i}, 'Patient.', '');
        filename = fullfile(basePaths{j}, patients{i}, ['pt' patientNumber], ['pt' patientNumber '_OptimumParameters.mat']);
        
        if exist(filename, 'file')
            load(filename, 'mu'); % Load the 'mu' variable
            disp(['Data loaded for ', patients{i}, ' under ', conditions{j}]);
            
            % Remove NaN values and ensure all values are positive
            mu_clean = mu(~isnan(mu) & mu > 0);

            % Further clean data by removing very small values close to zero
            mu_clean = mu_clean(mu_clean > 1e-10);

            % Apply logarithmic transformation if there are valid entries
            if ~isempty(mu_clean)
                mu_log_transformed = log(mu_clean);
                data{i, j} = mu_log_transformed;
                means(i, j) = mean(mu_log_transformed);
                medians(i, j) = median(mu_log_transformed);
            else
                disp(['No valid data to transform for ', patients{i}, ' under ', conditions{j}]);
                data{i, j} = [];
                means(i, j) = NaN;
                medians(i, j) = NaN;
            end
        else
            disp(['Warning: File not found for ', patients{i}, ' under ', conditions{j}]);
            data{i, j} = [];
            means(i, j) = NaN;
            medians(i, j) = NaN;
        end
    end
end

% Plot data for each patient under different conditions
colors = lines(numConditions);

for i = 1:numPatients
    figure;
    hold on;
    for j = 1:numConditions
        if ~isempty(data{i, j})
            jittered_x = j + 0.1 * randn(size(data{i, j}));
            scatter(jittered_x, data{i, j}, 5, 'filled', 'MarkerFaceColor', colors(j, :), 'MarkerEdgeColor', colors(j, :), 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
            % Calculate and plot mean and median
            mean_val = means(i, j);
            median_val = medians(i, j);
            plot(j, mean_val, 'kd', 'MarkerFaceColor', 'k'); % mean as black diamond
            plot(j, median_val, 'ks', 'MarkerFaceColor', 'r'); % median as red square
            % Display mean and median values
            text(j, mean_val + 1, sprintf('Mean: %.2f', mean_val), 'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', 'k');
            text(j, median_val - 1, sprintf('Median: %.2f', median_val), 'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', 'r');
        end
    end
    hold off;
    title(['Log Transformed Mutation Rates for ', patients{i}, ' under Different Conditions']);
    xlabel('Condition');
    ylabel('Log(Mutation Rate)');
    set(gca, 'XTick', 1:numConditions, 'XTickLabel', conditions);
    legend([conditions, 'Mean', 'Median'], 'Location', 'Best');
    grid on;
    ylim([-20, 5]); % Adjust y-axis limits for better focus on data range
    % Save figure
    saveas(gcf, fullfile(outputFolder, [patients{i} '_mu_scatter_plot.png']));
end

% Display mean and median values for all patients and conditions
disp('Mean and Median values for all patients and conditions:');
disp(table(patients', means, medians, 'VariableNames', {'Patient', 'Means', 'Medians'}));
