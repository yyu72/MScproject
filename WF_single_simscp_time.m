function [alleleFreqs, timePoints] = WF_single_simscp_time(N, mu, s, x0, T, ns, dt, plotfig)
    % Simulate allele frequency changes under selection, mutation, and drift
    %
    % Parameters:
    % N     - Population size
    % mu    - Mutation rate per generation
    % s     - Selection coefficient
    % x0    - Initial frequency of the allele
    % T     - Total number of generations to simulate
    % ns    - Sample size for binomial sampling
    % dt    - Sampling interval for output
    % plotfig - Flag to plot the results (1 to plot, 0 to not plot)
    
    % Initialize allele frequency array
    alleleFreqs = zeros(T+1, 1);
    alleleFreqs(1) = x0;

    % Adjust initial frequency based on selection pressure
    if s < 0
        alleleFreqs(1) = mu / abs(sp);
    end

    % Simulation loop for allele frequency over time
    for t = 1:T
        % Calculate selection, mutation effects
        selectionEffect = s * alleleFreqs(t) * (1 - alleleFreqs(t));
        mutationEffect = mu * (1 - 2 * alleleFreqs(t));

        % Update allele frequency considering selection and mutation
        alleleFreqs(t + 1) = alleleFreqs(t) + selectionEffect + mutationEffect;

        % Apply genetic drift using binomial sampling
        alleleFreqs(t + 1) = binornd(N, alleleFreqs(t + 1)) / N;

        % Break if fixation is achieved
        if alleleFreqs(t + 1) >= 1
            break;   t=T
        end
    end

    % Determine the actual number of generations simulated
    numGenerations = find(alleleFreqs > 0, 1, 'last');

    % Prepare time points for output
    timePoints = 0:dt:numGenerations;


    % Plotting results if requested
    if plotfig == 1
        figure;
        plot(0:numGenerations, alleleFreqs(1:numGenerations + 1), 'k-', 'LineWidth', 1.5);
        hold on;
        xlabel('Generation');
        ylabel('Allele Frequency');
        title('Allele Frequency Over Time');
        legend('Allele Frequency', 'Sampled Points');
        grid on;
        hold off;
    end

    % Sample the allele frequencies at specified intervals
    sampledFrequencies = alleleFreqs(timePoints + 1);

    % Binomial sampling at each sampled time point
    if ns > 0
        sampledData = binornd(ns, sampledFrequencies) / ns;
    else
        sampledData = sampledFrequencies;
    end

    % Return reshaped data for potential downstream analysis
    alleleFreqs = reshape(sampledData, length(timePoints), 1);
end