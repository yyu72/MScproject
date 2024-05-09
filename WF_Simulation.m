function [alleleFreqs, selectionCoeffs, types] = WF_Simulation(N, mu, sn, sp, fr_neg, fr_pos, L, T, dt)
    % Input Parameters:
    % N - Population size
    % mu - Mutation rate per generation per locus
    % sn - Mean strength of negative selection
    % sp - Mean strength of positive selection
    % fr_neg - Fraction of loci under negative selection
    % fr_pos - Fraction of loci under positive selection
    % L - Total number of loci
    % T - Number of generations to simulate
    % dt - Sampling interval for output

    % Initialize arrays
    alleleFreqs = zeros(T+1, L);    % Allele frequencies over time for each locus
    selectionCoeffs = zeros(L, 1);  % Selection coefficients for each locus
    types = strings(1, L);          % Type of selection for each locus

    % Set initial allele frequencies
    % Negative selection
    alleleFreqs(1, 1:round(fr_neg*L)) = gamrnd(2 * N * mu, 1/(2 * N * sn), [1, round(fr_neg*L)]);
    % Positive selection
    alleleFreqs(1, round(fr_neg*L)+1:round((fr_neg+fr_pos)*L)) = rand(1, round(fr_pos*L));
    % Neutral loci
    alleleFreqs(1, round((fr_neg+fr_pos)*L)+1:end) = rand(1, L - round((fr_neg+fr_pos)*L));

    % Ensure initial frequencies are within [0, 1]
    alleleFreqs(1, :) = min(max(alleleFreqs(1, :), 0), 1);

    % Assign selection coefficients and classify loci type
    for i = 1:L
        if i <= round(fr_neg*L)
            selectionCoeffs(i) = -abs(normrnd(sn, 0.1));  % Negative selection
            types(i) = "Negative";
        elseif i >= round((fr_neg+fr_pos)*L)
            selectionCoeffs(i) = abs(normrnd(sp, 0.1));   % Positive selection
            types(i) = "Positive";
        else
            selectionCoeffs(i) = 0;                       % Neutral
            types(i) = "Neutral";
        end
        
        % Check if locus is effectively neutral
        if 2 * N * abs(selectionCoeffs(i)) <= 1
            types(i) = "Neutral";  % Redefine as neutral if selection is weak
        end
    end

    % Simulation loop for T generations
    for t = 1:T
        for i = 1:L
            % Calculate selection and mutation effects
            selectionEffect = selectionCoeffs(i) * alleleFreqs(t, i) * (1 - alleleFreqs(t, i));
            mutationEffect = mu * (1 - 2 * alleleFreqs(t, i));

            % Update allele frequency
            alleleFreqs(t+1, i) = alleleFreqs(t, i) + selectionEffect + mutationEffect;

            % Apply genetic drift
            alleleFreqs(t+1, i) = binornd(N, alleleFreqs(t+1, i)) / N;

            % Ensure frequencies are within [0, 1]
            alleleFreqs(t+1, i) = min(max(alleleFreqs(t+1, i), 0), 1);
        end
    end

    % Plot results for a sample locus of each type
    plotLociResults(alleleFreqs, types, L);
end

function plotLociResults(alleleFreqs, types, L)
    % Find sample loci for each type for plotting
    neg_locus = find(types == "Negative", 1);
    pos_locus = find(types == "Positive", 1);
    neut_locus = find(types == "Neutral", 1);

    figure;
    subplot(1,3,1); % Negative selection
    plot(0:size(alleleFreqs,1)-1, alleleFreqs(:, neg_locus));
    title(sprintf('Negative Selection at Locus %d', neg_locus));
    xlabel('Generation');
    ylabel('Allele Frequency');

    subplot(1,3,2); % Positive selection
    plot(0:size(alleleFreqs,1)-1, alleleFreqs(:, pos_locus));
    title(sprintf('Positive Selection at Locus %d', pos_locus));
    xlabel('Generation');
    ylabel('Allele Frequency');

    subplot(1,3,3); % Neutral
    plot(0:size(alleleFreqs,1)-1, alleleFreqs(:, neut_locus));
    title(sprintf('Neutral at Locus %d', neut_locus));
    xlabel('Generation');
    ylabel('Allele Frequency');
end


