matFilePaths = {
    '/home/yi/Videos/mu0=-1e4/Patient.6/pt6/pt6_OptimumParameters.mat',
    '/home/yi/Videos/mu0=-1e6/Patient.6/pt6/pt6_OptimumParameters.mat'
};

analysisResults = struct();

for i = 1:length(matFilePaths)
    matFilePath = matFilePaths{i};
    
    data = load(matFilePath);
    
    if isfield(data, 'N')
        Ne = data.N;
  
        Ne = Ne(~isnan(Ne));  
        Ne = Ne(Ne <= 10^8);  
        
        if isempty(Ne)
            warning('Data empty：%s', matFilePath);
            continue;
        end
        

        medianNe = median(Ne);
        

        figure;
        edges = logspace(log10(min(Ne)), log10(max(Ne)), 50); 
        histogram(Ne, edges, 'DisplayStyle', 'bar');
        set(gca, 'XScale', 'log', 'YScale', 'log'); 
        title('Log Scale Distribution of Estimate N for Patient 6');
        xlabel('Ne');
        ylabel('Frequency');
        grid on;
        

        hold on;
        yLimits = ylim;  
        plot([medianNe, medianNe], yLimits, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Median Ne'); 
        scatter(medianNe, yLimits(1) + 0.05 * diff(yLimits), 100, 'r', 'filled', 'DisplayName', 'Median Ne');  
        legend('show');
        hold off;
        
        fprintf('中位数Ne for %s: %g\n', matFilePath, medianNe);
        
        analysisResults(i).filePath = matFilePath;
        analysisResults(i).medianNe = medianNe;
    else
        warning('文件中找不到变量N：%s', matFilePath);
    end
end


disp('Finished');
