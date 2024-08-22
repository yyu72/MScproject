basePath = '/home/yi/Documents/MScproject/data';


dirs = dir(basePath);
isSub = [dirs(:).isdir];  
subFolders = {dirs(isSub).name}';  
subFolders(ismember(subFolders, {'.', '..'})) = [];  


pattern = '^(Patient\.|Linked )';
matchedFolders = subFolders(cellfun(@(x) ~isempty(regexp(x, pattern, 'once')), subFolders));


diversityResults = struct();


colors = lines(length(matchedFolders)); 


for i = 1:length(matchedFolders)
    folderName = matchedFolders{i};
    matFilePath = fullfile(basePath, folderName, [folderName '_AlleleFrequencies+snp_table.mat']);


    if exist(matFilePath, 'file')
        data = load(matFilePath);
        if isfield(data, 'f') && isfield(data, 't') && size(data.t, 1) > 3
            f = data.f;  
            T = size(f, 3);  
            L = size(f, 1);  
            nnuc = size(f, 2);  
            
           
            sumDiversity = zeros(1, T);
            meanDiversity = zeros(1, T);
            
            
            for t = 1:T
                diversityAtTimePoint = zeros(L, 1);
                for l = 1:L
                    if any(f(l, 5:6, t) > 0)  
                        continue;  
                    end
                    
                   
                    for i = 1:nnuc
                        for j = i+1:nnuc
                            diversityAtTimePoint(l) = diversityAtTimePoint(l) + ...
                                2 * f(l, i, t) * f(l, j, t);
                        end
                    end
                end
                
               
                sumDiversity(t) = sum(diversityAtTimePoint);
                meanDiversity(t) = mean(diversityAtTimePoint(diversityAtTimePoint > 0));  
            end
            
           
            diversityResults.(strrep(folderName, '.', '_')).sumDiversity = sumDiversity;
            diversityResults.(strrep(folderName, '.', '_')).meanDiversity = meanDiversity;
            diversityResults.(strrep(folderName, '.', '_')).timePoints = data.t; 
        end
    end
end


save(fullfile(basePath, 'diversityResults.mat'), 'diversityResults');

maxTimePoints = 0;
for i = 1:length(matchedFolders)
    folderName = strrep(matchedFolders{i}, '.', '_');
    if isfield(diversityResults, folderName)
        maxTimePoints = max(maxTimePoints, max(diversityResults.(folderName).timePoints));
    end
end


figure;

subplot(2,1,1);
hold on;
title('Mean Nucleotide Diversity (0-50 Days)');
for i = 1:length(matchedFolders)
    folderName = strrep(matchedFolders{i}, '.', '_');
    if isfield(diversityResults, folderName)
        validIndices = diversityResults.(folderName).timePoints <= 50;
        plot(diversityResults.(folderName).timePoints(validIndices), ...
            diversityResults.(folderName).meanDiversity(validIndices), '-o', 'Color', colors(i, :), 'DisplayName', folderName);
    end
end
legend show;
xlabel('Time (Days)');
ylabel('Mean Diversity');
set(gca, 'YScale', 'log');
grid on;
hold off;

subplot(2,1,2);
hold on;
title(['Mean Nucleotide Diversity (50+ Days) up to ' num2str(maxTimePoints) ' Days']);
for i = 1:length(matchedFolders)
    folderName = strrep(matchedFolders{i}, '.', '_');
    if isfield(diversityResults, folderName)
        validIndices = diversityResults.(folderName).timePoints > 50;
        if any(validIndices)
            lastIndex = find(diversityResults.(folderName).timePoints <= 50, 1, 'last');
            if ~isempty(lastIndex)
                plot([diversityResults.(folderName).timePoints(lastIndex) diversityResults.(folderName).timePoints(find(validIndices, 1))], ...
                    [diversityResults.(folderName).meanDiversity(lastIndex) diversityResults.(folderName).meanDiversity(find(validIndices, 1))], ...
                    '--', 'Color', colors(i, :), 'HandleVisibility', 'off');
            end
            
            plot(diversityResults.(folderName).timePoints(validIndices), ...
                diversityResults.(folderName).meanDiversity(validIndices), '-*', 'Color', colors(i, :), 'DisplayName', folderName);
        end
    end
end
legend show;
xlabel('Time (Days)');
ylabel('Mean Diversity');
set(gca, 'YScale', 'log');
grid on;
hold off;
