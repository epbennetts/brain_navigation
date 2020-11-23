% %test params
% clear all;
% close all force;
% %see how to make SetTestParams() work without having to do area = %params.area, etc. 
% params = SetParams_AccVsNumGenes(); 
% 
% %"Translate" params (do this better)
% area = params.area;
% costFunction = params.costFunction;
% 
% cols = params.cols;
% prevBestGenes = params.prevBestGenes;
% noiseStDev = params.noiseStDev;
% numNoiseIterations = params.numNoiseIterations;
% numFolds = params.numFolds;
% %plotting
% %numplots = params.numplots;
% %range = params.range;
% %how many genes considered in this run:
% samples = 5;
% numplots = 5;
% range = 'top';
% sizeSampleSubset = 5;
% numGenes = 10;

function [indexOrder, geneNames_ranked, balAccuracies_ranked, confMatrices_ranked, trees_all_clean, balAcc_errors_ranked] = DT_classification_multiple(sizeSampleSubset, area, prevBestGenes, noiseStDev, numFolds, numNoiseIterations)
%
%Finds the gene within the sample that performs best at classifying the target areas with a DT, in combination with prev_best_genes.
%Ranks genes according to performance (balanced accuracy).
%Will also plot a histogram of the accuracies
%
%INPUTS:
%samples: size of the gene subset to complete the algorithm on (can choose samples < total #genes for testing)
%prev_best_genes: indices of the best genes obtained from previous iterations. 

%noise set
%here or in general?
rng('default');

%SET VARS
numGenesInDT = length(prevBestGenes) + 1;
[geneData, isTarget, geneNames, structInfo] = filter_nans(area);
%rows == #areas,  cols = #genes in the data
[rows, cols] = size(geneData);
classes = isTarget;
%empty arrays
balAccuracies_avg = zeros(sizeSampleSubset,1);
confMatrices_avg = zeros(sizeSampleSubset,4);
trees_all_clean = cell(sizeSampleSubset,1);
balAcc_errors = zeros(sizeSampleSubset,1);
predictedLabels_all = NaN(sizeSampleSubset,rows,numNoiseIterations);

% %test
fprintf('**SAMPLES = %d**\n', sizeSampleSubset)
fprintf('NEXT CLASSIFICATION ITERATION \n')
fprintf('(printed from DT_class) num genes considered in DT: %d \n\n', numGenesInDT)

%set cost function
costFunc = ComputeBalancedCostFunc(classes);
        
%MAKE FOLDS
partitions = cell(numNoiseIterations,1);
%MAKE NOISE MATRICES (here or in general?)
geneNoise = [];
%NaN(rows, numGenesInDT, numNoiseIterations);
for i = 1:numNoiseIterations
    geneNoiseCol = noiseStDev*randn(rows,1);    
    geneNoiseTemp = repmat(geneNoiseCol,1,numGenesInDT);
    geneNoise = cat(3, geneNoise, geneNoiseTemp);
end

for part = 1:numNoiseIterations
    %make stratified CV folds
    partitions{part} = cvpartition(classes,'KFold',numFolds,'Stratify',true);
end        

% Previously selected genes:
    prevGeneData = geneData(:, prevBestGenes);

%loop over genes or subset, classify and evaluate metrics
% samples == num genes
%parpool('local', 12)

confMatrices_iter = nan(numNoiseIterations,4);
balAccuracies_iter = nan(numNoiseIterations,1);
parfor i = 1:sizeSampleSubset
    %fprintf('NEXT GENE (%d) \n', i)
        
    confMatrices_iter = nan(numNoiseIterations,4);
    balAccuracies_iter = nan(numNoiseIterations,1);
    
    %SELECT GENE DATA
    % Set gene data for this iteration:
    geneCombo = [prevGeneData geneData(:,i)];
    
    for iter = 1:numNoiseIterations
        %predictedLabels_all(:,iter) = Nans(1,iter);
        
        %fprintf('noise iteration: %d (printed from DT_class) \n', iter);        
        noiseComboIter = geneNoise(:,:,iter);        
        
        % TRAINING AND TESTING
        %Iterating through CV folds, random noise in testing 
        [predictedLabels] = kFoldPredictNoisy(geneCombo, noiseComboIter, rows, classes, numGenesInDT, noiseStDev, costFunc, partitions{iter}, numFolds);
%         %save in external function because of parfor
%         doSaveLabels('Labels_LatestIter',predictedLabels);
         predictedLabels_all(i,:,iter) = predictedLabels;
        
        % Across the folds, predicted class labels are filled into predictedLabels
        % classes == real labels
        [confMatrices_iter(iter,:), balAccuracies_iter(iter)] = ComputeConfusion(classes,predictedLabels);
     
    end
    
    % COMPUTE AND SAVE IDEAL TREES (no CV or noise)
    tree_clean = fitctree(geneCombo, classes, 'MaxNumSplits', numGenesInDT, 'cost', costFunc, 'ClassNames', [1,0]);
    trees_all_clean{i,1} = tree_clean;
    
    %conf matrices
    confMatrix_avg = mean(confMatrices_iter, 1);
    confMatrices_avg(i,:) = confMatrix_avg;
    %bal accuracies
    balAccuracies_avg(i) = mean(balAccuracies_iter);
    %std deviations of the accuracies
    balAcc_error = std(balAccuracies_iter);
    balAcc_errors(i) = balAcc_error;
    
end

%sort accuracies
[balAccuracies_ranked, indexOrder] = metricSort(balAccuracies_avg, 'descend');
confMatrices_ranked = confMatrices_avg(indexOrder, :);
geneNames_ranked = geneNames(indexOrder);
balAcc_errors_ranked = balAcc_errors(indexOrder);

save('Workspace_Latest_Iteration.mat')



%-------------------------------------------------------------------------------
% PLOTTING (make external)
%-------------------------------------------------------------------------------

% %HISTOGRAMS
% numBins = 10;
% plotAccuracyHistogram(balAccuracies_ranked, numBins, numGenesInDT, area, sizeSampleSubset)
% 
%  
%  
% %plotting vars:
% prevGeneData = geneData(:, prevBestGenes);
% thresholds_all_clean = nan(samples, numGenes);
% %non/target indices (to plot them differently)
% targetIndices = isTarget;
% nonTargetIndices = find(~isTarget);
% target = geneData(targetIndices, :);
% nonTarget = geneData(nonTargetIndices, :);
% 
% %various plotting thresholds
% if strcmp(range,'top')
%     %top n
%     orderedRange = indexOrder(1:numplots)';
% elseif strcmp(range, 'middle')
%     middleStart = samples*0.2;
%     middleEnd = middleStart + numplots - 1;
%     %middle n
%     orderedRange = indexOrder(middleStart:middleEnd)';
% elseif strcmp(range, 'bottom')
%     bottomStart = samples - numplots + 1;
%     bottomEnd = samples;
%     %bottom n
%     orderedRange = indexOrder(bottomStart:bottomEnd)';
% end
% 
% %make plots for given accuracy range
% for i = orderedRange
%     
%     %-------------------------------------------------------------------------------
%     % COMPUTE THRESHOLDS FROM CLEAN DATA
%     %-------------------------------------------------------------------------------    
%     %fetch ideal tree (no CV or noise)
%     tree_clean = trees_all_clean{i};
%     %compute thresholds
%     thresholds_raw = tree_clean.CutPoint; %threshold with nans
%     thresholds = thresholds_raw(~isnan(thresholds_raw));
%     thresholds_pan = nan(1,numGenes);
%     thresholds_pan(1:length(thresholds)) = thresholds;
%     thresholds_all_clean(i,:) = thresholds_pan;
%     
%     %get predictor order
%     genesUsed = tree_clean.CutPredictor;    
%     
%     %-------------------------------------------------------------------------------
%     % PLOT TREE GRAPHS
%     %-------------------------------------------------------------------------------
%     
%     view(tree_clean,'Mode','graph')
%     figure;    
%     
%     %-------------------------------------------------------------------------------
%     % PLOT GENE EXPRESSIONS
%     %-------------------------------------------------------------------------------
%     
%     t1 = thresholds_all_clean(i,1);
%     
%     %------------------------
%     % 1D HISTOGRAMS
%     if (numGenes == 1)
%         histogram(target(:,i), 'Normalization', 'count', 'BinWidth', 0.1);
%         title(sprintf('%s : Gene %d with accuracy %g', geneNames{i}, i, balAccuracies(i)));
%         if strcmp(geneNames_ranked(i), geneNames(indexOrder(i)))
%             sprintf("Error: geneNames_ranked(i) is %s and geneNames(indexOrder(i) is %s.", geneNames_ranked{i}, geneNames{indexOrder(i)});
%         end
%         hold on;
%         histogram(nonTarget(:,i), 'Normalization', 'count', 'BinWidth', 0.1);
%         if ~isnan(t1)
%             xline(t1, '--r');
%         end
%         hold off;
%         legend({'target', '~target'});
%         
%     %------------------------    
%     % 2D PLOTS
%     elseif (numGenes == 2)
%         
%         %Plot gene expression data points
%         hold on;        
%         for a = 1:rows
%             %get area colours
%             colour = rgbconv(structInfo{a,2}{1});
%                      
%             %average deviation (i.e. estimate of the average absolute noise)
%             %noiseLevel == standard deviation
%             avgNoise = 0.7979*noiseStDev;
%             
%             %plot data points, formatting corresponding to area and class 
%             if isTarget(a)
%                 x = geneData(a,prevBestGenes);
%                 y = geneData(a,i);
%                 errorbar(x, y, avgNoise, 'both', '*', 'color', colour)
%             else
%                 x = geneData(a,prevBestGenes);
%                 y = geneData(a,i);
%                 errorbar(x, y, avgNoise, 'both', '.', 'color', colour)
%             end            
%         end
%         
%         %make square plot
%         axis equal
%         %axis limits
%         xlim([min(geneData(:,prevBestGenes)) - 0.1, max(geneData(:,prevBestGenes)) + 0.1])
%         ylim([min(geneData(:,i)) - 0.1, max(geneData(:,i)) + 0.1])        
%         
%         %plot thresholds
%         if strcmp(genesUsed{1,1}, 'x1')
%             xline(t1, '--');
%             if ~isnan(thresholds_all_clean(i,2))
%                 t2 = thresholds_all_clean(i,2);
%                 %this assumes that if numgenes = 2, prevGeneData is just 1 gene
%                 line([t1,max(prevGeneData)], [t2,t2], 'LineStyle', '--', 'Color', 'k')                
%             end
%             
%         elseif strcmp (genesUsed{1,1}, 'x2')
%             yline(t1, '--');
%             if ~isnan(thresholds_all_clean(i,2))
%                 t2 = thresholds_all_clean(i,2);
%                 %this assumes that if numgenes = 2, prevGeneData is just 1 gene
%                 line([t2,t2],[t1, max(prevGeneData)], 'LineStyle', '--', 'Color', 'k')
%             end
%         end
%         hold off;
%         title(sprintf('Genes %d and %d with accuracy %g', prevBestGenes, i, balAccuracies(i)));
%         %legend({'target', 'threshold', '~target'});
%         
%     end
end