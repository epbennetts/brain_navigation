% %dummy params
% clear all;
% close all force;
% params = SetDefaultParams();
% 
% cols = 19114;
% samples = 20;
% area = "Isocortex";
% %prev_best_genes = [4307, 636, 148]; %for noise
% prevBestGenes = [3725];
% noiseStDev = 1;


function [indexOrder, balAcc_ranked, confMatrices_ranked, geneNames_ranked, trees_all_clean] = DT_classification_multiple(samples, area, prevBestGenes, noiseStDev)
%Finds the gene within the sample that performs best at classifying the target areas with a DT, in combination with prev_best_genes.
%Ranks genes according to performance (balanced accuracy).
%Will also plot a histogram of the accuracies
%
%INPUTS:
%samples: size of the gene subset to complete the algorithm on (can choose samples < total #genes for testing)
%prev_best_genes: indices of the best genes obtained from previous iterations. 


%SET VARS
numGenes = length(prevBestGenes) + 1;
[genes, isTarget, geneNames, structInfo] = filter_nans(area);
%rows == #areas,  cols = #genes in the data
[rows, cols] = size(genes);
classes = isTarget;
%empty arrays
balAccuracies = zeros(samples,1);
confMatrices = zeros(samples,4);
trees_all_clean = cell(samples,1);

%loop over genes or subset, classify and evaluate metrics
% samples == num genes
for i = 1:samples
    
    numNoiseIterations = 15;
    confMatrices_iter = nan(numNoiseIterations,4);
    balAccuracies_iter = nan(numNoiseIterations,1);
    
    for iter = 1:numNoiseIterations
        
        %SELECT GENE DATA
        % Previously selected genes:
        prevGeneData = genes(:, prevBestGenes);
        % Set gene data for this iteration:
        geneCombo = [prevGeneData genes(:,i)];
        
        %MAKE FOLDS
        %set cost function
        costFunc = ComputeBalancedCostFunc(isTarget);        
        %make stratified CV folds
        numFolds = 10;
        partition = cvpartition(classes,'KFold',numFolds,'Stratify',true);
                
        % TRAINING AND TESTING
        %Iterating through CV folds, random noise in testing 
        [predictedLabels] = kFoldPredictNoisy(geneCombo, rows, classes, numGenes, noiseStDev, costFunc, partition, numFolds);
        
        % Across the folds, predicted class labels are filled into predictedLabels
        % classes == real labels
        [confMatrices_iter(iter,:), balAccuracies_iter(iter)] = ComputeConfusion(classes,predictedLabels);
     
    end
    
    % COMPUTE AND SAVE IDEAL TREES (no CV or noise)
    tree_clean = fitctree(geneCombo, classes, 'MaxNumSplits', numGenes, 'cost', costFunc, 'ClassNames', [1,0]);
    trees_all_clean{i,1} = tree_clean;
    
    %conf matrices
    confMatrix = mean(confMatrices_iter, 1);
    confMatrices(i,:) = confMatrix;
    %bal accuracies
    balAccuracies(i) = mean(balAccuracies_iter);
    
end

%sort accuracies
[balAcc_ranked, indexOrder] = metricSort(balAccuracies, 'descend');
confMatrices_ranked = confMatrices(indexOrder, :);
geneNames_ranked = geneNames(indexOrder);




%-------------------------------------------------------------------------------
% PLOTTING (make external)
%-------------------------------------------------------------------------------

%ACCURACY HISTOGRAMS
figure();
histogram(balAcc_ranked, 'Normalization', 'count', 'NumBins', 10);
title(sprintf('Balanced Accuracy for %g genes in %s (samples = %g)', numGenes, area, samples))
xlabel('Balanced accuracy');
ylabel('Counts');


%plotting vars
range = 'top';
numplots = 5;
thresholds_all_clean = nan(samples, numGenes);

%various plotting thresholds
plots = numplots;
middleStart = samples*0.2;
middleEnd = middleStart + plots - 1;
bottomStart = samples - plots + 1;
bottomEnd = samples;

if strcmp(range,'top')
    %top n
    orderedRange = indexOrder(1:plots)';
elseif strcmp(range, 'middle')
    %middle n
    orderedRange = indexOrder(middleStart:middleEnd)';
elseif strcmp(range, 'bottom')
    %bottom n
    orderedRange = indexOrder(bottomStart:bottomEnd)';
end

%make various plots in a certain accuracy range
targetIndices = isTarget;
nonTargetIndices = find(~isTarget);
target = genes(targetIndices, :);
nonTarget = genes(nonTargetIndices, :);
prevGeneData = genes(:, prevBestGenes);

for i = orderedRange
    
    %-------------------------------------------------------------------------------
    % COMPUTE THRESHOLDS FROM CLEAN DATA
    %-------------------------------------------------------------------------------    
    %fetch ideal tree (no CV or noise)
    tree_clean = trees_all_clean{i};
    %compute thresholds
    thresholds_raw = tree_clean.CutPoint; %threshold with nans
    thresholds = thresholds_raw(~isnan(thresholds_raw));
    thresholds_pan = nan(1,numGenes);
    thresholds_pan(1:length(thresholds)) = thresholds;
    thresholds_all_clean(i,:) = thresholds_pan;
    
    %get predictor order
    genesUsed = tree_clean.CutPredictor;    
    
    %-------------------------------------------------------------------------------
    % PLOT TREE GRAPHS
    %-------------------------------------------------------------------------------
    
    view(tree_clean,'Mode','graph')
    figure;    
    
    %-------------------------------------------------------------------------------
    % PLOT GENE EXPRESSIONS
    %-------------------------------------------------------------------------------
    
    t1 = thresholds_all_clean(i,1);
    
    %------------------------
    % 1D HISTOGRAMS
    if (numGenes == 1)
        histogram(target(:,i), 'Normalization', 'count', 'BinWidth', 0.1);
        title(sprintf('%s : Gene %d with accuracy %g', geneNames{i}, i, balAccuracies(i)));
        if strcmp(geneNames_ranked(i), geneNames(indexOrder(i)))
            sprintf("Error: geneNames_ranked(i) is %s and geneNames(indexOrder(i) is %s.", geneNames_ranked{i}, geneNames{indexOrder(i)});
        end
        hold on;
        histogram(nonTarget(:,i), 'Normalization', 'count', 'BinWidth', 0.1);
        if ~isnan(t1)
            xline(t1, '--r');
        end
        hold off;
        legend({'target', '~target'});
        
    %------------------------    
    % 2D PLOTS
    elseif (numGenes == 2)
        
        %Plot gene expression data points
        hold on;        
        for a = 1:rows
            %get area colours
            colour = rgbconv(structInfo{a,2}{1});
                     
            %average deviation (i.e. estimate of the average absolute noise)
            %noiseLevel == standard deviation
            avgNoise = 0.7979*noiseStDev;
            
            %plot data points, formatting corresponding to area and class 
            if isTarget(a)
                x = genes(a,prevBestGenes);
                y = genes(a,i);
                errorbar(x, y, avgNoise, 'both', '*', 'color', colour)
            else
                x = genes(a,prevBestGenes);
                y = genes(a,i);
                errorbar(x, y, avgNoise, 'both', '.', 'color', colour)
            end            
        end
        
        %axis limits
        xlim([min(genes(:,prevBestGenes)) - 0.1, max(genes(:,prevBestGenes)) + 0.1])
        ylim([min(genes(:,i)) - 0.1, max(genes(:,i)) + 0.1])
        
        %plot thresholds
        if strcmp(genesUsed{1,1}, 'x1')
            xline(t1, '--');
            if ~isnan(thresholds_all_clean(i,2))
                t2 = thresholds_all_clean(i,2);
                %this assumes that if numgenes = 2, prevGeneData is just 1 gene
                line([t1,max(prevGeneData)], [t2,t2], 'LineStyle', '--', 'Color', 'k')                
            end
            
        elseif strcmp (genesUsed{1,1}, 'x2')
            yline(t1, '--');
            if ~isnan(thresholds_all_clean(i,2))
                t2 = thresholds_all_clean(i,2);
                %this assumes that if numgenes = 2, prevGeneData is just 1 gene
                line([t2,t2],[t1, max(prevGeneData)], 'LineStyle', '--', 'Color', 'k')
            end
        end
        hold off;
        title(sprintf('Genes %d and %d with accuracy %g', prevBestGenes, i, balAccuracies(i)));
        %legend({'target', 'threshold', '~target'});
        
    end
end
