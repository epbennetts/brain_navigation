%TO DO
%later: have clean genes matrix saved somewhere (instead of processing every time)
%generalise to a function that uses n num of genes and just goes through
%everything without plotting
%generalise further to any area


%clear all;
close all force;

cols = 19114;
samples = 10;
area = "Isocortex";
%prev_best_genes = [4307, 636, 148]; %for noise
prev_best_genes = [17];


%function [indexOrder, bal_acc_ranked, geneNames_ranked, thresholds_all] = DT_classification_multiple(samples, area, prev_best_genes)
%Input:
%numgenes: Number of genes to use in each DT
%samples: number of "samples" to complete the algorithm on (can choose samples < cols for
%testing)
%prev_best_genes: array of the indices of the best genes for each
%iteration of increasing numgenes. i.e. for numgenes = 1, don't need
%anything, for numgenes = 2, need indexOrder(1) of the previous iteration,
%for numgenes = 3, need array of size 2 with the 2 prev ones, etc.
%ideally this would eventually be done recursively within the function
%rather than manually outside
%Function:
%The function will train a tree with those genes and test it
%against the real classes. Will also plot the accuracies histogram
%Returns:
%indexOrder: gene index sorted by decreasing accuracy
%accuracies_ranked: Accuracies from highest to lowest
%geneNames_ranked: Gene names ordered by highest to lowest accuracy
%thresholds_all: DT thresholds for every tree

numgenes = length(prev_best_genes) + 1;

[genes, isTarget, geneNames, structInfo] = filter_nans(area);
[rows, cols] = size(genes);
classes = isTarget;

%num of trees we want to plot (not in order)
%plotting = 1:10;

%set vars and empty arrays
bal_accuracies = zeros(samples,1);
f1_score = zeros(samples,1);
precision = zeros(samples,1);
conf_matrices = zeros(samples,4);
thresholds_all = NaN(samples, numgenes);
%genes_used = cell(samples, 2);

%make set of noisy data
%this will be noise of max ~scale st. deviations (because data is z-scored) 
scale = 1;
noise = scale*(rand(size(genes)) - 0.5); 
genes_noisy = genes + noise;

% Previously selected genes as baseline:
prevGeneData = genes(:, prev_best_genes);
prevGeneDataNoisy = genes_noisy(:, prev_best_genes);

%loop, make trees, classify and evaluate
for i = 1:samples
    
    % Set gene data for this iteration:
    gene_combo = [prevGeneData genes(:,i)];
    gene_combo_noisy = [prevGeneDataNoisy genes_noisy(:,i)];
    
    %set cost function
    %matrix order is set by 'ClassNames'
    cost = size(isTarget(isTarget == 0), 1)/size(isTarget(isTarget == 1), 1);
    cost_f = [0, 1; cost, 0];
    
    %train 
    tree_ideal = fitctree(gene_combo, classes, 'MaxNumSplits', numgenes, 'cost', cost_f, 'ClassNames', [1,0]);
    %tree_CV = fitctree(gene_combo, classes, 'MaxNumSplits', numgenes, 'CrossVal','on', 'cost', cost_f, 'ClassNames', [1,0]);
    tree_CV_noise = fitctree(gene_combo_noisy, classes, 'MaxNumSplits', numgenes, 'CrossVal','on', 'cost', cost_f, 'ClassNames', [1,0]);
    
    % Store all thresholds for later:
    %with CV: save 10 thresholds, then average?
    thresholds_raw = tree_ideal.CutPoint; %threshold with NaNs
    thresholds = thresholds_raw(~isnan(thresholds_raw));
    threshold_pan = nan(1,numgenes);
    threshold_pan(1:length(thresholds)) = thresholds;
    thresholds_all(i,:) = threshold_pan;
    genes_used = tree_ideal.CutPredictor;   
 
    %test CV
%     %labels_ideal = predict(tree_ideal, gene_combo); 
%     %labels_noisy = predict(tree_ideal, gene_combo_noisy);   %how well does the ideal tree do with noisy data?
    [labels_CV, score_CV, cost_CV] = kfoldPredict(tree_CV_noise);
    
    %conf matrices -- 
    tp = sum(labels_CV==1 & isTarget==1);
    fn = sum(labels_CV==0 & isTarget==1);
    fp = sum(labels_CV==1 & isTarget==0);
    tn = sum(labels_CV==0 & isTarget==0);    
    conf_matrices(i,:) = [tp, fp, fn, tn];
    
    %calc balanced accuracy
    bal_acc = (tp/(tp+fp) + tn/(tn+fn))/2;
    %if there are no predicted targets this will give nan, so:
    if tp == 0
        %accuracy = (0 + tn/(tn+fn))/2;
        bal_acc = 0;
    end
    bal_accuracies(i) = bal_acc;
    
    %f1 score
    prec = tp/(tp+fp);
    recall = tp/(tp+fn);
    %f1_score(i) = 2*prec*recall/(prec+recall)*100;
    %precision(i) = prec;
    
end

%sort accuracies
bal_accuracies(isnan(bal_accuracies)) = -Inf;
[bal_acc_ranked, indexOrder] = sort(bal_accuracies, 'descend');

%sort F1
% f1_score(isnan(f1_score)) = -Inf;
% [f1_ranked, indexOrder] = sort(f1_score, 'descend');
% conf_matrices_sorted = conf_matrices(indexOrder, :);

% %sort precision
% precision(isnan(precision)) = -Inf;
% [prec_ranked, indexOrder] = sort(precision, 'descend');
% conf_matrices_sorted = conf_matrices(indexOrder, :);


% indexOrder ~ [4,2,3,1]
% X = [5,6,2,1]; X(indexOrder) => [1,6,2,5]

geneNames_ranked = geneNames(indexOrder);

% %all genes ranked
% geneNames_ranked = strings(samples,1);
% for i = 1:size(indexOrder,1)
%     geneIndex = indexOrder(i);
%     genesStr = string(geneNames(geneIndex));
%     geneNames_ranked(i) = genesStr;
% end


%plot the accuracies overall
figure();
histogram(bal_acc_ranked, 'Normalization', 'count', 'NumBins', 10);
title(sprintf('Balanced Accuracy for %g genes in %s', numgenes, area))
xlabel('Balanced accuracy');
ylabel('Counts');


%%***#1 best gene ***
%%doesn't really make sense for cases when many genes have same top accuracy (usually
%%numgenes > 1)
%bestGeneIndex = indexOrder(1);
%bestGene = string(geneNames(bestGeneIndex))


%plotting
range = 'top';
numplots = 5;

%various plotting thresholds
plots = numplots;
middleStart = samples*0.2;
middleEnd = middleStart + plots - 1;
bottomStart = samples - plots + 1;
bottomEnd = samples;

if strcmp(range,'top')
    %top n
    ordered_range = indexOrder(1:plots)';
elseif strcmp(range, 'middle')
    %middle n
    ordered_range = indexOrder(middleStart:middleEnd)';
elseif strcmp(range, 'bottom')
    %bottom n
    ordered_range = indexOrder(bottomStart:bottomEnd)';
end

%make various plots in a certain accuracy range
targetIndices = find(isTarget);
nonTargetIndices = find(~isTarget);
target = genes(isTarget, :);
nonTarget = genes(nonTargetIndices, :);
prevGeneData = genes(:, prev_best_genes);
for i = ordered_range
    
    % Set gene data for this iteration:
    gene_combo = [prevGeneData genes(:,i)];
    
    %set cost function
    %matrix order is set by 'ClassNames'
    cost = size(isTarget(isTarget == 0), 1)/size(isTarget(isTarget == 1), 1);
    cost_f = [0, 1; cost, 0];
    
    %train 
    tree_ideal = fitctree(gene_combo, classes, 'MaxNumSplits', numgenes, 'cost', cost_f, 'ClassNames', [1,0]);
    
    %rename threshold
    t1 = thresholds_all(i,1);
    
    %predictor order
    genes_used = tree_ideal.CutPredictor; 
    
    %view trees, plotsview, etc.
    view(tree_ideal,'Mode','graph')
    figure;    
    if (numgenes == 1)
        histogram(target(:,i), 'Normalization', 'count', 'BinWidth', 0.1);
        title(sprintf('%s : Gene %d with accuracy %g', geneNames{i}, i, bal_accuracies(i)));
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
        
    elseif (numgenes == 2)
        hold on;
        for ro = 1:rows
            colour = rgbconv(structInfo{a,2}{1});
            
%             %do this but with colours
%             x = genes(:,prev_best_genes);
%             x_err = noise(:,prev_best_genes);
%             y = genes(:,i);
%             y_err = noise(:,i);
%             errorbar(x, y, y_err, y_err, x_err, x_err, 'o')
%             
            
            if isTarget(a)                
                x = genes(a,prev_best_genes);
                y = genes(a,i);
                x_err = noise(a,prev_best_genes);
                y_err = noise(a,i);
                errorbar(x, y, y_err, y_err, x_err, x_err, '*', 'color', colour)
            else
                x = genes(a,prev_best_genes);
                y = genes(a,i);
                x_err = noise(a,prev_best_genes);
                y_err = noise(a,i);
                errorbar(x, y, y_err, y_err, x_err, x_err, '.', 'color', colour)
            end
            
            
%             if isTarget(a)
%                 plot(genes(a,prev_best_genes), genes(a,i), '*', 'color', colour);
%             else
%                 plot(genes(a,prev_best_genes), genes(a,i), '.', 'color', colour);
%             end
        end
        xlim([min(genes(:,prev_best_genes)) - 0.1, max(genes(:,prev_best_genes)) + 0.1])
        ylim([min(genes(:,i)) - 0.1, max(genes(:,i)) + 0.1])
        
        %plot thresholds
        if strcmp(genes_used{1,1}, 'x1')
            xline(t1, '--');
            if ~isnan(thresholds_all(i,2))
                t2 = thresholds_all(i,2);
                %this assumes that if numgenes = 2, prevGeneData is just 1 gene
                line([t1,max(prevGeneData)], [t2,t2], 'LineStyle', '--', 'Color', 'k')
                
            end
        elseif strcmp (genes_used{1,1}, 'x2')
            yline(t1, '--');
            if ~isnan(thresholds_all(i,2))
                t2 = thresholds_all(i,2);
                %this assumes that if numgenes = 2, prevGeneData is just 1 gene
                line([t2,t2],[t1, max(prevGeneData)], 'LineStyle', '--', 'Color', 'k')                
            end
        end
        hold off;
        title(sprintf('Genes %d and %d with accuracy %g', prev_best_genes, i, bal_accuracies(i)));
        %legend({'target', 'threshold', '~target'});
       
    end    
    end

