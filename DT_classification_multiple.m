%TO DO
%later: have clean genes matrix saved somewhere (instead of processing every time)
%generalise to a function that uses n num of genes and just goes through
%everything without plotting
%generalise further to any area


% clear all;
% close all force;
% 
% samples = 19114;
% area = "Isocortex";
% prev_best_genes = [];%[1237];%, 7];%, 48];


function [indexOrder, accuracies_ranked, geneNames_ranked, thresholds_all] = DT_classification_multiple(samples, area, prev_best_genes)
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

% COULD BE EASIER?:
% isTarget
% targetIndices = find(isTarget);
% notTargetIndices = find(~isTarget);

[genes, targetIndices, target, nonTarget, classes, geneNames] = filter_nans(area);
[rows, cols] = size(genes);

%num of trees we want to plot (not in order)
%plotting = 1:10;
%plotting = [];

%set vars and empty arrays
accuracies = zeros(samples,1);
f1_score = zeros(samples,1);
conf_matrices = zeros(samples,4);
thresholds_all = NaN(samples, numgenes);
%genes_used = cell(samples, 2);

% Previously selected genes as baseline:
prevGeneData = genes(:, prev_best_genes);

%loop, make trees, classify and evaluate
for i = 1:samples
    
    % Set gene data for this iteration:
    gene_combo = [prevGeneData genes(:,i)];
    
    %train
    tree = fitctree(gene_combo, classes, 'MaxNumSplits', numgenes);
    
    %check first predictor (should be best gene: x1)
    %doesn't exactly work as doesn't always show second split?
    %gene_used = tree.CutPredictor;
    %genes_used{i,1} = gene_used{1};
    %genes_used{i,2} = gene_used{2};
    
    % Store all thresholds for later:
    thresholds_raw = tree.CutPoint; %threshold with NaNs
    thresholds = thresholds_raw(~isnan(thresholds_raw));
    threshold_pan = nan(1,numgenes);
    threshold_pan(1:length(thresholds)) = thresholds;
    thresholds_all(i,:) = threshold_pan;
    %gene_index_used = tree.CutPredictorIndex;
    
    %test
    labels = predict(tree, gene_combo);
    %threshold_old = edge(tree, gene_combo, classes);
    
    % calc old accuracy (for reference)
    %        correct = zeros(rows,1);
    %         for j = 1:rows
    %             correct(j) = (strcmp(labels{j}, classes(j)));
    %         end
    %         accuracies_old(i) = mean(correct,1);
    
    %calc confusion matrix
    %targets
    
    % MAYBE CLEARER WITH BINARY FORMULATION:
    % tp = sum(labels==1 & isTarget==1);
    % fn = sum(labels==1 & isTarget==0);
    
    % labels_t is the set of predictions for areas that are targets
    labels_t = labels(targetIndices == 1);
    tp = sum(strcmp(labels_t,'target')); % much faster if binary labels
    fn = sum(strcmp(labels_t,'~target'));
    %non-targets
    labels_nt = labels(targetIndices == 0);
    fp = sum(strcmp(labels_nt,'target'));
    tn = sum(strcmp(labels_nt,'~target'));
    conf_matrices(i,:) = [tp, fp, fn, tn];
    
    %calc balanced accuracy
    accuracies(i) = (tp/(tp+fp) + tn/(tn+fn))/2;
    prec = tp/(tp+fp);
    recall = tp/(tp+fn);
    f1_score(i) = 2*prec*recall/(prec+recall);
    %if there are no predicted targets this will give nan, so:
    % if isnan(accuracy)
    %     accuracy = (0 + tn/(tn+fn))/2;
    % end
end


%sort accuracies
accuracies(isnan(accuracies)) = -Inf;
[accuracies_ranked, indexOrder] = sort(accuracies, 'descend');

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
histogram(accuracies_ranked);
title(sprintf('Accuracy for %g genes in %s', numgenes, area))
xlabel('accuracy');
ylabel('counts');


%%***#1 best gene ***
%%doesn't really make sense for cases when many genes have same top accuracy (usually
%%numgenes > 1)
%bestGeneIndex = indexOrder(1);
%bestGene = string(geneNames(bestGeneIndex))


% 
% %plotting
% range = 'top';
% numplots = 5;
% 
% 
% %various plotting thresholds
% plots = numplots;
% middleStart = samples*0.2;
% middleEnd = middleStart + plots - 1;
% bottomStart = samples - plots + 1;
% bottomEnd = samples;
% 
% if strcmp(range,'top')
%     %top n
%     ordered_range = indexOrder(1:plots)';
% elseif strcmp(range, 'middle')
%     %middle n
%     ordered_range = indexOrder(middleStart:middleEnd)';
% elseif strcmp(range, 'bottom')
%     %bottom n
%     ordered_range = indexOrder(bottomStart:bottomEnd)';
% end
% 
% %make various plots in a certain accuracy range
% prevGeneData = genes(:, prev_best_genes);
% for i = ordered_range
%     
%     % Set gene data for this iteration:
%     gene_combo = [prevGeneData genes(:,i)];
%     
%     %train
%     tree = fitctree(gene_combo, classes, 'MaxNumSplits', numgenes);
%     
%     %thresholds
%     t1 = thresholds_all(i,1);
%     if(numgenes > 1 && ~isnan(thresholds_all(i,2)))
%         t2 = thresholds_all(i,2);
%         if(numgenes > 2 && ~isnan(thresholds_all(i,3)))
%             t3 = thresholds_all(i,3);
%         end
%     end
%     
%     %view trees, plots, etc.
%     view(tree,'Mode','graph')
%     figure;    
%     if (numgenes == 1)
%         histogram(target(:,i), 'Normalization', 'count', 'BinWidth', 0.1);
%         title(sprintf('%s : Gene %d with accuracy %g', geneNames{i}, i, accuracies(i)));
%         if strcmp(geneNames_ranked(i), geneNames(indexOrder(i)))
%             sprintf("Error: geneNames_ranked(i) is %s and geneNames(indexOrder(i) is %s.", geneNames_ranked{i}, geneNames{indexOrder(i)});
%         end
%         hold on;
%         histogram(nonTarget(:,i), 'Normalization', 'count', 'BinWidth', 0.1);
%         if ~isnan(thresholds_all(i,1))
%             xline(t1, '--r');
%         end        
%         hold off;
%         legend({'target', '~target'});
%         
%     elseif (numgenes == 2)
%         plot(target(:,prev_best_genes),target(:,i), '.b');
%         title(sprintf('Genes %d and %d with accuracy %g', prev_best_genes, i, accuracies(i)));
%         hold on;
%         plot(nonTarget(:,prev_best_genes), nonTarget(:,i), '.r');
%         %disp("xline:");
%         %disp(i);
%         %disp(t1);
%         xline(t1, '--');
%         if ~isnan(thresholds_all(i,2))
%             yline(t2, '--');
%         end
%         hold off;
%         legend({'target', '~target'});
%     elseif numgenes == 3
%         %unsure
%     end
%     
% end

