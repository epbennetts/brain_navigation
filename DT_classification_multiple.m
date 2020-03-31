%TO DO
%later: have clean genes matrix saved somewhere (instead of processing every time)
close all force;
clear vars;

[genes, target, nonTarget, classes, geneNames] = filter_nans;
[rows, cols] = size(genes);

%num of trees we want to check
samples = 100;
plotting = 41:50;
losses = zeros(1,samples);
margins = zeros(rows,samples);
%should have multiple thresholds
%needs to be fixed
threshold_old = zeros(samples);

%loop, make trees, classify and evaluate
accuracies = zeros(samples,1);
bestGeneIndex_fixed = 1237;
a = bestGeneIndex_fixed;
genes_used = cell(samples, 2);

%edit: 1:samples exc a
for i = 1:samples
    %should we repeat first one as well for comparison??
%     if i == a
%         continue;
%     end
    
    b = i;    
    gene_pair = [genes(:,a) genes(:,b)]; %gene = genes(:,i:i+1);
    
    %train
    tree = fitctree(gene_pair, classes, 'MaxNumSplits', 2);
    
    %check predictor index (should be index of best gene)
    gene_used = tree.CutPredictor;    
    genes_used{i,1} = gene_used{1};
    genes_used{i,2} = gene_used{2};
    %check e.g. gene 7
    if i == 7
        thresholds = tree.CutPoint;
    end
    %gene_index_used = tree.CutPredictorIndex;
    
    %test
    labels = predict(tree,gene_pair);
    threshold_old = edge(tree, gene_pair, classes);
    
    % calc accuracy
    correct = zeros(rows,1);
    for j = 1:rows
        correct(j) = (strcmp(labels{j}, classes(j)));
    end
    accuracies(a) = mean(correct,1);
    
    %if i is in the range of elements we want to plot:
         if ismembertol(i,plotting)
            view(tree,'Mode','graph')
            figure;
            %better to do a scatter plot I think...
            plot(target(:,a),target(:,b), '.b');
            title(sprintf('Genes %d and %d with accuracy %g', a, b, accuracies(a)));
            hold on;
            plot(nonTarget(:,a),nonTarget(:,b), '.r');
            xline(threshold_old, '--');
            yline(threshold_old, '--');
            %but don't we have 2 thresholds
            %better way to do this is with the other functions
            hold off;
            legend({'target', '~target'});
         end
end
 
% %sort accuracies
% [ranked_accuracies, indexOrder] = sort(accuracies, 'descend');
% numGenes = size(ranked_accuracies, 1);
% numGenesArray = 1:numGenes;
% 
% %all genes ranked
% genenames_ranked = strings(samples,1);
% for i = 1:size(indexOrder,1)
%     geneIndex = indexOrder(i);
%     genes = string(geneNames(geneIndex));
%     genenames_ranked(i) = genes;
% end
% 
% %plot the accuracies overall
% figure;
% plot(numGenesArray, ranked_accuracies, '.');
% xlabel('num genes');
% ylabel('accuracy');
% 
% 
% % %testing vars
% % count = [];
% % counter = 0;
% 
% %it's treating i as an entire vector when it's a column vector, not actually looping 
% %various plotting thresholds
% plots = 5;
% middleStart = samples*0.2;
% middleEnd = middleStart + plots;
% bottomStart = samples - plots;
% bottomEnd = samples;
% %top 5
% ordered_range = indexOrder(1:plots)';
% %middle 5
% ordered_range = indexOrder(middleStart:middleEnd)';
% %bottom 5
% ordered_range = indexOrder(bottomStart:bottomEnd)';
% 
% threshold = [];
% for i = ordered_range
% %     counter = counter +1;
% %     count = [count i];
%     
%     genes = genes(:,i);
%     tree = fitctree(genes, classes, 'MaxNumSplits', 1);
%     threshold = edge(tree, genes, classes);
%     threshold = [threshold threshold];
%     
%     %view trees and corresponding distributions
%     view(tree,'Mode','graph')
%     figure;
%     histogram(target(:,i), 'Normalization', 'count', 'BinWidth', 0.1);
%     title(sprintf('%s : Gene %d with accuracy %g', geneNames{i}, i, accuracies(i)));
%     if strcmp(genenames_ranked(i), geneNames(indexOrder(i)))
%        sprintf("Error: genenames_ranked(i) is %s and geneNames(indexOrder(i) is %s.", genenames_ranked(i), geneNames{indexOrder(i)}); 
%     end
%     hold on;
%     histogram(notTarget(:,i), 'Normalization', 'count', 'BinWidth', 0.1);
%     xline(threshold, '--r');
%     hold off;
%     legend({'target', '~target'});
%        
% end
% 
% 
% %***#1 best gene *** 
% bestGeneIndex = indexOrder(1);
% bestGene = string(geneNames(bestGeneIndex))


