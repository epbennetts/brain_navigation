%TO DO
%later: have clean genes matrix saved somewhere (instead of processing every time)
%generalise to a function that uses n num of genes and just goes through
%everything without plotting
%generalise further to any area

%numgenes = 3;
%samples = 100;

function [indexOrder, accuracies_ranked, genenames_ranked, thresholds_all] = DT_classification_multiple(numgenes, samples, area)
%Input: 
%Number of genes to use in each DT and number of
%"samples" to complete the algorithm on (can choose samples < cols for
%testing)
%Function:
%The function will train a tree with those genes and test it
%against the real classes. Will also plot the accuracies histogram
%Returns: 
%indexOrder: gene index sorted by decreasing accuracy
%accuracies_ranked: Accuracies from highest to lowest
%genenames_ranked: Gene names ordered by highest to lowest accuracy
%thresholds_all: DT thresholds for every tree


[genes, targetIndices, classes, geneNames] = filter_nans(area);
[rows, cols] = size(genes);

%num of trees we want to plot (not in order)
%plotting = 1:10;
plotting = [];

%set vars and empty arrays
accuracies = zeros(samples,1);
genes_used = cell(samples, 2);
conf_matrices = zeros(samples,4);
thresholds_all = NaN(samples, numgenes);

if (numgenes > 1)
    bestGeneIndex_global = 1237;
    a = bestGeneIndex_global;
    if (numgenes > 2)
        bestGeneIndex_second = 7;
        b = bestGeneIndex_second;
    end
end


%loop, make trees, classify and evaluate
for i = 1:samples
    
    %set gene pair
    if numgenes == 1
        gene_combo = genes(:,i);
    elseif numgenes == 2
        gene_combo = [genes(:,a) genes(:,i)]; %gene = genes(:,i:i+1);
    elseif numgenes == 3
        gene_combo = [genes(:,a) genes(:,b) genes(:,i)];
    end
    
    %train
    tree = fitctree(gene_combo, classes, 'MaxNumSplits', numgenes);
    
    %check first predictor (should be best gene: x1)
    %doesn't exactly work as doesn't always show second split?
    %gene_used = tree.CutPredictor;
    %genes_used{i,1} = gene_used{1};
    %genes_used{i,2} = gene_used{2};
    
    %thresholds
    thresholds_raw = tree.CutPoint; %threshold with NaNs
    thresholds = [];
    %delete nans
    for j = 1:size(thresholds_raw,1)
        thr = thresholds_raw(j);
        if ~isnan(thr)
            thresholds = [thresholds thr];
            if j>1
            end
        end
    end
    %1st threshold
    if ~isnan(thresholds)
        t1 =  thresholds(1,1);
        thresholds_all(i,1) = t1;
    end
    %2nd & 3rd threshold
    if(size(thresholds,2) == 2)
        t2 = thresholds(1,2);
        thresholds_all(i,2) = t2;
    elseif(size(thresholds,2) == 3)
        t3 = thresholds(1,3);
        thresholds_all(i,3) = t3;
    else
        sprintf("There were more than 3 thresholds");
    end
    %gene_index_used = tree.CutPredictorIndex;
    correct = zeros(rows,1);
    
    %test
    labels = predict(tree, gene_combo);
    %threshold_old = edge(tree, gene_combo, classes);
    
    % calc accuracy
    for j = 1:rows
        correct(j) = (strcmp(labels{j}, classes(j)));
    end
    accuracies(i) = mean(correct,1);
    
    %calc confusion matrix
    %targets
    labels_t = labels(targetIndices == 1);
    tp = size(labels_t(strcmp(labels_t,"target")),1);
    fn = size(labels_t(strcmp(labels_t,"~target")),1);
    %non-targets
    labels_nt = labels(targetIndices == 0);
    fp = size(labels_nt(strcmp(labels_nt,"target")),1);
    tn = size(labels_nt(strcmp(labels_nt,"~target")),1);
    conf_mat = [tp, fp, fn, tn];
    conf_matrices(i,:) = conf_mat;
    
    
    %     %if i is in the range of elements we want to plot:
    %     if ismembertol(i,plotting)
    %         disp("Plotting from first iteration")
    %         view(tree,'Mode','graph')
    %         figure;
    %         %better to do a scatter plot I think...
    %         plot(target(:,a),target(:,b), '.b');
    %         title(sprintf('Genes %d and %d with accuracy %g', a, b, accuracies_old(i)));
    %         hold on;
    %         plot(nonTarget(:,a),nonTarget(:,b), '.r');
    %         xline(t1, '--');
    %         if(size(thresholds) > 1)
    %             yline(t2, '--');
    %         %elseif i == 7
    %         %    yline(-0.5956, '--');
    %         end
    %         %but don't we have 2 thresholds
    %         %better way to do this is with the other functions
    %         hold off;
    %         legend({'target', '~target'});
    %     end
end

%disp(size(accuracies_old(accuracies_old > 0.991),1));
%disp(size(accuracies_old(accuracies_old < 0.991),1));

%sort accuracies
[accuracies_ranked, indexOrder] = sort(accuracies, 'descend');
numGenesTotal = size(accuracies_ranked, 1);
numGenesArray = 1:numGenesTotal;

%all genes ranked
genenames_ranked = strings(samples,1);
for i = 1:size(indexOrder,1)
    geneIndex = indexOrder(i);
    genesStr = string(geneNames(geneIndex));
    genenames_ranked(i) = genesStr;
end



%plot the accuracies overall
figure();
histogram(accuracies_ranked);
title(sprintf('Accuracy for %g genes in %s', numgenes, area))
xlabel('accuracy');
ylabel('counts');




% 
% %various plotting thresholds
% plots = 5;
% middleStart = samples*0.2;
% middleEnd = middleStart + plots - 1;
% bottomStart = samples - plots + 1;
% bottomEnd = samples;
% %top 5
% ordered_range = indexOrder(1:plots)';
% %middle 5
% %ordered_range = indexOrder(middleStart:middleEnd)';
% %bottom 5
% %ordered_range = indexOrder(bottomStart:bottomEnd)';
% 
% %make various plots in a certain accuracy range
% for i = ordered_range
%     %set gene pair
%     if numgenes == 1
%         gene_combo = genes(:,i);
%     elseif numgenes == 2
%         gene_combo = [genes(:,a) genes(:,i)]; %gene = genes(:,i:i+1);
%     elseif numgenes == 3
%         gene_combo = [genes(:,a) genes(:,b) genes(:,i)];
%     end
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
%     if numgenes == 1
%         histogram(target(:,i), 'Normalization', 'count', 'BinWidth', 0.1);
%         title(sprintf('%s : Gene %d with accuracy %g', geneNames{i}, i, accuracies(i)));
%         if strcmp(genenames_ranked(i), geneNames(indexOrder(i)))
%             sprintf("Error: genenames_ranked(i) is %s and geneNames(indexOrder(i) is %s.", genenames_ranked(i), geneNames{indexOrder(i)});
%         end
%         hold on;
%         histogram(nonTarget(:,i), 'Normalization', 'count', 'BinWidth', 0.1);
%         xline(t1, '--r');
%         hold off;
%         legend({'target', '~target'});
%     elseif numgenes == 2
%         plot(target(:,a),target(:,i), '.b');
%         title(sprintf('Genes %d and %d with accuracy %g', a, i, accuracies(i)));
%         hold on;
%         plot(nonTarget(:,a),nonTarget(:,i), '.r');
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

%%***#1 best gene ***
%%doesn't really make sense for this one cos so many have same accuracy
%bestGeneIndex = indexOrder(1);
%bestGene = string(geneNames(bestGeneIndex))

end

