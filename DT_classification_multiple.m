%TO DO
%later: have clean genes matrix saved somewhere (instead of processing every time)
close all force;
clear vars;

[genes, targetIndices, target, nonTarget, classes, geneNames] = filter_nans;
[rows, cols] = size(genes);

%num of trees we want to check
samples = cols;
%plotting = 1:10;
plotting = [];
losses = zeros(1,samples);
margins = zeros(rows,samples);
%should have multiple thresholds
%needs to be fixed
threshold_old = zeros(samples);

%loop, make trees, classify and evaluate
accuracies_old = zeros(samples,1);
bestGeneIndex_global = 1237;
a = bestGeneIndex_global;
genes_used = cell(samples, 2);
conf_matrices = zeros(samples,4);


thresholds_all = NaN(samples,2);
for i = 1:samples
    
    %set gene pair
    b = i;
    gene_pair = [genes(:,a) genes(:,b)]; %gene = genes(:,i:i+1);
    
    %train
    tree = fitctree(gene_pair, classes, 'MaxNumSplits', 2);
    
    %check first predictor (should be best gene: x1)
    gene_used = tree.CutPredictor;
    genes_used{i,1} = gene_used{1};
    genes_used{i,2} = gene_used{2};
    
    %thresholds
    thresholds_raw = tree.CutPoint; %threshold with NaNs
    thresholds = [];    
    %delete nans
    for j = 1:size(thresholds_raw,1)
        thr = thresholds_raw(j);
        if ~isnan(thr)
            thresholds = [thresholds thr];
            if j>1 
                %disp(thr);
                %disp(thresholds);
            end
        end
    end
    %1st threshold
    t1 =  thresholds(1,1);    
    thresholds_all(i,1) = t1;
    %2nd threshold (doesn't count any extra thresholds)
    if(size(thresholds,2) > 1)
        t2 = thresholds(1,2);
        %disp("t2 = ");
        %disp(t2);
        thresholds_all(i,2) = t2;
        if(size(thresholds,2) > 2)
            sprintf("There were more thresholds");
        end
    end
    %gene_index_used = tree.CutPredictorIndex;
    
    %test
    labels = predict(tree, gene_pair);
    threshold_old = edge(tree, gene_pair, classes);
    
    % calc accuracy
    correct = zeros(rows,1);
    for j = 1:rows
        correct(j) = (strcmp(labels{j}, classes(j)));
    end
    accuracies_old(i) = mean(correct,1);
    
    %calc confusion matrix
    conf_mat = zeros(1,4);
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
    
%     for i = 1:rows
%         %if: target
%         if strcmp(classes(j),"target")
%             
%         end
%     end
%     
%     tp = 
%     tn = 
%     fp = 
%     fn = 
%     
    %if i is in the range of elements we want to plot:
    if ismembertol(i,plotting)
        disp("Plotting from first iteration")
        view(tree,'Mode','graph')
        figure;
        %better to do a scatter plot I think...
        plot(target(:,a),target(:,b), '.b');
        title(sprintf('Genes %d and %d with accuracy %g', a, b, accuracies_old(i)));
        hold on;
        plot(nonTarget(:,a),nonTarget(:,b), '.r');
        xline(t1, '--');
        if(size(thresholds) > 1)
            yline(t2, '--');
        %elseif i == 7
        %    yline(-0.5956, '--');
        end
        %but don't we have 2 thresholds
        %better way to do this is with the other functions
        hold off;
        legend({'target', '~target'});
    end
end

disp(size(accuracies_old(accuracies_old > 0.991),1));
disp(size(accuracies_old(accuracies_old < 0.991),1));

%sort accuracies
[ranked_accuracies, indexOrder] = sort(accuracies_old, 'descend');
numGenes = size(ranked_accuracies, 1);
numGenesArray = 1:numGenes;

%all genes ranked
genenames_ranked = strings(samples,1);
for i = 1:size(indexOrder,1)
    geneIndex = indexOrder(i);
    genesStr = string(geneNames(geneIndex));
    genenames_ranked(i) = genesStr;
end

%plot the accuracies overall
figure;
histogram(ranked_accuracies);
title('Accuracy for 2 genes')
ylabel('accuracy');

% figure;
% plot(numGenesArray, ranked_accuracies, '.');
% title('Accuracy for 2 genes')
% xlabel('num genes');
% ylabel('accuracy');


%various plotting thresholds
plots = 5;
middleStart = samples*0.2;
middleEnd = middleStart + plots - 1;
bottomStart = samples - plots + 1;
bottomEnd = samples;
%top 5
%ordered_range = indexOrder(1:plots)';
%middle 5
%ordered_range = indexOrder(middleStart:middleEnd)';
%bottom 5
ordered_range = indexOrder(bottomStart:bottomEnd)';

%can't do like this because won't iterate
%indices = indexOrder(ordered_range);
for i = ordered_range
    b = i;
    gene_pair = [genes(:,a) genes(:,b)]; %gene = genes(:,i:i+1);
    
    %train
    tree = fitctree(gene_pair, classes, 'MaxNumSplits', 2);
    
    %thresholds
    t1 = thresholds_all(i,1);
    if ~isnan(thresholds_all(i,2))
        t2 = thresholds_all(i,2);
    end
    
    %view trees, plots, etc.
    view(tree,'Mode','graph')
    figure;
    plot(target(:,a),target(:,b), '.b');
    title(sprintf('Genes %d and %d with accuracy %g', a, b, accuracies_old(i)));
    hold on;
    plot(nonTarget(:,a),nonTarget(:,b), '.r');
    %disp("xline:");
    %disp(i);
    %disp(t1);
    xline(t1, '--');
    if ~isnan(thresholds_all(i,2))
        yline(t2, '--');
    end
    %but don't we have 2 thresholds
    %better way to do this is with the other functions
    hold off;
    legend({'target', '~target'});
    
end

%%***#1 best gene ***
%%doesn't really make sense for this one cos so many have same accuracy
%bestGeneIndex = indexOrder(1);
%bestGene = string(geneNames(bestGeneIndex))