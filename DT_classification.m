%TO DO
%put all cols into DT
%rank genes based on largest margin
%later: have clean genes matrix saved somewhere (instead of processing every time)
close all force;
clear vars;

[genes, target, notTarget, classes, geneNames] = filter_nans;
[rows, cols] = size(genes);

%scv = cvpartition(classes,'KFold', 10, 'Stratify',true);

%num of trees we want to check
samples = cols;
plotting = 1:5;
losses = zeros(1,samples);
margins = zeros(rows,samples);
thresholds = zeros(samples);

% tree = fitctree(genes(:,1), classes, 'MaxNumSplits', 1);
% view(tree,'Mode','graph')


accuracies = zeros(samples,1);
for i = 1:samples
    gene = genes(:,i);
    
    %train
    tree = fitctree(gene, classes, 'MaxNumSplits', 1);
    
    %test
    labels = predict(tree,gene);
    
    % calc accuracy
    correct = zeros(rows,1);
    for j = 1:rows
        correct(j) = (strcmp(labels{j}, classes(j)));
    end
    accuracies(i) = mean(correct,1);
    
    %calc margin
    %margin = mean(kfoldMargin(tree));
    %margins(:,i) = margin;
    
    %if i is in the range of elements we want to plot:
        if ismembertol(i,plotting)
            view(tree,'Mode','graph')
            figure;
            histogram(target(:,i), 'Normalization', 'count', 'BinWidth', 0.1);
            title(sprintf('Gene %d with margin ', i));
            hold on;
            histogram(notTarget(:,i), 'Normalization', 'count', 'BinWidth', 0.1);
            hold off;
            legend({'target', '~target'});
        end
end

%fix losses/accuracy first
[sorted_accuracies, indexOrder] = sort(accuracies, 'descend');
numGenes = size(sorted_accuracies, 1);
numGenesArray = 1:numGenes;


%all genes ranked
genes_ranked = strings(samples,1);
for i = 1:size(indexOrder,1)
    geneIndex = indexOrder(i);
    gene = string(geneNames(geneIndex));
    genes_ranked(i) = gene;
end


%view: losses, margins


%plot the accuracies overall
figure;
plot(numGenesArray, sorted_accuracies, '.');
xlabel('num genes');
ylabel('accuracy');

%***#1 best gene *** 
bestGeneIndex = indexOrder(1);
bestGene = string(geneNames(bestGeneIndex))


