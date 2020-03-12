%TO DO
%put all cols into DT
%rank genes based on largest margin
%later: have clean genes matrix saved somewhere (instead of processing every time)
close all force;
clear vars;

[genes, target, notTarget, classes, geneNames] = filter_nans;
[rows, cols] = size(genes);

scv = cvpartition(classes,'KFold', 10, 'Stratify',true);

%num of trees we want to check
samples = 50;
%plotting = 10;
plotting = 1:5;
losses = zeros(1,samples);
margins = zeros(rows,samples);
thresholds = zeros(samples);

for i = 1:samples
    gene = genes(:,i);
    tree = fitctree(gene, classes, 'Crossval', 'on', 'CVPartition', scv, 'MaxNumSplits', 1);
    L = kfoldLoss(tree);
    %kfoldMargin seems to create an entire column of margins... Each
    %element's value depending on what fold it was in
    margin = mean(kfoldMargin(tree));
    losses(i) = L;
    margins(:,i) = margin;
    %if i is in the range of elements we want to plot:
    if ismembertol(i,plotting)
        view(tree.Trained{1},'Mode','graph')
        figure;
        histogram(target(:,i), 'Normalization', 'count', 'BinWidth', 0.1);
        title(sprintf('Gene %d with margin %g', i, margin));
        hold on;
        histogram(notTarget(:,i), 'Normalization', 'count', 'BinWidth', 0.1);
        hold off;
        legend({'target', '~target'});
    end
end

[sorted_losses, indexOrder] = sort(losses, 'ascend');
numGenes = size(sorted_losses, 2);
numGenesArray = 1:numGenes;

% plot(numGenesArray, sorted_losses, '.');
% xlabel('num genes');
% ylabel('kfoldLoss');

%***#1 best gene *** (for now)
bestGeneIndex = indexOrder(1);
bestGene = string(geneNames(bestGeneIndex))


%this is wrong
% for i = indexOrder(1:10)
%     view(tree.Trained{1},'Mode','graph')
%     figure;
%     histogram(target(:,i), 'Normalization', 'count', 'BinWidth', 0.1);
%     title(sprintf('Gene %d with threshold ...', i));
%     hold on;
%     histogram(notTarget(:,i), 'Normalization', 'count', 'BinWidth', 0.1);
%     hold off;
%     legend({'target','~target'});
% end


%all genes ranked
genes_ranked = strings(1,samples);
for i = 1:size(indexOrder,2)
    geneIndex = indexOrder(i);
    gene = string(geneNames(geneIndex));
    genes_ranked(i) = gene;
end

%view: losses, margins

