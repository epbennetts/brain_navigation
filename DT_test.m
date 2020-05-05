close all force;
clear all;

%-------------------------------------------------------------------------------
% Parameters:
area = 'Isocortex';

% Loads in the data, and sets up targets/nontargets for the chosen area:
[genes, isTarget, classes, geneNames] = filter_nans(area);
[rows, cols] = size(genes);

numgenes = 3;
samples = cols;
%-------------------------------------------------------------------------------

%gene data
prev_best_genes = [1237 7];
prevGeneData = genes(:, prev_best_genes);
gene_combo = [prevGeneData genes(:,48)];

%train
tree = fitctree(gene_combo, classes, 'MaxNumSplits', numgenes);

%test
labels = predict(tree, gene_combo);

%view
view(tree,'Mode','graph')

