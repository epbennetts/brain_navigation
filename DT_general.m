close all force;
clear vars;

%-------------------------------------------------------------------------------
% Parameters:
numgenes = 3;
area = 'Isocortex';
%-------------------------------------------------------------------------------

% Loads in th data, and sets up targets/nontargets for the chosen area:
[genes, targetIndices, target, nonTarget, classes, geneNames] = filter_nans(area);
[rows, cols] = size(genes);
samples = cols;

% Initialize arrays
top_accuracy = zeros(numgenes,1);
best_genes = zeros(numgenes,1);
best_gene_names = zeros(numgenes,1);

% Loop through one gene at a time
for n = 1:numgenes
    [indexOrder, accuracies_ranked, genenames_ranked, thresholds_all] = ...
                    DT_classification_multiple(n, samples, area, best_genes(1:n-1));
    best_genes(n) = indexOrder(1);
    best_gene_names(n) = genenames_ranked(1);
    top_accuracy(n) = accuracies_ranked(1);
end

%-------------------------------------------------------------------------------
% Plot:
numgenes_array = 1:numgenes;
figure();
plot(numgenes_array, top_accuracy,'.-b');
title(sprintf('Accuracy vs # genes in %s', area));
xlabel('Num genes used in DTs')
ylabel('Accuracy')
set(gca,'xtick',0:numgenes)
grid on;
