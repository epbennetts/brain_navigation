close all force;
clear vars;

%arrays
top_accuracy = [];
best_genes = [];
best_gene_names = []; 

%vars
area = "Isocortex";
[genes, targetIndices, target, nonTarget, classes, geneNames] = filter_nans(area);
[rows, cols] = size(genes);
samples = cols;
numgenes = 3;
for n = 1:numgenes
    [indexOrder, accuracies_ranked, genenames_ranked, thresholds_all] = DT_classification_multiple(n, samples, area, best_genes);
    best_genes = [best_genes indexOrder(1)];
    best_gene_names = [best_gene_names genenames_ranked(1)];
    top_accuracy = [top_accuracy accuracies_ranked(1)];
end
numgenes_array = 1:numgenes;
figure();
plot(numgenes_array, top_accuracy,'.-b');
title(sprintf('Accuracy vs # genes in %s', area));
xlabel('Num genes used in DTs')
ylabel('Accuracy')
set(gca,'xtick',0:numgenes)
grid on;

