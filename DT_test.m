%close all force;
%clear all;

%-------------------------------------------------------------------------------
% Parameters:
area = 'Isocortex';

% Loads in the data, and sets up targets/nontargets for the chosen area:
[genes, isTarget, geneNames, structInfo] = filter_nans(area);
[rows, cols] = size(genes);
classes = isTarget;

%More Parameters
%numgenes = 4;
%samples = cols;

%set prev genes
prev_best_genes = [10];
prevGeneData = genes(:, prev_best_genes);
%gene index I want to examine:
i = 123;
gene_combo = [prevGeneData genes(:,i)];
%-------------------------------------------------------------------------------

%train
tree = fitctree(gene_combo, classes, 'MaxNumSplits', numgenes);

%test
labels = predict(tree, gene_combo);

%view
view(tree,'Mode','graph')

%make plots
%(only works if thresholds are still in workspace)
for a = 1:rows
    hold on;
    if isTarget(a)
        plot(genes(a,prev_best_genes), genes(a,i), '*', 'Color', 'b');
    else
        plot(genes(a,prev_best_genes), genes(a,i), '.',  'Color', 'r');
    end
end
xlim([min(genes(:,prev_best_genes)) - 0.1, max(genes(:,prev_best_genes)) + 0.1])
ylim([min(genes(:,i)) - 0.1, max(genes(:,i)) + 0.1])
x = genes(:,prev_best_genes);
x_err = noise(:,prev_best_genes);
y = genes(:,i);
y_err = noise(:,i);
errorbar(x, y, y_err, y_err, x_err, x_err, 'o')
%disp("xline:");
%disp(i);
%disp(t1);
xline(t1, '--');
if ~isnan(thresholds_all(i,2))
    t2 = thresholds_all(i,2);
    %this assumes that if numgenes = 2, prevGeneData is just 1 gene
    line([min(prevGeneData),t1],[t2,t2], 'LineStyle', '--', 'Color', 'k')
    
end
hold off;
title(sprintf('Genes %d and %d with accuracy %g', prev_best_genes, i, f1_score(i)));


