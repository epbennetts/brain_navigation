
function [] = DT_plot_multiple(genes, geneNames, target, nonTarget, thresholds_all, samples, indexOrder, numplots, range)
%inputs: thresholds, ..., numplots, range
%does: plots the gene expression plots with thresholds for the given 

%-------------------------------------------------------------------------------
% ACCURACY HISTOGRAMS
%-------------------------------------------------------------------------------
figure();
histogram(balAcc_ranked, 'Normalization', 'count', 'NumBins', 10);
title(sprintf('Balanced Accuracy for %g genes in %s (samples = %g)', numGenes, area, samples))
xlabel('Balanced accuracy');
ylabel('Counts');

%-------------------------------------------------------------------------------
% GENE EXPRESSION PLOTS
%-------------------------------------------------------------------------------
% %Don't think need this:
% %plotting vars
% prevGeneData = genes(:, prevBestGenes);
% thresholds_all_clean = nan(samples, numGenes);
% %non/target indices (to plot them differently)
% targetIndices = isTarget;
% nonTargetIndices = find(~isTarget);
% target = genes(targetIndices, :);
% nonTarget = genes(nonTargetIndices, :);

%various plotting thresholds
plots = numplots;
middleStart = samples*0.2;
middleEnd = middleStart + plots - 1;
bottomStart = samples - plots + 1;
bottomEnd = samples;

if strcmp(range,'top')
    %top n
    ordered_range = indexOrder(1:plots)';
elseif strcmp(range, 'middle')
    %middle n
    ordered_range = indexOrder(middleStart:middleEnd)';
elseif strcmp(range, 'bottom')
    %bottom n
    ordered_range = indexOrder(bottomStart:bottomEnd)';
end

%make various plots in a certain accuracy range
prevGeneData = genes(:, prev_best_genes);
for i = ordered_range
    
    % Set gene data for this iteration:
    gene_combo = [prevGeneData genes(:,i)];
    
    %train
    tree = fitctree(gene_combo, classes, 'MaxNumSplits', numgenes);
    
    %thresholds
    t1 = thresholds_all(i,1);
    if(numgenes > 1 && ~isnan(thresholds_all(i,2)))
        t2 = thresholds_all(i,2);
        if(numgenes > 2 && ~isnan(thresholds_all(i,3)))
            t3 = thresholds_all(i,3);
        end
    end
    
    %view trees, plots, etc.
    view(tree,'Mode','graph')
    figure;    
    if (numgenes == 1)
        histogram(target(:,i), 'Normalization', 'count', 'BinWidth', 0.1);
        title(sprintf('%s : Gene %d with accuracy %g', geneNames{i}, i, accuracies(i)));
        if strcmp(geneNames_ranked(i), geneNames(indexOrder(i)))
            sprintf("Error: geneNames_ranked(i) is %s and geneNames(indexOrder(i) is %s.", geneNames_ranked{i}, geneNames{indexOrder(i)});
        end
        hold on;
        histogram(nonTarget(:,i), 'Normalization', 'count', 'BinWidth', 0.1);
        if ~isnan(thresholds_all(i,1))
            xline(t1, '--r');
        end        
        hold off;
        legend({'target', '~target'});
        
    elseif (numgenes == 2)
        %for a = 1:rows
            hold on;
            plot(target(:,prev_best_genes),target(:,i), '.b');
            title(sprintf('Genes %d and %d with accuracy %g', prev_best_genes, i, accuracies(i)));
            plot(nonTarget(:,prev_best_genes), nonTarget(:,i), '.r');
            %disp("xline:");
            %disp(i);
            %disp(t1);
            xline(t1, '--');
            if ~isnan(thresholds_all(i,2))
                yline(t2, '--');
            end
            legend({'target', '~target'});
        %end
    end    
end