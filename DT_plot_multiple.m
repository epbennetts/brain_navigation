
function [] = DT_plot_multiple(genes, geneNames, target, nonTarget, thresholds_all, numplots, range)
%inputs: thresholds, ..., numplots, range
%does: plots the gene expression plots with thresholds for the given 

%various plotting thresholds
plots = numplots;
middleStart = samples*0.2;
middleEnd = middleStart + plots - 1;
bottomStart = samples - plots + 1;
bottomEnd = samples;

if strcmp(range,top)
    %top n
    ordered_range = indexOrder(1:plots)';
elseif strcmp(range, middle)
    %middle n
    ordered_range = indexOrder(middleStart:middleEnd)';
elseif strcmp(range, bottom)
    %bottom n
    ordered_range = indexOrder(bottomStart:bottomEnd)';
end

%make various plots in a certain accuracy range
for i = ordered_range
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
        if strcmp(genenames_ranked(i), geneNames(indexOrder(i)))
            sprintf("Error: genenames_ranked(i) is %s and geneNames(indexOrder(i) is %s.", genenames_ranked(i), geneNames{indexOrder(i)});
        end
        hold on;
        histogram(nonTarget(:,i), 'Normalization', 'count', 'BinWidth', 0.1);
        xline(t1, '--r');
        hold off;
        legend({'target', '~target'});
        
    elseif (numgenes == 2)
        plot(target(:,a),target(:,i), '.b');
        title(sprintf('Genes %d and %d with accuracy %g', a, i, accuracies(i)));
        hold on;
        plot(nonTarget(:,a),nonTarget(:,i), '.r');
        %disp("xline:");
        %disp(i);
        %disp(t1);
        xline(t1, '--');
        if ~isnan(thresholds_all(i,2))
            yline(t2, '--');
        end
        hold off;
        legend({'target', '~target'});
    elseif numgenes == 3
        %unsure
    end
    
end
end
