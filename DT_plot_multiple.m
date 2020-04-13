close all force;
numgenes = 2;

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

for i = ordered_range
    %set gene pair
    b = i;
    c = i+1;
    if numgenes == 1
        gene_combo = genes(:,i);;
    elseif numgenes == 2
        gene_combo = [genes(:,a) genes(:,b)]; %gene = genes(:,i:i+1);
    elseif numgenes == 3
        gene_combo = [genes(:,a) genes(:,b) genes(:,c)]; 
    end    
    
    %train
    tree = fitctree(gene_pair, classes, 'MaxNumSplits', numgenes);
    
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
    if numgenes == 1
        %blah
    elseif numgenes == 2
        plot(target(:,a),target(:,b), '.b');
        title(sprintf('Genes %d and %d with accuracy %g', a, b, accuracies(i)));
        xlabel('gene 1 expression')
        ylabel('gene 2 expression')
        hold on;
        plot(nonTarget(:,a),nonTarget(:,b), '.r');
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
