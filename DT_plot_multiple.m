close all force;

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
    title(sprintf('Genes %d and %d with accuracy %g', a, b, accuracies(i)));
    hold on;
    plot(nonTarget(:,a),nonTarget(:,b), '.r');
    xline(t1, '--');
    if ~isnan(thresholds_all(i,2))
        yline(t2, '--');
    end
    %but don't we have 2 thresholds
    %better way to do this is with the other functions
    hold off;
    legend({'target', '~target'});
    
end
