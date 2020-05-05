close all force;
clear all;

%target = rand(1000,1)+1;
%nonTarget = rand(100000,1);

%make dummy data
target = normrnd(3,1,100, 1);
nonTarget = normrnd(0,1,10000, 1);
allData = [target; nonTarget];

%set classes
classes = zeros(size(target,1)+size(nonTarget,1),1);
classes(1:size(target,1),:) = 1;

%set costs
cost_f = [0, 100; 1, 0];

%train
tree = fitctree(allData, classes, 'MaxNumSplits', 1, 'cost', cost_f, 'ClassNames', [1,0]);



% Store threshold 
threshold_raw = tree.CutPoint; %threshold with NaNs
threshold = threshold_raw(~isnan(threshold_raw));


%view
view(tree,'Mode','graph')

%plot histograms
histogram(target, 'Normalization', 'count', 'BinWidth', 0.2);
title('Gene expression of dummy targets and non-targets');
hold on;
histogram(nonTarget, 'Normalization', 'count', 'BinWidth', 0.1);
xticks(-5:1:6)
xline(threshold, '--r', 'LineWidth', 1);
hold off;
legend({'target', '~target'});
