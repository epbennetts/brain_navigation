close all force;
clear all;

%target = rand(1000,1)+1;
%nonTarget = rand(100000,1);

%make dummy data
targetCentre = 3;  nonTargetCentre = 0;
targetSize = 100;  nonTargetSize = 10000;
target = normrnd(targetCentre,1,targetSize, 1);
nonTarget = normrnd(nonTargetCentre,1,nonTargetSize, 1);
allData = [target; nonTarget];

%set classes
classes = zeros(size(target,1)+size(nonTarget,1),1);
classes(1:size(target,1),:) = 1;

%set cost function (balanced)
costFunc = ComputeBalancedCostFunc(classes);

%train
tree = fitctree(allData, classes, 'MaxNumSplits', 1, 'cost', costFunc, 'ClassNames', [1,0]);

% Store threshold 
threshold_raw = tree.CutPoint; %threshold with NaNs
threshold = threshold_raw(~isnan(threshold_raw));

%test
predictedLabels = predict(tree, allData);

%calc confusion mat & bal accuracy
[confMat,balAcc] = ComputeConfusion(classes,predictedLabels);



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
