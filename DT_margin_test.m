close all force;
clear all;

set(0,'DefaultAxesColorOrder','default');
%target = rand(1000,1)+1;
%nonTarget = rand(100000,1);

%make dummy data 1
targetCentre = 8;  
nonTargetCentre = 0;
targetSize = 1000;  
nonTargetSize = 10000;
target1 = normrnd(targetCentre,1,targetSize, 1);
nonTarget1 = normrnd(nonTargetCentre,1,nonTargetSize, 1);
allData1 = [target1; nonTarget1];

%make dummy data 2
sep = 2;
nonTarget2 = nonTarget1-sep;
target2 = target1 + sep;
allData2 = [target2; nonTarget2];

%set classes
classes = zeros(size(target1,1)+size(nonTarget1,1),1);
classes(1:size(target1,1),:) = 1;

%set cost function (balanced)
costFunc = ComputeBalancedCostFunc(classes);

%train
tree1 = fitctree(allData1, classes, 'MaxNumSplits', 1, 'cost', costFunc, 'ClassNames', [1,0]);
tree2 = fitctree(allData2, classes, 'MaxNumSplits', 1, 'cost', costFunc, 'ClassNames', [1,0]);

% Store thresholds
%unbalanced
threshold_raw = tree1.CutPoint; %threshold with NaNs
threshold1 = threshold_raw(~isnan(threshold_raw));
%balanced
threshold_raw = tree2.CutPoint; %threshold with NaNs
threshold2 = threshold_raw(~isnan(threshold_raw));

%test
predictedLabels = predict(tree2, allData2);

%calc confusion mat & bal accuracy
[confMat,balAcc] = ComputeConfusion(classes,predictedLabels);


%view
view(tree1,'Mode','graph')
view(tree2,'Mode','graph')

%PLOT 1
%gene expr histograms
%UNBALANCED
figure()
hold on;
%title('Gene expression');
histogram(target1, 'Normalization', 'count', 'BinWidth', 0.2);
histogram(nonTarget1, 'Normalization', 'count', 'BinWidth', 0.1);
%xticks(-5:1:6)
%thresholds
xline(threshold1, '--r', 'LineWidth', 1);
%xline(threshold_bal, '--g', 'LineWidth', 1);
%legend
legend({'target','~target','threshold'}, 'Location','northeast');
title('Small margin')
xlabel('gene expression')
ylabel('counts')
hold off;

%PLOT 2 
%BALANCED
%-- histograms but with target modified to look like it has similar counts
%(Except it looks as if there are 2x as many targets because we're repeating
%the exact same examples...)
%gene expr histograms
figure()
hold on;
%title('Gene expression');
histogram(target2, 'Normalization', 'count', 'BinWidth', 0.2);
histogram(nonTarget2, 'Normalization', 'count', 'BinWidth', 0.1);
%xticks(-5:1:6)
%thresholds
%xline(threshold_unbal, '--r', 'LineWidth', 1);
xline(threshold2, '--g', 'LineWidth', 1);
%legend
legend({'target','~target','threshold'}, 'Location','northeast');
title('Large margin')
xlabel('gene expression')
ylabel('counts')
hold off;
