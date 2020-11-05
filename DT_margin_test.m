close all force;
clear all;

%target = rand(1000,1)+1;
%nonTarget = rand(100000,1);

%make dummy data
targetCentre = 10;  
nonTargetCentre = 0;
targetSize = 1000;  
nonTargetSize = 10000;
target = normrnd(targetCentre,1,targetSize, 1);
nonTarget = normrnd(nonTargetCentre,1,nonTargetSize, 1);
allData = [target; nonTarget];
%for visualising target in histogram:
modifier = 2; %because otherwise looks like there are more targets... Justify or just make 100% random set
x_repetitions = nonTargetSize/targetSize/modifier; %assumes target is smallest
y_repetitions = 1;
target_magnified = repmat(target, x_repetitions, y_repetitions);


%set classes
classes = zeros(size(target,1)+size(nonTarget,1),1);
classes(1:size(target,1),:) = 1;

%set cost function (balanced)
costFunc = ComputeBalancedCostFunc(classes);

%train
tree_unbal = fitctree(allData, classes, 'MaxNumSplits', 1, 'ClassNames', [1,0]);
tree_bal = fitctree(allData, classes, 'MaxNumSplits', 1, 'cost', costFunc, 'ClassNames', [1,0]);

% Store thresholds
%unbalanced
threshold_raw = tree_unbal.CutPoint; %threshold with NaNs
threshold_unbal = threshold_raw(~isnan(threshold_raw));
%balanced
threshold_raw = tree_bal.CutPoint; %threshold with NaNs
threshold_bal = threshold_raw(~isnan(threshold_raw));

%test
predictedLabels = predict(tree_bal, allData);

%calc confusion mat & bal accuracy
[confMat,balAcc] = ComputeConfusion(classes,predictedLabels);


%view
view(tree_unbal,'Mode','graph')
view(tree_bal,'Mode','graph')

%PLOT 1
%gene expr histograms
figure()
hold on;
%title('Gene expression');
histogram(target, 'Normalization', 'count', 'BinWidth', 0.2);
histogram(nonTarget, 'Normalization', 'count', 'BinWidth', 0.1);
%xticks(-5:1:6)
%thresholds
%xline(threshold_unbal, '--r', 'LineWidth', 1);
xline(threshold_bal, '--g', 'LineWidth', 1);
%legend
legend({'target','~target','threshold'}, 'Location','northwest');
xlabel('gene expression')
ylabel('counts')
hold off;

%PLOT 2 -- histograms but with target modified to look like it has similar counts
%(Except it looks as if there are 2x as many targets because we're repeating
%the exact same examples...)
%gene expr histograms
figure()
hold on;
%title('Gene expression');
histogram(target_magnified, 'Normalization', 'count', 'BinWidth', 0.2);
histogram(nonTarget, 'Normalization', 'count', 'BinWidth', 0.1);
%xticks(-5:1:6)
%thresholds
%xline(threshold_unbal, '--r', 'LineWidth', 1);
xline(threshold_bal, '--g', 'LineWidth', 1);
%legend
legend({'target','~target','threshold'}, 'Location','northwest');
xlabel('gene expression')
ylabel('counts')
hold off;
