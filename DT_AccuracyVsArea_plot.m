%-------------------------------------------------------------------------------
% Plot accuracy vs area:
figure();
accuracies = bestAccs_ALL(10,:).*100;
X = categorical(areaNames);
Y = accuracies;

bar(X,Y);
title(sprintf('Accuracy vs area'));
xlabel('Area')
ylabel('Balanced Accuracy (%)')
%set(gca,'xtick', 0:numgenes)
grid on;
hold on;
