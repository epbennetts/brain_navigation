%-------------------------------------------------------------------------------
% Plot accuracy vs area:
figure();
X = categorical(areaNames);
Y = accuracies;
bar(X,Y);
title(sprintf('Accuracy vs area'));
xlabel('Area')
ylabel('Balanced Accuracy (%)')
%set(gca,'xtick', 0:numgenes)
grid on;
hold on;
