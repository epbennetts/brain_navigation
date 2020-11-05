function [] = DT_plot_accuracies(maxNumGenesInDT, top_accuracies, area)
%Accuracies vs numgenes

%Plot accuracy increase:
numgenes_array = 1:maxNumGenesInDT;
figure();
plot(numgenes_array, top_accuracies,'.-b');
title(sprintf('Balanced Accuracy vs genes in %s', area));
xlabel('Number of genes used in Decision Tree')
ylabel('Balanced Accuracy (%)')
set(gca,'xtick', 0:maxNumGenesInDT)
grid on;
hold on;

end