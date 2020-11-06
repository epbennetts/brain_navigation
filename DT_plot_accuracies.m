function [] = DT_plot_accuracies(maxNumGenesInDT, top_accuracies, top_accs_stdevs, area)
%Accuracies vs numgenes

%Plot accuracy increase:
numgenes_array = 1:maxNumGenesInDT;
figure();
errorbar(numgenes_array, top_accuracies, top_accs_stdevs,'.-b');
title(sprintf('Balanced Accuracy vs genes in %s', area));
xlabel('Number of genes used in Decision Tree')
ylabel('Balanced Accuracy (%)')
set(gca,'xtick', 0:maxNumGenesInDT)
xlim([0.5 10.5])
ylim([0 1])
grid on;
hold on;

end