function [] = DT_plot_accuracies(numgenes, top_accuracies, area)
%Accuracies vs numgenes

%Plot accuracy increase:
numgenes_array = 1:numgenes;
figure();
plot(numgenes_array, top_accuracies,'.-b');
title(sprintf('Accuracy vs # genes in %s', area));
xlabel('Num genes used in DTs')
ylabel('Accuracy')
set(gca,'xtick', 0:numgenes)
grid on;

end