close all force;
clear vars;
[rows, cols] = size(genes);

%vars
numgenes = 3;
samples = cols;
area = "Isocortex";

top_accuracy = [];
for n = 1:numgenes
    [indexOrder, ranked_accuracies, genenames_ranked] = DT_classification_multiple(n, samples, area);
    acc = ranked_accuracies(1);
    top_accuracy = [top_accuracy acc];
end

figure();
plot(1:numgenes, top_accuracy,'.-b');
title(sprintf('Accuracy vs # genes in %s', area));
xlabel('Num genes used in DTs')
ylabel('Accuracy')
set(gca,'xtick',0:numgenes)
grid on;
