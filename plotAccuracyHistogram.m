function [] = plotAccuracyHistogram(balAccuracies_ranked, numBins, numGenesInDT, area, sizeGeneSubset)
%plots the histogram of accuracies for n genes considered in the DT in a
%particular brain area

%ACCURACY HISTOGRAMs
figure();
histogram(balAccuracies_ranked, 'Normalization', 'count', 'NumBins', numBins);
title(sprintf('Balanced Accuracy for %g genes in %s (samples = %g)', numGenesInDT, area, sizeGeneSubset))
xlabel('Balanced accuracy');
ylabel('Counts');