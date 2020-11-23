%function[] = DT_AccuracyVsNoiseLevel_plot(noiseLevelSamples, topAccuracies, area)

% Plot accuracy vs noise:
figure();
errorbar(noiseLevelSamples, topAccuracies.*100, topAccs_errors,'.-b');
%plot(noiseLevelSamples, topAccuracies.*100,'.-b');
title(sprintf('Accuracy vs noise level in %s (samples: %d, iters: %d, folds: %d, numGenes: %d)', area, sizeSampleSubset, numNoiseIterations, numFolds, numGenes_temp));
xlabel('Noise level used (\sigma)')
ylabel('Balanced Accuracy (%)')
%set(gca,'xtick', 0:numgenes)
grid on;
