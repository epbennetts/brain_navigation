%function[] = DT_AccuracyVsNoiseLevel_plot(noiseLevelSamples, topAccuracies, area)

% Plot accuracy vs noise:
figure();
errorbar(noiseLevelSamples, topAccuracies.*100, topAccs_errors.*100,'.-b');
%plot(noiseLevelSamples, topAccuracies.*100,'.-b');
%title(sprintf('Accuracy vs noise level in %s (samples: %d, iters: %d, folds: %d, numGenes: %d)', area, sizeSampleSubset, numNoiseIterations, numFolds, numGenes_temp));
xlabel('Noise level (\sigma)')
ylabel('Balanced Accuracy (%)')
xlim([-0.25 20.25])
ylim([45,101])
%set(gca,'xtick', 0:numgenes)
grid on;
