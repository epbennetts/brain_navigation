%function[] = DT_AccuracyVsNoiseLevel_plot(noiseLevelSamples, topAccuracies, area)

% Plot accuracy vs noise:
figure();
plot(noiseLevelSamples, topAccuracies.*100,'.-b');
title(sprintf('Accuracy vs noise level in %s', area));
xlabel('Noise level used (\sigma)')
ylabel('Balanced Accuracy (%)')
%set(gca,'xtick', 0:numgenes)
grid on;