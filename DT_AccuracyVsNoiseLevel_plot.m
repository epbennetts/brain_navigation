%function[] = DT_AccuracyVsNoiseLevel_plot(noiseLevelSamples, topAccuracies, area)

% Plot accuracy vs noise:
figure();
plot(noiseLevelSamples, topAccuracies,'.-b');
title(sprintf('Accuracy vs noise level in %s', area));
xlabel('Noise level used (Std Devs)')
ylabel('Accuracy')
%set(gca,'xtick', 0:numgenes)
grid on;