clear all;
close all force;

%-------------------------------------------------------------------------------
% Loads in the data, and sets up targets/nontargets for the chosen area:
params = SetTestParams(); 
area = params.area;
[genes, isTarget, classes, geneNames] = filter_nans(area);
[rows, cols] = size(genes);

%see how to make SetTestParams() work without having to do area = params.area, etc. 
params = SetTestParams(); 
%-------------------------------------------------------------------------------
% Parameters:
%-------------------------------------------------------------------------------
%"Translate" params (do this better)
costFunction = params.costFunction;

prevBestGenes = params.prevBestGenes;
noiseStDev = params.noiseStDev;
numNoiseIterations = params.numNoiseIterations;
numFolds = params.numFolds;
%-------------------------------------------------------------------------------
% Extra Parameters:
geneSamples = 10;
numgenes = 1;
noiseLevelSamples = [0:50:200];
noiseLevelSamples = [noiseLevelSamples 1000 10000];
numNoiseSamples = size(noiseLevelSamples,2);
%-------------------------------------------------------------------------------

% Initialize arrays
topAccuracy = NaN(numNoiseSamples,1);

% Loop through one gene at a time
for n = 1:numNoiseSamples
    disp(n)
    disp(noiseLevelSamples(n))
    [indexOrder, balAcc_ranked, confMatrices_ranked, geneNames_ranked, trees_all_clean] = ...
                    DT_classification_multiple(geneSamples, area, prevBestGenes, noiseLevelSamples(n), numFolds, numNoiseIterations);
    topAccuracy(n) = balAcc_ranked(1);
end

%-------------------------------------------------------------------------------
% Plot accuracy vs noise:
figure();
plot(noiseLevelSamples, topAccuracy,'.-b');
title(sprintf('Accuracy vs noise level in %s', area));
xlabel('Noise level used (Std Devs)')
ylabel('Accuracy')
%set(gca,'xtick', 0:numgenes)
grid on;

