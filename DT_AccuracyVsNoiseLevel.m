clear all;
close all force;

%-------------------------------------------------------------------------------
% Loads in the data, and sets up targets/nontargets for the chosen area:
params = SetTestParams(); 
area = 'Isocortex';
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

AccuracyVsNoise_filename = params.AccuracyVsNoise_filename;

%-------------------------------------------------------------------------------
% Extra Parameters:
sizeSampleSubset = params.sizeSampleSubset;
numGenesInDT = 1;
noiseLevelSamples = params.noiseLevelSamples;
numNoiseSamples = size(noiseLevelSamples,2);
%-------------------------------------------------------------------------------

% Initialize arrays
topAccuracies = NaN(numNoiseSamples,1);

% Loop through one gene at a time
for n = 1:numNoiseSamples
    disp("Current noise iteration:")
    disp(n)
    disp(noiseLevelSamples(n))
    [indexOrder, geneNames_ranked, balAcc_ranked, confMatrices_ranked, trees_all_clean] = ...
                    DT_classification_multiple(sizeSampleSubset, area, prevBestGenes, noiseLevelSamples(n), numFolds, numNoiseIterations);
    topAccuracies(n) = balAcc_ranked(1);
    
    %save, in case terminates:
    save(AccuracyVsNoise_filename)
end

save(AccuracyVsNoise_filename)
