clear all;
close all force;

%-------------------------------------------------------------------------------
% Loads in the data, and sets up targets/nontargets for the chosen area:
params = SetParams_AccVsNoise(); 
area = params.area;
genes = filter_nans(area);
[rows, cols] = size(genes);

%-------------------------------------------------------------------------------
% Parameters (!!)
%Check in params file: 
%-------------------------------------------------------------------------------
%size sample subset
%area
%(prev best genes) --> numGenes in DT 
%-------------------------------------------------------------------------------
%"Translate" params (do this better)
%general fixed
costFunction = params.costFunction;
%general variable
sizeSampleSubset = params.sizeSampleSubset;
AccuracyVsNoise_filename = params.AccuracyVsNoise_filename;
numNoiseIterations = params.numNoiseIterations;
numFolds = params.numFolds;
prevBestGenes = params.prevBestGenes;
%just for this script
noiseLevelSamples = params.noiseLevelSamples;
numNoiseSamples = size(noiseLevelSamples,2);
%-------------------------------------------------------------------------------

% Initialize arrays
topAccuracies = NaN(numNoiseSamples,1);

% Loop through different noise levels
for n = 1:numNoiseSamples
    fprintf('(Printed from general) NOISE LEVEL iteration: %d \n', n);
    fprintf('Noise lev: %d \n', noiseLevelSamples(n));
    [indexOrder, geneNames_ranked, balAcc_ranked, confMatrices_ranked, trees_all_clean, balAcc_stDevs_ranked] = ...
                    DT_classification_multiple(sizeSampleSubset, area, prevBestGenes, noiseLevelSamples(n), numFolds, numNoiseIterations);
    topAccuracies(n) = balAcc_ranked(1);
    
    %save, in case terminates:
    save(AccuracyVsNoise_filename)
end

save(AccuracyVsNoise_filename)

WarnWave = [sin(1:.6:400), sin(1:.7:400), sin(1:.4:400)];
Audio = audioplayer(WarnWave, 22050);
play(Audio);