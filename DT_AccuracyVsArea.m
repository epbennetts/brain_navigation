clear all;
close all force;

%-------------------------------------------------------------------------------
% Loads in the data, and sets up targets/nontargets for the chosen area:
params = SetParams_AccVsArea(); 
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
AccuracyVsAreas_filename = params.AccuracyVsAreas_filename;
noiseStDev = params.noiseStDev;
numNoiseIterations = params.numNoiseIterations;
numFolds = params.numFolds;
prevBestGenes = params.prevBestGenes;
%just for this script

%-------------------------------------------------------------------------------

%load area info
load('AllenGeneDataset_19419_Ben.mat', '-mat');
areaNames_doubledUp = structInfo{:,5}; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
areaNames_redundant = structInfo{:,5};
areaIndices = [1];
areaNames = {};

for i = 1:(size(areaNames_redundant,1)-1)
    if ~strcmp(areaNames_redundant{i}, areaNames_redundant{i+1})
        areaNames = [areaNames; areaNames_redundant{i}];
        areaIndices = [areaIndices; i+1];
    end
end
numAreas = size(areaNames,1); 
% Initialize arrays
accuracies = NaN(numAreas,1);

% Loop through areas
for n = 1:numAreas
    disp(areaNames(n))
%     if (n>1) && (strcmp(areas_cell{n},areas_cell{n-1}))
%         continue;
%     end
    % Load in the data, and set up targets/nontargets for the chosen area:
    area_curr = areaNames{n};
    [genes, isTarget, classes, geneNames] = filter_nans(area_curr);
    [rows, cols] = size(genes);
    
    %run algorithm on this area
    [indexOrder, geneNames_ranked, balAcc_ranked, confMatrices_ranked, trees_all_clean] = ...
    DT_classification_multiple(sizeSampleSubset, area_curr, prevBestGenes, noiseStDev, numFolds, numNoiseIterations); 
    
    accuracies(n) = balAcc_ranked(1);
    
end

%-------------------------------------------------------------------------------
% Plot accuracy vs area:
figure();
X = categorical(areaNames);
Y = accuracies;
bar(X,Y);
title(sprintf('Accuracy vs area'));
xlabel('Area')
ylabel('Accuracy')
%set(gca,'xtick', 0:numgenes)
grid on;

