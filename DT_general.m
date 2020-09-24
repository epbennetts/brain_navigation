%Finds best genes and accuracies for 1, 2, ... n genes considered in the DT at a
%time. Allows to see how accuracy changes as we increase the number of genes
%in the DT. 

clear vars;
close all force;

% Loads in the data, and sets up targets/nontargets for the chosen area:
%(see how to make SetTestParams() work without having to do area = params.area, etc).
params = SetTestParams();  
area = params.area;
[genes, isTarget, classes, geneNames] = filter_nans(area);
[rows, cols] = size(genes);

%-------------------------------------------------------------------------------
% Parameters (!!)
%-------------------------------------------------------------------------------
prevBestGenes = [];
sizeGeneSubset = cols;
numGenesInDT = 10;
%-------------------------------------------------------------------------------


% More Parameters:
%"Translate" params (do this better)
costFunction = params.costFunction;

noiseStDev = params.noiseStDev;
numNoiseIterations = params.numNoiseIterations;
numFolds = params.numFolds;


% Initialize arrays
top_accuracies = NaN(numGenesInDT,1);
best_geneIndices = NaN(numGenesInDT,1);
best_gene_names = strings(numGenesInDT,1);

% Greedy approach: Add 1 more gene into the DT each time 
for n = 1:numGenesInDT
    disp('num genes considered at once in the DT is:')
    disp(n)
    [indexOrder, geneNames_ranked, balAcc_ranked, confMatrices_ranked, trees_all_clean] = ...
                    DT_classification_multiple(sizeGeneSubset, area, prevBestGenes, noiseStDev, numFolds, numNoiseIterations);
    best_geneIndices(n) = indexOrder(1);
    prevBestGenes = [prevBestGenes best_geneIndices(n)];
    best_gene_names(n) = geneNames_ranked{1};
    top_accuracies(n) = balAcc_ranked(1);
    %stopping criterion: the accuracy decreases
    if (n > 1) && (top_accuracies(n) < top_accuracies(n-1)) 
        break;
    end
end

save('AccuracyVsNumGenes.mat')