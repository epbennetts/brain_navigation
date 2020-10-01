%Finds best genes and accuracies for 1, 2, ... n genes considered in the DT at a
%time. Allows to see how accuracy changes as we increase the number of genes
%in the DT. 

clear vars;
close all force;

% Loads in the data, and sets up targets/nontargets for the chosen area:
%(see how to make SetTestParams() work without having to do area = params.area, etc).
params = SetParams_AccVsNumGenes();  
area = params.area;
[genes, isTarget, classes, geneNames] = filter_nans(area);
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
AccuracyVsNumGenes_filename = params.AccuracyVsNumGenes_filename;
noiseStDev = params.noiseStDev;
numNoiseIterations = params.numNoiseIterations;
numFolds = params.numFolds;
prevBestGenes = params.prevBestGenes;
%just for this script
maxNumGenesInDT = params.maxNumGenesInDT;
%-------------------------------------------------------------------------------

% Initialize arrays
top_accuracies = NaN(maxNumGenesInDT,1);
best_geneIndices = NaN(maxNumGenesInDT,1);
best_gene_names = strings(maxNumGenesInDT,1);

% Greedy approach: Add 1 more gene into the DT each time 
for n = 1:maxNumGenesInDT
    fprintf('(printed from general) Num genes considered at once in the DT is: %d \n', n)
    
    [indexOrder, geneNames_ranked, balAcc_ranked, confMatrices_ranked, trees_all_clean] = ...
                    DT_classification_multiple(sizeSampleSubset, area, prevBestGenes, noiseStDev, numFolds, numNoiseIterations);
    best_geneIndices(n) = indexOrder(1);
    prevBestGenes = [prevBestGenes best_geneIndices(n)];
    best_gene_names(n) = geneNames_ranked{1};
    top_accuracies(n) = balAcc_ranked(1);
    
    %just in case program doesn't finish:
    save(AccuracyVsNumGenes_filename)
    
    %stopping criterion: the accuracy decreases twice
    if (n > 2) && (top_accuracies(n) < top_accuracies(n-1)) && (top_accuracies(n-1))
        break;
    end
end

save(AccuracyVsNumGenes_filename)