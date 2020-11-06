%Finds best genes and accuracies for 1, 2, ... n genes considered in the DT at a
%time. Allows to see how accuracy changes as we increase the number of genes
%in the DT.

clear all;
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
AccuracyVsNumGenes_filename_lighter = params.AccuracyVsNumGenes_filename_lighter;
noiseStDev = params.noiseStDev;
numNoiseIterations = params.numNoiseIterations;
numFolds = params.numFolds;
prevBestGenes = params.prevBestGenes;
%just for this script
maxNumGenesInDT = params.maxNumGenesInDT;
stoppingCrit = params.stoppingCrit;
%-------------------------------------------------------------------------------

% Initialize arrays
top_accuracies = NaN(maxNumGenesInDT,1);
top_accs_stdevs = NaN(maxNumGenesInDT,1);
best_geneIndices = NaN(maxNumGenesInDT,1);
best_gene_names = strings(maxNumGenesInDT,1);

geneNames_ranked_all = cell(sizeSampleSubset, maxNumGenesInDT);
balAcc_ranked_all = zeros(sizeSampleSubset, maxNumGenesInDT);
confMatrices_ranked_all = NaN(sizeSampleSubset, 4, maxNumGenesInDT);
balAcc_ranked = NaN(sizeSampleSubset, maxNumGenesInDT);

% Greedy approach: Add 1 more gene into the DT each time
for n = 1:maxNumGenesInDT
    fprintf('(printed from general) Num genes considered at once in the DT is: %d \n', n)
    
    [indexOrder, geneNames_ranked, balAcc_ranked, confMatrices_ranked, trees_all_clean, balAcc_stDevs_ranked] = ...
        DT_classification_multiple(sizeSampleSubset, area, prevBestGenes, noiseStDev, numFolds, numNoiseIterations);
    best_geneIndices(n) = indexOrder(1);
    prevBestGenes = [prevBestGenes best_geneIndices(n)];
    best_gene_names(n) = geneNames_ranked{1};
    top_accuracies(n) = balAcc_ranked(1);
    top_accs_stdevs(n) = balAcc_stDevs_ranked(1);
    
    geneNames_ranked_all(:,n) = geneNames_ranked;
    %confMatrices_ranked_all(:,:,n) = confMatrices_ranked(sizeSampleSubset, 4, maxNumGenesInDT);
    balAcc_ranked_all(:,n) = balAcc_ranked;
    
    %just in case program doesn't finish:
    AccuracyVsNumGenes_filename_temp = sprintf('AccuracyVsNumGenes_%s_%d_%d_temp.mat',params.area,params.sizeSampleSubset,params.maxNumGenesInDT);
    save(AccuracyVsNumGenes_filename_temp);
    
    %stopping criterion: stop if the accuracy decreases either once or twice
    if (stoppingCrit == 1)
        if (n > 1) && (top_accuracies(n) < top_accuracies(n-1))
            break;
        end
    elseif (stoppingCrit == 2)
        if (n > 2) && (top_accuracies(n) < top_accuracies(n-1)) && (top_accuracies(n-1) < top_accuracies(n-2))
            break;
        end
    end
    
end

Main_Results = struct('best_gene_names', best_gene_names, 'maxNumGenesInDT', maxNumGenesInDT, 'top_accuracies', top_accuracies, 'area', area);
save(AccuracyVsNumGenes_filename)
%save(AccuracyVsNumGenes_filename_lighter, best_gene_names,maxNumGenesInDT, top_accuracies, area)