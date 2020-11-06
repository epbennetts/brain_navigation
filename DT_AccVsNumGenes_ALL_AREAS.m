%Finds best genes and accuracies for 1, 2, ... n genes considered in the DT at a
%time. Allows to see how accuracy changes as we increase the number of genes
%in the DT.

clear all;
close all force;

%******************************CHANGE!
%GET PARAMS FROM A SPECIALISED SETPARAMS***************
%FILENAME
%no areas, prevBest Genes , etc.
%


% Loads in the data, and sets up targets/nontargets for the chosen area:
%(see how to make SetTestParams() work without having to do area = params.area, etc).
params = SetParams_AccVsNumGenes_ALL();
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
AccuracyVsNumGenes_ALL_filename = params.AccuracyVsNumGenes_ALL_filename;
AccuracyVsNumGenes_ALL_filename_lighter = params.AccuracyVsNumGenes_ALL_filename_lighter;
AccVsNumGenes_ALL_filename_temp = params.AccuracyVsNumGenes_ALL_filename_temp; 
noiseStDev = params.noiseStDev;
numNoiseIterations = params.numNoiseIterations;
numFolds = params.numFolds;
%just for this script
maxNumGenesInDT = params.maxNumGenesInDT;
stoppingCrit = params.stoppingCrit;
%-------------------------------------------------------------------------------

%get areas
load('AllenGeneDataset_19419_Ben.mat', 'structInfo');
areaNames_doubledUp = structInfo{:,5};
areaNames = getDistinctAreas(areaNames_doubledUp);
numAreas = size(areaNames,1); 

%Initialize overall vars
bestAccs_ALL = NaN(maxNumGenesInDT,numAreas);
bestAccs_stdevs_ALL = NaN(maxNumGenesInDT,numAreas);
bestGeneIndices_ALL = NaN(maxNumGenesInDT,numAreas);
bestGeneNames_ALL = strings(maxNumGenesInDT,numAreas);

for a = 1:numAreas 
    area = areaNames(a);
    [genes, isTarget, classes, geneNames] = filter_nans(area);
    [rows, cols] = size(genes);
    
    fprintf('\n\n-----AREA: %s-----\n', area{1})
    
    %Initialize vars
    prevBestGenes = [];
    
    % Initialize arrays
    bestAccs = NaN(maxNumGenesInDT,1);
    bestAccs_stdevs = NaN(maxNumGenesInDT,1);
    bestGeneIndices = NaN(maxNumGenesInDT,1);
    bestGeneNames = strings(maxNumGenesInDT,1);
    
    geneNames_ranked_all = cell(sizeSampleSubset, maxNumGenesInDT);
    balAcc_ranked_all = zeros(sizeSampleSubset, maxNumGenesInDT);
    confMatrices_ranked_all = NaN(sizeSampleSubset, 4, maxNumGenesInDT);
    balAcc_ranked = NaN(sizeSampleSubset, maxNumGenesInDT);
    
    % Greedy approach: Add 1 more gene into the DT each time
    for n = 1:maxNumGenesInDT
        realNumGenesInDT = length(prevBestGenes) + 1;
        fprintf('(printed from general) Num genes considered at once in the DT is: %d = %d \n', n, realNumGenesInDT)
        
        [indexOrder, geneNames_ranked, balAcc_ranked, confMatrices_ranked, trees_all_clean, balAcc_stDevs_ranked] = ...
            DT_classification_multiple(sizeSampleSubset, area, prevBestGenes, noiseStDev, numFolds, numNoiseIterations);
        bestGeneIndices(n) = indexOrder(1);
        prevBestGenes = [prevBestGenes bestGeneIndices(n)];
        bestGeneNames(n) = geneNames_ranked{1};
        bestAccs(n) = balAcc_ranked(1);
        bestAccs_stdevs(n) = balAcc_stDevs_ranked(1);
        
        geneNames_ranked_all(:,n) = geneNames_ranked;
        %confMatrices_ranked_all(:,:,n) = confMatrices_ranked(sizeSampleSubset, 4, maxNumGenesInDT);
        balAcc_ranked_all(:,n) = balAcc_ranked;
        
        %just in case program doesn't finish:
        save(AccVsNumGenes_ALL_filename_temp);
        
        %stopping criterion: stop if the accuracy decreases either once or twice
        if (stoppingCrit == 1)
            if (n > 1) && (bestAccs(n) < bestAccs(n-1))
                break;
            end
        elseif (stoppingCrit == 2)
            if (n > 2) && (bestAccs(n) < bestAccs(n-1)) && (bestAccs(n-1) < bestAccs(n-2))
                break;
            end
        end
        
    end
    
    bestAccs_ALL(:,a) = bestAccs;
    bestAccs_stdevs_ALL(:,a) = bestAccs_stdevs;
    bestGeneIndices_ALL(:,a) = bestGeneIndices;
    bestGeneNames_ALL(:,a) = bestGeneNames;
    
end

%Main_Results = struct('best_gene_names', bestGeneNames, 'maxNumGenesInDT', maxNumGenesInDT, 'top_accuracies', bestAccs, 'area', area);
save(AccuracyVsNumGenes_ALL_filename)
%save(AccuracyVsNumGenes_filename_lighter, best_gene_names,maxNumGenesInDT, top_accuracies, area)

WarnWave = [sin(1:.6:400), sin(1:.7:400), sin(1:.4:400)];
Audio = audioplayer(WarnWave, 22050);
play(Audio);