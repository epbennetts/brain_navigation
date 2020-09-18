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

noiseStDev = params.noiseStDev;
numNoiseIterations = params.numNoiseIterations;
numFolds = params.numFolds;
%-------------------------------------------------------------------------------
% Extra Parameters:
prevBestGenes = [];
samples = 100;
numgenes = 10;
%-------------------------------------------------------------------------------

% Initialize arrays
top_accuracy = NaN(numgenes,1);
best_genes = NaN(numgenes,1);
best_gene_names = strings(numgenes,1);

% Loop through one gene at a time
for n = 1:numgenes
    disp(n)
    [indexOrder, balAcc_ranked, confMatrices_ranked, geneNames_ranked, trees_all_clean] = ...
                    DT_classification_multiple(samples, area, prevBestGenes, noiseStDev, numFolds, numNoiseIterations);
    best_genes(n) = indexOrder(1);
    prevBestGenes = [prevBestGenes best_genes(n)];
    best_gene_names(n) = geneNames_ranked{1};
    top_accuracy(n) = balAcc_ranked(1);
    if (n > 1) && (top_accuracy(n) < top_accuracy(n-1)) 
        break;
    end
end

%-------------------------------------------------------------------------------
% Plot accuracy increase:
numgenes_array = 1:numgenes;
figure();
plot(numgenes_array, top_accuracy,'.-b');
title(sprintf('Accuracy vs # genes in %s', area));
xlabel('Num genes used in DTs')
ylabel('Accuracy')
set(gca,'xtick', 0:numgenes)
grid on;


%Plot thresholds
numplots = 5;
range = 'top';
%work on this function
%DT_plot_multiple(genes, geneNames, numgenes, best_genes, isTarget, thresholds_all, samples, indexOrder, numplots, range)

%end of program
WarnWave = [sin(1:.6:400), sin(1:.7:400), sin(1:.4:400)];
Audio = audioplayer(WarnWave, 22050);
play(Audio);