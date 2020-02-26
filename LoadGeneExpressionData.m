function [geneData, geneInfo, structInfo] = LoadGeneExpressionData()
    % Use: [geneData, geneInfo, structInfo] = LoadGeneExpressionData()
    % Loads gene data from a file 

    dataFile = 'AllenGeneDataset_19419.mat';
    load(dataFile, 'GeneExpData', 'geneInfo', 'structInfo');

    % Use combination (z-score) sections to estimate energy-based expression:
    geneData = GeneExpData.combZ.('energy');

end