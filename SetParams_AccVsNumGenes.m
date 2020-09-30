function params = SetParams_AccVsNumGenes()
%Set vars

params.cols = 19114;
%MAIN PARAMS (!!)
params.sizeSampleSubset = cols; %!
params.area = 'Isocortex';
params.prevBestGenes = []; 


%TREE
params.costFunction = 'balanced';
%CV
params.numFolds = 10;
%noise
params.noiseStDev = 1;
params.numNoiseIterations = 5;
%num genes
params.numGenesInDT = 10;

%FILE NAME
filename = sprintf('AccuracyVsNumGenes_%s_%d.mat',params.area,params.sizeSampleSubset);
params.AccuracyVsNumGenes_filename = filename;



% %PLOTTING
% params.numplots = 5;
% params.range = 'top';

end