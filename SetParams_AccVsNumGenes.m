function params = SetParams_AccVsNumGenes()
%Set vars
cols = 19114; %STATIC
params.cols = cols; %STATIC
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
params.maxNumGenesInDT = 10;
params.prevBestGenes = [];

%FILE NAME
filename = sprintf('AccuracyVsNumGenes_%s_%d.mat',params.area,params.sizeSampleSubset);
params.AccuracyVsNumGenes_filename = filename;

%STOPPING CRIT
%i.e. stop adding genes after accuracy decreases either 1ce or 2ce
params.stoppingCrit = 2;

% %PLOTTING
% params.numplots = 5;
% params.range = 'top';

end