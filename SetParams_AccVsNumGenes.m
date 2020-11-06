function params = SetParams_AccVsNumGenes()
%Set vars
cols = 19114; %STATIC
params.cols = cols; %STATIC
%MAIN PARAMS (!!)
params.sizeSampleSubset = 10; %!VARIABLE!
params.area = 'Isocortex';


%TREE
params.costFunction = 'balanced';
%CV
params.numFolds = 10;
%noise
params.noiseStDev = 1;
params.numNoiseIterations = 5;

%--------------------------------------------------------------------
%JUST FOR THIS SCRIPT
%--------------------------------------------------------------------
%num genes
params.prevBestGenes = [];
params.maxNumGenesInDT = 10;
%STOPPING CRIT
%i.e. stop adding genes after accuracy decreases either 1ce or 2ce
%0 means no stopping criterion
params.stoppingCrit = 0;

%FILE NAME
filename = sprintf('AccuracyVsNumGenes_%s_%d_%dgenes_%diters.mat',params.area,params.sizeSampleSubset,params.maxNumGenesInDT,params.numNoiseIterations);
params.AccuracyVsNumGenes_filename = filename;

filename_lighter = sprintf('AccuracyVsNumGenes_%s_%d_%dgenes_%diters_lighter.mat',params.area,params.sizeSampleSubset,params.maxNumGenesInDT,params.numNoiseIterations);
params.AccuracyVsNumGenes_filename_lighter = filename_lighter;

% %PLOTTING
% params.numplots = 5;
% params.range = 'top';

end