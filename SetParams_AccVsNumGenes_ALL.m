function params = SetParams_AccVsNumGenes_ALL()
%Set vars
cols = 19114; %STATIC
params.cols = cols; %STATIC
%MAIN PARAMS (!!)
params.sizeSampleSubset = cols; %!VARIABLE!

%TREE
params.costFunction = 'balanced';
%CV
params.numFolds = 10;
%noise
params.noiseStDev = 1;
params.numNoiseIterations = 10;

%--------------------------------------------------------------------
%JUST FOR THIS SCRIPT
%--------------------------------------------------------------------
%num genes
params.maxNumGenesInDT = 10;
%STOPPING CRIT
%i.e. stop adding genes after accuracy decreases either 1ce or 2ce
%0 means no stopping criterion
params.stoppingCrit = 0;

%FILE NAME
filename = sprintf('AccVsNumGenes_ALL_%d_%s_%dnoiselev_%dgenes_%diters.mat',params.sizeSampleSubset,date,params.noiseStDev,params.maxNumGenesInDT,params.numNoiseIterations);
params.AccuracyVsNumGenes_ALL_filename = filename;

filename_lighter = sprintf('AccVsNumGenes_ALL_%d_%dgenes_%dnoiselev_%diters_lighter.mat',params.sizeSampleSubset,params.maxNumGenesInDT,params.noiseStDev,params.numNoiseIterations);
params.AccuracyVsNumGenes_ALL_filename_lighter = filename_lighter;

params.AccuracyVsNumGenes_ALL_filename_temp = sprintf('AccVsNumGenes_ALL_%d_%dgenes_%dnoiselev_%diters_temp.mat',params.sizeSampleSubset,params.maxNumGenesInDT,params.noiseStDev,params.numNoiseIterations);
% %PLOTTING
% params.numplots = 5;
% params.range = 'top';

end