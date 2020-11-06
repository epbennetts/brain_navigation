function params = SetParams_AccVsNumGenes_ALL()
%acc vs numgenes for all/multiple areas

%Set vars
cols = 19114; %STATIC
params.cols = cols; %STATIC
%MAIN PARAMS (!!)----------------------------
params.sizeSampleSubset = 10; %!
%%%%%%%%%%params.areas = 'Cerebellar Nuclei';


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
params.stoppingCrit = 2;

%FILE NAME
filename = sprintf('AccuracyVsNumGenes_allAreas_%d_%d.mat',params.sizeSampleSubset,params.maxNumGenesInDT);
params.AccuracyVsNumGenes_AllAreas_filename = filename;

filename_lighter = sprintf('AccuracyVsNumGenes_allAreas_%d_%d_lighter.mat',params.sizeSampleSubset,params.maxNumGenesInDT);
params.AccuracyVsNumGenes_AllAreas_filename_lighter = filename_lighter;

% %PLOTTING
% params.numplots = 5;
% params.range = 'top';

end