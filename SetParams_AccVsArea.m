function params = SetParams_AccVsArea()
%Set vars
cols = 19114; %STATIC
params.cols = cols; %STATIC
%MAIN PARAMS (!!)
params.sizeSampleSubset = cols; %!
params.prevBestGenes = []; 


%TREE
params.costFunction = 'balanced';
%CV
params.numFolds = 10;
%noise
params.noiseStDev = 1;
params.numNoiseIterations = 5;
%num genes
params.maxNumGenesInDT = 1;

%FILE NAME
filename = sprintf('AccuracyVsAreas_%d_%d.mat',params.sizeSampleSubset,params.maxNumGenesInDT);
params.AccuracyVsAreas_filename = filename;

%--------------------------------------------------------------------
%JUST FOR THIS SCRIPT
%--------------------------------------------------------------------


% %PLOTTING
% params.numplots = 5;
% params.range = 'top';

end