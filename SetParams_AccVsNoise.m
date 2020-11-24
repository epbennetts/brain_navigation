function params = SetParams_AccVsNoise()
%Set vars
cols = 19114; %STATIC
params.cols = cols; %STATIC
%MAIN PARAMS (!!)
params.sizeSampleSubset = cols; %!VARIABLE!
params.area = 'Isocortex';
params.prevBestGenes = [];


%TREE
params.costFunction = 'balanced';
%CV
params.numFolds = 10;
%noise
params.numNoiseIterations = 10;
%num genes
params.maxNumGenesInDT = 1;

%--------------------------------------------------------------------
%JUST FOR THIS SCRIPT
%--------------------------------------------------------------------
i_start = 0;
i_step = 0.25;
i_end = 4;
params.noiseLevelSamples = [i_start:i_step:i_end]; %[0:0.5:5];
numNoisePoints = size(params.noiseLevelSamples,2);
%numNoiseSamples = size(params.noiseLevelSamples,2);

%FILE NAME
today=date;
filename = sprintf('AccuracyVsNoise_%s_%d_%s_%dto%d_%dpoints.mat',params.area,params.sizeSampleSubset,today,i_start, i_end,numNoisePoints);
params.AccuracyVsNoise_filename = filename;

% %PLOTTING
% params.numplots = 5;
% params.range = 'top';

end