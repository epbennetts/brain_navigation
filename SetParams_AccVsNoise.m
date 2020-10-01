function params = SetParams_AccVsNoise()
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
%noise samples
params.maxNumGenesInDT = 1;
params.noiseLevelSamples = [0:0.5:5]; %[0:0.2:3.5]

%FILE NAME
filename = sprintf('AccuracyVsNoise_%s_%d.mat',params.area,params.sizeSampleSubset);
params.AccuracyVsNoise_filename = filename;



% %PLOTTING
% params.numplots = 5;
% params.range = 'top';

end