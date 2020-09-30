function params = SetParams_AccVsNoise()
%Set vars

params.cols = 19114;
%MAIN PARAMS (!!)
params.sizeSampleSubset = 1000; %!
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
params.noiseLevelSamples = [0:0.5:3.5]; %[0:0.2:3.5]

%FILE NAME
filename = sprintf('AccuracyVsNoise_%s_%d.mat',params.area,params.sizeSampleSubset);
params.AccuracyVsNoise_filename = filename;



% %PLOTTING
% params.numplots = 5;
% params.range = 'top';

end