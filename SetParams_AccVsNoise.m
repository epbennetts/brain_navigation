function params = SetParams_AccVsNoise()
%Set vars

%DATA
params.area = 'Isocortex';
params.cols = 19114;
params.prevBestGenes = [];

%TREE
params.costFunction = 'balanced';
%CV
params.numFolds = 10;
%noise
params.noiseStDev = 1;
params.numNoiseIterations = 5;
params.noiseLevelSamples = [0:0.5:3.5]; %[0:0.2:3.5]

% %PLOTTING
% params.numplots = 5;
% params.range = 'top';

end