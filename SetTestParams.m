function params = SetTestParams()

%DATA
params.area = 'Isocortex';
params.cols = 19114;
params.prevBestGenes = [];%[12782];

%TREE
params.costFunction = 'balanced';
%CV
params.numFolds = 10;
%noise
params.noiseStDev = 0.75;
params.numNoiseIterations = 30;

%PLOTTING
params.numplots = 5;
params.range = 'top';

end