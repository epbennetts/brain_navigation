function params = SetTestParams()

params.area = 'Isocortex';
params.numFolds = 10;
params.numNoiseRepeats = 1;
params.costFunction = 'balanced';

params.cols = 19114;
params.samples = 20;
params.area = "Isocortex";
params.prevBestGenes = [];
params.noiseStDev = 1;
params.numNoiseIterations = 15;
params.numFolds = 10;
%plotting
params.numplots = 5;
params.range = 'top';

end