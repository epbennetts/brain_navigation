function params = SetTestParams()
%Set vars

%DATA
params.area = 'Isocortex';
params.cols = 19114;
params.prevBestGenes = []; %[12782];

%TREE
params.costFunction = 'balanced';
%CV
params.numFolds = 10;
%noise
params.noiseStDev = 1;
params.numNoiseIterations = 20;

% %PLOTTING
% params.numplots = 5;
% params.range = 'top';

end