
function [predictedLabels] = kFoldPredictNoisy(geneCombo, noiseComboIter, numAreas, classes, numgenes, noiseStDev, costFunc, partition, numFolds)
%Classifies the gene data using k-fold CV, with noise in the test set of each fold
%Returns: An array of predicted labels/classes.
%-------------------------------------------------------------------------------

%iterate through the k folds
predictedLabels = nan(numAreas,1);
for k = 1:numFolds
    %fprintf('CV fold: %d \n', k);
    
    % SET TRAINING AND TESTING VARS FOR THIS FOLD:
    %indices
    trainIndices = training(partition, k);
    testIndices = test(partition, k);
    %gene data
    geneDataTrain = geneCombo(trainIndices,:);
    geneDataTest = geneCombo(testIndices,:);
    %noise
    geneNoiseTest = noiseComboIter(testIndices,:);
    %class labels
    classLabelsTrain = classes(trainIndices);
    %classLabelsTest = classes(testIndices);
        
    
    %TRAINING
    % Get the training bits and train the model
    treeFoldTrain = fitctree(geneDataTrain, classLabelsTrain, 'MaxNumSplits', numgenes, 'cost', costFunc, 'ClassNames', [1,0]);
    
    %TESTING
    % Add noise to test gene-expression data
    %geneDataTestNoisy = geneDataTest + noiseStDev*(randn(size(geneDataTest)));
    geneDataTestNoisy = geneDataTest + geneNoiseTest; %load from outer func
    % Store predictions from the trained model on the noisy test data
    predictedLabels(testIndices) = predict(treeFoldTrain, geneDataTestNoisy);
    
end