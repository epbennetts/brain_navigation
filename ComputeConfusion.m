function [confusionMatrix,balAcc] = ComputeConfusion(realLabels,predictedLabels)
% Evaluate real labels against predicted labels
%-------------------------------------------------------------------------------

% Construct confusion matrices:
tp = sum(predictedLabels==1 & realLabels==1);
fn = sum(predictedLabels==0 & realLabels==1);
fp = sum(predictedLabels==1 & realLabels==0);
tn = sum(predictedLabels==0 & realLabels==0);

% Put the confusion matrix together
confusionMatrix = [tp, fp, fn, tn];

% Compute balanced accuracy:
if tp == 0 || (tp+fp) == 0
    %accuracy = (0 + tn/(tn+fn))/2;
    balAcc = 0;
else
    balAcc = (tp/(tp+fp) + tn/(tn+fn))/2;
end



end
