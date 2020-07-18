function [costFunc] = ComputeBalancedCostFunc(ClassNames)
% Evaluate balanced cost function
%-------------------------------------------------------------------------------

%calculate num non/targets
numNonTargets = size(ClassNames(ClassNames == 0),1);
numTargets = size(ClassNames(ClassNames == 1),1);

%costs
costTarget = numNonTargets/numTargets;
costNonTarget = 1;

%cost function
%The order (of target, nonTarget) is determined by the order of classes in the
%DT (needs to be consistent: target first)
costFunc = [0, costTarget; costNonTarget, 0];
end

