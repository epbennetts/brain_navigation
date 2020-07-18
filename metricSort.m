function [metric_ranked, indexOrder] = metricSort(metric, direction)
% Sorts our desired accuracy metric (or other), putting the nans at the end
%
%string *direction* can be 'ascend' or 'descend'
%Usage: [f1_ranked, indexOrder] = metricSort(f1_score, 'descend')
%-------------------------------------------------------------------------------

%put nans at end in descending case
if strcmp(direction, 'descend')
    metric(isnan(metric)) = -Inf;
end

%sort accuracies
[metric_ranked, indexOrder] = sort(metric, direction);

end