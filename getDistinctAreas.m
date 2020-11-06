function areaNames = getDistinctAreas(areaNames_doubledUp)
%get area names without double-ups

%areaIndices = [1];
areaNames = {};
areaNames = [areaNames; areaNames_doubledUp{1}];

for i = 1:(size(areaNames_doubledUp,1)-1)
    if ~strcmp(areaNames_doubledUp{i}, areaNames_doubledUp{i+1})
        areaNames = [areaNames; areaNames_doubledUp{i+1}];
        %areaIndices = [areaIndices; i+1];
    end
end