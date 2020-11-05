% a=2;
% b=5;
% save('TestData.mat')
% 
% %view workspace: whos

areaNames_redundant = structInfo{:,5};
areaNames = {};
areaIndices = [1];
for i = 1:(size(areaNames_redundant,1)-1)
    if ~strcmp(areaNames_redundant{i}, areaNames_redundant{i+1})
        areaNames = [areaNames; areaNames_redundant{i}];
        areaIndices = [areaIndices; i+1];
    end
end

