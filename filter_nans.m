function [genes, isTarget, geneNames, rowInfo] = filter_nans(area)

%load new dataset (have to check which section)
load('AllenGeneDataset_19419_Ben.mat');
genes = GeneExpData.combZ.energy; %for example
colInfo = geneInfo;
rowInfo = structInfo;

%calc how many nans in rows and cols
propNansRows = mean(isnan(genes),2);
propNansCols = mean(isnan(genes),1);

%histogram of what comes out of each of these commands.
% figure(1);
% histogram(propNansRows, 'Normalization', 'probability');
% title('Proportion of nans in rows');
% figure(2);
% histogram(propNansCols,'Normalization', 'probability', 'BinLimits',[0,0.3]);
% title('Proportion of nans in columns')

%proportion of rows or cols that are allowed to be nans: 
max_nans_rows = 0.1;
max_nans_cols = 0.1;

%cols
%this returns an array of size cols with 0/1 values
doRemove = (propNansCols > max_nans_cols);
%remove the arrays at indices with 1
genes(:,doRemove) = [];
colInfo(doRemove,:) = [];

%rows
doRemove = (propNansRows > max_nans_rows);
genes(doRemove,:) = [];
rowInfo(doRemove,:) = [];
areas = table2array(rowInfo(:,5));

[rows, cols] = size(genes);

%should probably also have a min amount of vals in target area (e.g.
%isocortex) --> later

%calc average for target and ~target elements
isTarget = (strcmp(area, areas));
target = genes(isTarget, :);
nonTarget = genes(~isTarget, :);

%set remaining nans to average of their area (target or not target)
avgTarget = mean(~isnan(target),1);
avgNotTarget = mean(~isnan(nonTarget),1);

for i = 1:rows
    for j = 1:cols
        if isnan(genes(i,j))
            %doReplace = (isnan(genes(i,:)));
            if isTarget(i) == 1
                genes(i,j) = avgTarget(j);
            elseif isTarget(i) == 0
                genes(i,j) = avgNotTarget(j);
            else
                print("Error: element does not belong to target or non-target.");
            end
        end
    end
end

nans_left = genes(isnan(genes));  % --> none left
geneNames = table2array(colInfo(:,1));

end
