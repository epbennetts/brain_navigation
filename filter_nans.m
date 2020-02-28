%TO DO 
%put cols into DT

clear vars;
%load new dataset (have to check which section)
load("C:\Users\elobe\Dropbox\1_HONOURS\code\data\AllenGeneDataset_19419_Ben.mat", '-mat')
genes = GeneExpData.combZ.energy; %for example
[rows, cols] = size(genes);
colLabels = geneInfo;
rowLabels = structInfo;

%calc how many nans in rows and cols
propNansRows = mean(isnan(genes),2);
propNansCols = mean(isnan(genes),1);

%histogram of what comes out of each of these commands.
figure(1);
histogram(propNansRows, 'Normalization', 'probability');
title('Proportion of nans in rows');
figure(2);
histogram(propNansCols,'Normalization', 'probability', 'BinLimits',[0,0.3]);
title('Proportion of nans in columns')

max_nans_rows = 0.1;
max_nans_cols = 0.1;

%cols
%this returns an array of size cols with 0/1 values
doRemove = (propNansCols > max_nans_cols);
%remove the arrays at indices with 1
genes(:,doRemove) = [];
colLabels(doRemove,:) = [];

%rows
doRemove = (propNansRows > max_nans_rows);
genes(doRemove,:) = [];
rowLabels(doRemove,:) = [];
areas = table2array(rowLabels(:,5));

[rows, cols] = size(genes);

%should probably also have a min amount of vals in target area (e.g.
%isocortex) --> later

%calc average of target and ~target elements
area = "Isocortex";
targetIndices = (strcmp(area, areas));
notTargetIndices = (targetIndices == 0);
target = genes(targetIndices, :);
notTarget = genes(notTargetIndices, :);

avgTarget = mean(~isnan(target),1);
avgNotTarget = mean(~isnan(notTarget),1);

%set remaining nans to average of their area (target or not target)
for i = 1:rows
    for j = 1:cols
        if isnan(genes(i,j))
            %doReplace = (isnan(genes(i,:)));
            if targetIndices(i) == 1
                genes(i,j) = avgTarget(j);
            elseif targetIndices(i) == 0
                genes(i,j) = avgNotTarget(j);
            else
                print("Error: element does not belong to target or not-target.")
            end
        end
    end
end

nans_left = genes(isnan(genes));  % --> none left
