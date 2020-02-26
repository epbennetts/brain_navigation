clear vars;
%load new dataset (have to check which section)
load("C:\Users\elobe\Dropbox\1_HONOURS\code\data\AllenGeneDataset_19419_Ben.mat", '-mat')
genes = GeneExpData.combZ.energy; %for example
[rows, cols] = size(genes);
colLabels = geneInfo;
rowLabels = structInfo;

% %nans in each row
% nans_in_rows = zeros(1,rows);
% for i = 1:rows
% row = genes(i,:);
% [x, nans_row] = size(row(isnan(row)));
% nans_in_rows(i) = nans_row;
% end
% 
% %nans in each col
% nans_in_cols = zeros(1,cols);
% for i = 1:cols
% col = genes(:,i)';
% nancols = col(isnan(col));
% [x, nans_col] = size(nancols);
% nans_in_cols(i) = nans_col;
% end

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

%should probably also have a min amount of vals in target area (e.g.
%isocortex) --> later




 
% %deleting rows
% %while there are still more rows:
% %check if row has too many nans
% %if so, delete (next row will replace current one)
% %if not, go to next row
% i = 1;
% while i <= length(genes(:,1))
%     curr_row = genes(i,:);
%     if length(curr_row(isnan(curr_row)))/cols*100 > max_nans_rows
%         genes(i,:) = [];
%     else
%         i = i + 1;
%     end
% end
% rows = length(genes(:,1));
% 
% %deleting cols
% j = 1;
% while j <= length(genes(1,:))
%     curr_col = genes(:,j);
%     if length(curr_col(isnan(curr_col)))/rows*100 > max_nans_rows
%         genes(:,j) = [];
%     else
%         j = j + 1;
%     end
% end
% cols = length(genes(1,:));
% 
% %calculate average of each gene(col) and set remaining NaNs to that.
% for j = 1:cols
%     gene = genes(:,j);
%     valid_vals_in_col = gene(~isnan(gene));
%     avg = mean(valid_vals_in_col);
% 
%     for i = 1:rows
%         if isnan(genes(i,j))
%             genes(i,j) = avg;
%         end
%     end
% end
% 
% %nans_left = genes(isnan(genes));  % --> none left
