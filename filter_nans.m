clear vars;
%load new dataset (have to check which section)
load("C:\Users\elobe\Dropbox\1_HONOURS\code\data\AllenGeneDataset_19419_Ben.mat", '-mat')
genes = GeneExpData.combZ.energy; %for example
[rows, cols] = size(genes);

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
figure(2);
histogram(propNansCols,'Normalization', 'probability', 'BinLimits',[0,0.3]);



% %doRemove = (propNansCols>0.01);
% 
% 
% %see how many we would be deleting when we
% %delete anything with more than 1%, 5%, 10%, 15% nans
% perc_vals = [1,5,10,15,20,30,40,50];
% 
% rows_deleted = zeros(length(perc_vals));
% perc_rows_deleted = zeros(length(perc_vals));
% for i = 1:length(perc_vals)
%     [x,r_d] = size(nans_in_rows(nans_in_rows > perc_vals(i)));
%     p_r_d = r_d/rows*100;
%     rows_deleted(i) = r_d;
%     perc_rows_deleted(i) = p_r_d;
% end
% 
% cols_deleted = zeros(length(perc_vals));
% perc_cols_deleted = zeros(length(perc_vals));
% for i = 1:length(perc_vals)
%     [x,c_d] = size(nans_in_cols(nans_in_cols > perc_vals(i)));
%     p_c_d = c_d/cols*100;
%     cols_deleted(i) = c_d;
%     perc_cols_deleted(i) = p_c_d;
% end
% 
% %view: cols_deleted and perc_cols_deleted
% %decide which threshold to use: e.g. max 15% Nans in rows and max 10% NaNs
% %in cols
% %this means deleting 5% of genes and 34% of areas... Latter is maybe not ideal but good for now
% max_nans_rows = 15;
% max_nans_cols = 10;
% 
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
