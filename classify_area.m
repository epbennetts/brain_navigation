clear vars;
close all;
load("C:\Users\elobe\Documents\1_HONOURS\Gene_data_and_mat_scripts\GeneDataFile.mat", '-mat')

%looks like geneData is already Z-scored
%ZgeneData = zscore(geneData);

%instead:
%loop through last column of structInfo 
%make array of 1/0 depending if is in Isocortex or not
%Don't have a better way of doing this atm but should find one!
%This is also dependent on the dataset structure being as it is (Isocortex all in one place)
classes = zeros(1,213);
%how to extract the contents of this table to make it a string?
for i = 1:213
    area = table2array(structInfo(i,5));
    if strcmp(area,'Isocortex') == 1 %never happening
        classes(i) = 1;
    end
end


%class = (structInfo(:,5));
%tree = fitctree(geneData,class);


% %loop over areas
% for i = 1:5
%     gene_expr = zeros(1,19419);
%     
%     %loop over genes
%     for j = 1:19419
%         % how different is this gene expressed in Isocortex vs ~Isocortex? 
%         %make arrays of first rows
%             if isnan(ZgeneData(i,j)) 
%                 ZgeneData(i,j) = 0;
%             end
%             gene_expr(j) = ZgeneData(i,j);
%            % gene_entropy(j) = 
%             %else it is NaN
%             %in which case we just leave it at 0    
%             
%     end
%     figure()
%     plot(1:19419,gene_expr)
%     %should be plotting expression of one gene in all areas 
%     %instead of all genes in 1 area
% 
% end
